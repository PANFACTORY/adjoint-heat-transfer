import numpy as np


def forward(
    u_tm1: np.ndarray,
    u_t: np.ndarray,
    qx: np.ndarray,
    qy: np.ndarray,
    k: np.ndarray,
    nx: int,
    ny: int,
    dx: float,
    dy: float,
    dt: float,
    u0: float,
    q0: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    # Compute x flux
    for j in range(ny):
        # inner faces
        for i in range(1, nx):
            k_face = 2 * k[i - 1, j] * k[i, j] / (k[i - 1, j] + k[i, j])
            qx[i, j] = -k_face * (u_tm1[i, j] - u_tm1[i - 1, j]) / dx

        # boundaries
        qx[0, j] = -k[0, j] * (u_tm1[0, j] - u0) / (0.5 * dx)
        qx[nx, j] = -k[nx - 1, j] * (u0 - u_tm1[nx - 1, j]) / (0.5 * dx)

    # Compute y flux
    for i in range(nx):
        # inner faces
        for j in range(1, ny):
            k_face = 2 * k[i, j - 1] * k[i, j] / (k[i, j - 1] + k[i, j])
            qy[i, j] = -k_face * (u_tm1[i, j] - u_tm1[i, j - 1]) / dy

        # boundaries
        qy[i, 0] = q0
        qy[i, ny] = -k[i, ny - 1] * (u0 - u_tm1[i, ny - 1]) / (0.5 * dy)

    # Update cell-centered value
    for i in range(nx):
        for j in range(ny):
            u_t[i, j] = u_tm1[i, j] - dt * (
                (qx[i + 1, j] - qx[i, j]) / dx + (qy[i, j + 1] - qy[i, j]) / dy
            )

    return u_t, qx, qy


def reverse_scatter(
    au_tm1: np.ndarray,
    au_t: np.ndarray,
    aqx: np.ndarray,
    aqy: np.ndarray,
    ak: np.ndarray,
    u_tm1: np.ndarray,
    k: np.ndarray,
    nx: int,
    ny: int,
    dx: float,
    dy: float,
    dt: float,
    u0: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # Update cell-centered value
    for i in range(nx):
        for j in range(ny):
            # u_tp1[i, j] = u_t[i, j] - dt * (
            #     (qx[i + 1, j] - qx[i, j]) / dx + (qy[i, j + 1] - qy[i, j]) / dy
            # )
            au_tm1[i, j] += au_t[i, j]
            aqx[i + 1, j] += -dt / dx * au_t[i, j]
            aqx[i, j] += dt / dx * au_t[i, j]
            aqy[i, j + 1] += -dt / dy * au_t[i, j]
            aqy[i, j] += dt / dy * au_t[i, j]

    # Compute y flux
    for i in range(nx):
        # boundaries
        # qy[i, 0] = q0
        # qy[i, ny] = -k[i, ny - 1] * (u0 - u_t[i, ny - 1]) / (0.5 * dy)
        au_tm1[i, ny - 1] += k[i, ny - 1] / (0.5 * dy) * aqy[i, ny]

        ak[i, ny - 1] += -(u0 - u_tm1[i, ny - 1]) / (0.5 * dy) * aqy[i, ny]

        # inner faces
        for j in range(1, ny):
            k_face = 2 * k[i, j - 1] * k[i, j] / (k[i, j - 1] + k[i, j])
            # qy[i, j] = -k_face * (u_t[i, j] - u_t[i, j - 1]) / dy
            au_tm1[i, j] += -k_face / dy * aqy[i, j]
            au_tm1[i, j - 1] += k_face / dy * aqy[i, j]

            ak_face = -(u_tm1[i, j] - u_tm1[i, j - 1]) / dy * aqy[i, j]
            ak[i, j - 1] += 2 * (k[i, j] / (k[i, j - 1] + k[i, j])) ** 2 * ak_face
            ak[i, j] += 2 * (k[i, j - 1] / (k[i, j - 1] + k[i, j])) ** 2 * ak_face

    # Compute x flux
    for j in range(ny):
        # boundaries
        # qx[0, j] = -k[0, j] * (u_t[0, j] - u0) / (0.5 * dx)
        # qx[nx, j] = -k[nx - 1, j] * (u0 - u_t[nx - 1, j]) / (0.5 * dx)
        au_tm1[0, j] += -k[0, j] / (0.5 * dx) * aqx[0, j]
        au_tm1[nx - 1, j] += k[nx - 1, j] / (0.5 * dx) * aqx[nx, j]

        ak[0, j] += -(u_tm1[0, j] - u0) / (0.5 * dx) * aqx[0, j]
        ak[nx - 1, j] += -(u0 - u_tm1[nx - 1, j]) / (0.5 * dx) * aqx[nx, j]

        # inner faces
        for i in range(1, nx):
            k_face = 2 * k[i - 1, j] * k[i, j] / (k[i - 1, j] + k[i, j])
            # qx[i, j] = -k_face * (u_t[i, j] - u_t[i - 1, j]) / dx
            au_tm1[i, j] += -k_face / dx * aqx[i, j]
            au_tm1[i - 1, j] += k_face / dx * aqx[i, j]

            ak_face = -(u_tm1[i, j] - u_tm1[i - 1, j]) / dx * aqx[i, j]
            ak[i - 1, j] += 2 * (k[i, j] / (k[i - 1, j] + k[i, j])) ** 2 * ak_face
            ak[i, j] += 2 * (k[i - 1, j] / (k[i - 1, j] + k[i, j])) ** 2 * ak_face

    return au_tm1, aqx, aqy, ak


def reverse_gather(
    au_tm1: np.ndarray,
    au_t: np.ndarray,
    aqx: np.ndarray,
    aqy: np.ndarray,
    ak: np.ndarray,
    u_tm1: np.ndarray,
    k: np.ndarray,
    nx: int,
    ny: int,
    dx: float,
    dy: float,
    dt: float,
    u0: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    # Update cell-centered value
    for i in range(nx):
        for j in range(ny):
            au_tm1[i, j] += au_t[i, j]

    for i in range(nx + 1):
        for j in range(ny):
            if i > 0:
                aqx[i, j] += -dt / dx * au_t[i - 1, j]
            if i < nx:
                aqx[i, j] += dt / dx * au_t[i, j]

    for i in range(nx):
        for j in range(ny + 1):
            if j > 0:
                aqy[i, j] += -dt / dy * au_t[i, j - 1]
            if j < ny:
                aqy[i, j] += dt / dy * au_t[i, j]

    # Compute y flux on
    for i in range(nx):
        # boundaries
        au_tm1[i, ny - 1] += k[i, ny - 1] / (0.5 * dy) * aqy[i, ny]

        ak[i, ny - 1] += -(u0 - u_tm1[i, ny - 1]) / (0.5 * dy) * aqy[i, ny]

        # inner faces
        for j in range(ny):
            if j > 0:
                k_face = 2 * k[i, j - 1] * k[i, j] / (k[i, j - 1] + k[i, j])
                au_tm1[i, j] += -k_face / dy * aqy[i, j]
            if j < ny - 1:
                k_face = 2 * k[i, j] * k[i, j + 1] / (k[i, j] + k[i, j + 1])
                au_tm1[i, j] += k_face / dy * aqy[i, j + 1]

            if j > 0:
                ak_face = -(u_tm1[i, j] - u_tm1[i, j - 1]) / dy * aqy[i, j]
                ak[i, j] += 2 * (k[i, j - 1] / (k[i, j - 1] + k[i, j])) ** 2 * ak_face
            if j < ny - 1:
                ak_face = -(u_tm1[i, j + 1] - u_tm1[i, j]) / dy * aqy[i, j + 1]
                ak[i, j] += 2 * (k[i, j + 1] / (k[i, j] + k[i, j + 1])) ** 2 * ak_face

    # Compute x flux
    for j in range(ny):
        # boundaries
        au_tm1[0, j] += -k[0, j] / (0.5 * dx) * aqx[0, j]
        au_tm1[nx - 1, j] += k[nx - 1, j] / (0.5 * dx) * aqx[nx, j]

        ak[0, j] += -(u_tm1[0, j] - u0) / (0.5 * dx) * aqx[0, j]
        ak[nx - 1, j] += -(u0 - u_tm1[nx - 1, j]) / (0.5 * dx) * aqx[nx, j]

        # inner faces
        for i in range(nx):
            if i > 0:
                k_face = 2 * k[i - 1, j] * k[i, j] / (k[i - 1, j] + k[i, j])
                au_tm1[i, j] += -k_face / dx * aqx[i, j]
            if i < nx - 1:
                k_face = 2 * k[i, j] * k[i + 1, j] / (k[i, j] + k[i + 1, j])
                au_tm1[i, j] += k_face / dx * aqx[i + 1, j]

            if i > 0:
                ak_face = -(u_tm1[i, j] - u_tm1[i - 1, j]) / dx * aqx[i, j]
                ak[i, j] += 2 * (k[i - 1, j] / (k[i - 1, j] + k[i, j])) ** 2 * ak_face
            if i < nx - 1:
                ak_face = -(u_tm1[i + 1, j] - u_tm1[i, j]) / dx * aqx[i + 1, j]
                ak[i, j] += 2 * (k[i + 1, j] / (k[i, j] + k[i + 1, j])) ** 2 * ak_face

    return au_tm1, aqx, aqy, ak
