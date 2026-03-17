import numpy as np


def forward(
    u_old: np.ndarray,
    u_new: np.ndarray,
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
) -> np.ndarray:
    # Compute x flux
    for j in range(ny):
        # inner faces
        for i in range(1, nx):
            k_face = k[i - 1, j] * k[i, j] / (k[i - 1, j] + k[i, j])
            qx[i, j] = -k_face * (u_old[i, j] - u_old[i - 1, j]) / dx

        # boundaries
        qx[0, j] = -k[0, j] * (u_old[0, j] - u0) / (0.5 * dx)
        qx[nx, j] = -k[nx - 1, j] * (u0 - u_old[nx - 1, j]) / (0.5 * dx)

    # Compute y flux
    for i in range(nx):
        # inner faces
        for j in range(1, ny):
            k_face = k[i - 1, j] * k[i, j] / (k[i - 1, j] + k[i, j])
            qy[i, j] = -k_face * (u_old[i, j] - u_old[i, j - 1]) / dy

        # boundaries
        qy[i, 0] = q0
        qy[i, ny] = -k[i, ny - 1] * (u0 - u_old[i, ny - 1]) / (0.5 * dy)

    # Update cell-centered value
    for i in range(nx):
        for j in range(ny):
            u_new[i, j] = u_old[i, j] - dt * (
                (qx[i + 1, j] - qx[i, j]) / dx + (qy[i, j + 1] - qy[i, j]) / dy
            )

    return u_new
