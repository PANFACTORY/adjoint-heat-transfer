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
) -> np.ndarray:
    # Compute qx in domain
    k_facex = 2 * k[: nx - 1, :] * k[1:, :] / (k[: nx - 1, :] + k[1:, :])
    qx[1:nx, :] = -k_facex * (u_tm1[1:, :] - u_tm1[: nx - 1, :]) / dx

    # Compute qx on boundary
    qx[0, :] = -k[0, :] * (u_tm1[0, :] - u0) / (0.5 * dx)
    qx[nx, :] = -k[nx - 1, :] * (u0 - u_tm1[nx - 1, :]) / (0.5 * dx)

    # Compute qy in domain
    k_facey = 2 * k[:, : ny - 1] * k[:, 1:] / (k[:, : ny - 1] + k[:, 1:])
    qy[:, 1:ny] = -k_facey * (u_tm1[:, 1:] - u_tm1[:, : ny - 1]) / dy

    # Compute qy on boundary
    qy[:, 0] = q0
    qy[:, ny] = -k[:, ny - 1] * (u0 - u_tm1[:, ny - 1]) / (0.5 * dy)

    # Update cell-centered value
    u_t = u_tm1 - dt * ((qx[1:, :] - qx[:nx, :]) / dx + (qy[:, 1:] - qy[:, :ny]) / dy)
    return u_t, qx, qy


def reverse(
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
    au_tm1 += au_t

    aqx[1 : nx + 1, :] += -dt / dx * au_t
    aqx[:nx, :] += dt / dx * au_t

    aqy[:, 1 : ny + 1] += -dt / dx * au_t
    aqy[:, :ny] += dt / dy * au_t

    # Compute y flux
    au_tm1[:, ny - 1] += k[:, ny - 1] / (0.5 * dy) * aqy[:, ny]
    ak[:, ny - 1] += -(u0 - u_tm1[:, ny - 1]) / (0.5 * dy) * aqy[:, ny]

    k_face = 2 * k[:, : ny - 1] * k[:, 1:ny] / (k[:, : ny - 1] + k[:, 1:ny])
    au_tm1[:, 1:ny] += -k_face / dy * aqy[:, 1:ny]
    au_tm1[:, : ny - 1] += k_face / dy * aqy[:, 1:ny]

    ak_face = -(u_tm1[:, 1:ny] - u_tm1[:, : ny - 1]) / dy * aqy[:, 1:ny]
    ak[:, 1:ny] += 2 * (k[:, : ny - 1] / (k[:, : ny - 1] + k[:, 1:ny])) ** 2 * ak_face
    ak[:, : ny - 1] += 2 * (k[:, 1:ny] / (k[:, : ny - 1] + k[:, 1:ny])) ** 2 * ak_face

    # Compute x flux
    au_tm1[0, :] += -k[0, :] / (0.5 * dx) * aqx[0, :]
    au_tm1[nx - 1, :] += k[nx - 1, :] / (0.5 * dx) * aqx[nx, :]

    ak[0, :] += -(u_tm1[0, :] - u0) / (0.5 * dx) * aqx[0, :]
    ak[nx - 1, :] += -(u0 - u_tm1[nx - 1, :]) / (0.5 * dx) * aqx[nx, :]

    k_face = 2 * k[: nx - 1, :] * k[1:nx, :] / (k[: nx - 1, :] + k[1:nx, :])
    au_tm1[1:nx, :] += -k_face / dx * aqx[1:nx, :]
    au_tm1[: nx - 1, :] += k_face / dx * aqx[1:nx, :]

    ak_face = -(u_tm1[1:nx, :] - u_tm1[: nx - 1, :]) / dx * aqx[1:nx, :]
    ak[: nx - 1, :] += 2 * (k[1:nx, :] / (k[: nx - 1, :] + k[1:nx, :])) ** 2 * ak_face
    ak[1:nx, :] += 2 * (k[: nx - 1, :] / (k[: nx - 1, :] + k[1:nx, :])) ** 2 * ak_face

    return au_tm1, aqx, aqy, ak
