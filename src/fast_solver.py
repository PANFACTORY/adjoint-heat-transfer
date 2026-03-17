import numpy as np


def forward(
    u_old: np.ndarray,
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
    k_facex = k[: nx - 1, :] * k[1:, :] / (k[: nx - 1, :] + k[1:, :])
    qx[1:nx, :] = -k_facex * (u_old[1:, :] - u_old[: nx - 1, :]) / dx

    # Compute qx on boundary
    qx[0, :] = -k[0, :] * (u_old[0, :] - u0) / (0.5 * dx)
    qx[nx, :] = -k[nx - 1, :] * (u0 - u_old[nx - 1, :]) / (0.5 * dx)

    # Compute qy in domain
    k_facey = k[:, : ny - 1] * k[:, 1:] / (k[:, : ny - 1] + k[:, 1:])
    qy[:, 1:ny] = -k_facey * (u_old[:, 1:] - u_old[:, : ny - 1]) / dy

    # Compute qy on boundary
    qy[:, 0] = q0
    qy[:, ny] = -k[:, ny - 1] * (u0 - u_old[:, ny - 1]) / (0.5 * dy)

    # Update cell-centered value
    return u_old - dt * ((qx[1:, :] - qx[:nx, :]) / dx + (qy[:, 1:] - qy[:, :ny]) / dy)
