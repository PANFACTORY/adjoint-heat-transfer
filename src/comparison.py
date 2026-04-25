import matplotlib.pyplot as plt
import numpy as np

from fast_solver import forward, reverse

lx, ly = 100.0, 100.0
nx, ny, nt = 100, 100, 10000
dx, dy, dt = lx / nx, ly / ny, 0.1

u0, q0, qv = 0.0, 1e-1, 0.0
k_min, k_max = 0.01, 1

dgamma = 1e-5

u = np.zeros((nt, nx, ny))
qx, qy = np.zeros((nx + 1, ny)), np.zeros((nx, ny + 1))
au_tm1, au_t, ak = np.zeros((nx, ny)), np.zeros((nx, ny)), np.zeros((nx, ny))
aqx, aqy = np.zeros((nx + 1, ny)), np.zeros((nx, ny + 1))

i, j = np.arange(nx), np.arange(ny)
is_heat_source = np.abs(i - 0.5 * nx) < 0.1 * nx
db = (
    np.full(ny, True),
    np.full(ny, True),
    np.full(nx, False),
    np.full(nx, True),
)
bu = (np.full(ny, u0), np.full(ny, u0), np.full(nx, u0), np.full(nx, u0))
bq = (
    np.full(ny, 0.0),
    np.full(ny, 0.0),
    np.where(is_heat_source, q0, 0.0),
    np.full(nx, 0.0),
)

comp_points = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90]

# FDM
for j0 in comp_points:
    print(f"j={j0}", end="\r")

    ak[nx // 2, j0] = 0

    for d in [-1, 1]:
        i, j = np.arange(nx)[:, None], np.arange(ny)[None, :]
        gamma = np.where((np.abs(i - nx / 2) < nx / 4) & (j < ny / 2), 0.9, 0.1)
        gamma[nx // 2, j0] += d * dgamma
        k = (k_max - k_min) * gamma + k_min

        # Forward analysis
        u[0, :, :] = 0
        for t in range(1, nt):
            u[t], qx, qy = forward(
                u[t - 1], u[t], qx, qy, k, nx, ny, dx, dy, dt, db, bu, bq
            )
            u[t] += qv * dt

            ak[nx // 2, j0] += d * np.sum(u[t, is_heat_source, 0]) / lx / ly / nt

    ak[nx // 2, j0] /= 2 * dgamma

df_fdm = (k_max - k_min) * ak

# Adjoint
u[0, :, :] = 0
i, j = np.arange(nx)[:, None], np.arange(ny)[None, :]
gamma = np.where((np.abs(i - nx / 2) < nx / 4) & (j < ny / 2), 0.9, 0.1)
k = (k_max - k_min) * gamma + k_min
ak[:, :] = 0
for t in range(1, nt):
    if t % 100 == 0:
        print(f"Forward analysis: t={t}", end="\r")

    u[t], qx, qy = forward(u[t - 1], u[t], qx, qy, k, nx, ny, dx, dy, dt, db, bu, bq)
    u[t] += qv * dt

for t in range(nt - 1, 0, -1):
    if t % 100 == 0:
        print(f"Reverse analysis: t={t}", end="\r")

    au_tm1.fill(0)
    aqx.fill(0)
    aqy.fill(0)

    au_t[is_heat_source, 0] += 1 / lx / ly / nt

    au_tm1, aqx, aqy, ak = reverse(
        au_tm1, au_t, aqx, aqy, ak, u[t - 1], k, nx, ny, dx, dy, dt, db, bu
    )
    au_tm1, au_t = au_t, au_tm1

df_adj = (k_max - k_min) * ak

# Plot result
plt.figure()
plt.plot(comp_points, df_adj[nx // 2, comp_points], "o", label="Adjoint")
plt.plot(comp_points, df_fdm[nx // 2, comp_points], "x--", label="FDM")

plt.xlabel("y")
plt.ylabel("Design sensitivity")
plt.legend()

plt.show()
