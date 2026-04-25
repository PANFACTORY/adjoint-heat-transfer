import matplotlib.pyplot as plt
import numpy as np

from fast_solver import forward, reverse

# from reference_solver import forward
# from reference_solver import reverse_gather as reverse

lx, ly = 100.0, 100.0
nx, ny, nt = 100, 100, 10000
dx, dy, dt = lx / nx, ly / ny, 0.1

u0, q0, qv = 0.0, 1e-1, 0.0
k_min, k_max = 0.01, 1

i, j = np.arange(nx)[:, None], np.arange(ny)[None, :]
gamma = np.where((np.abs(i - nx / 2) < nx / 4) & (j < ny / 2), 0.9, 0.1)
k = (k_max - k_min) * gamma + k_min

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

# Forward analysis
for t in range(1, nt):
    if t % 100 == 0:
        print(f"Forward analysis: t={t}", end="\r")

    u[t], qx, qy = forward(u[t - 1], u[t], qx, qy, k, nx, ny, dx, dy, dt, db, bu, bq)
    u[t] += qv * dt

# Reverse analysis
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

# Design sensitivity
df = (k_max - k_min) * ak

# Show result
fig, (ax_k, ax_u, ax_au, ax_ak) = plt.subplots(1, 4, figsize=(12, 4))

im_k = ax_k.imshow(k.T, origin="lower", cmap="bwr")
ax_k.set_title("Design variable")
fig.colorbar(im_k, ax=ax_k)
ax_k.axis("off")

im_u = ax_u.imshow(u[nt - 1].T, origin="lower", cmap="inferno")
ax_u.set_title("Temperature")
fig.colorbar(im_u, ax=ax_u)
ax_u.axis("off")

im_au = ax_au.imshow(au_tm1.T, origin="lower", cmap="inferno_r")
ax_au.set_title("Adjoint Temperature")
fig.colorbar(im_au, ax=ax_au)
ax_au.axis("off")

im_ak = ax_ak.imshow(df.T, origin="lower", cmap="jet")
ax_ak.set_title("Design sensitivity")
fig.colorbar(im_ak, ax=ax_ak)
ax_ak.axis("off")

plt.tight_layout()
plt.show()
