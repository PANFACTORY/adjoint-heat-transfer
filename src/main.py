import matplotlib.pyplot as plt
import numpy as np

from fast_solver import forward

# from reference_solver import forward

lx, ly = 100.0, 100.0
nx, ny = 100, 100
dx, dy = lx / nx, ly / ny

u0, q0 = 0.0, 1e-1
nt = 100000
dt = 0.1

u_old, u_new = np.zeros((nx, ny)), np.zeros((nx, ny))
qx, qy = np.zeros((nx + 1, ny)), np.zeros((nx, ny + 1))
k = np.ones((nx, ny))

for t in range(nt):
    if t % 100 == 0:
        print(f"\rt={t}", end="")

    u_new = forward(u_old, qx, qy, k, nx, ny, dx, dy, dt, u0, q0)
    # u_new = forward(u_old, u_new, qx, qy, k, nx, ny, dx, dy, dt, u0, q0)
    u_old, u_new = u_new, u_old

plt.imshow(u_old.T, origin="lower", cmap="viridis")
plt.colorbar()
plt.show()
