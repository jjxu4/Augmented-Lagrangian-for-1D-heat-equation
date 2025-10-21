import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Parameters
alpha = 0.01      # diffusion coefficient
L = 1.0           # spatial domain length
T = 1.0           # total time
Nx = 60           # number of spatial steps
Nt = 800          # number of time steps
dx = L / Nx
dt = T / Nt

# Stability check for explicit scheme
if alpha * dt / dx**2 > 0.5:
    print("Warning: stability condition may be violated!")

# Discretization grids
x = np.linspace(0, L, Nx + 1)
t = np.linspace(0, T, Nt + 1)

# Initialize temperature field and source term
y = np.zeros((Nt + 1, Nx + 1))
f = np.zeros_like(y)

# Gaussian source term parameters
A = 20.0
x0, t0 = 0.5, 0.5
sigma_x, sigma_t = 0.05, 0.05

# Build source term f(x,t)
for n in range(Nt + 1):
    f[n, :] = A * np.exp(-((x - x0)**2) / (2 * sigma_x**2)) * \
              np.exp(-((t[n] - t0)**2) / (2 * sigma_t**2))

# Time stepping loop (explicit Euler)
for n in range(Nt):
    for i in range(1, Nx):
        y[n + 1, i] = y[n, i] + dt * (
            alpha * (y[n, i + 1] - 2 * y[n, i] + y[n, i - 1]) / dx**2 + f[n, i]
        )
    # Dirichlet boundary conditions
    y[n + 1, 0] = 0
    y[n + 1, -1] = 0

# Create mesh for plotting
X, Tt = np.meshgrid(x, t)

# 3D plots
fig = plt.figure(figsize=(12, 8))

# Temperature field
ax1 = fig.add_subplot(1, 2, 1, projection='3d')
ax1.plot_surface(X, Tt, y, cmap='viridis')
ax1.set_title("Temperature field y(x,t)")
ax1.set_xlabel("x")
ax1.set_ylabel("t")
ax1.set_zlabel("y")

# Source term
ax2 = fig.add_subplot(1, 2, 2, projection='3d')
ax2.plot_surface(X, Tt, f, cmap='plasma')
ax2.set_title("Gaussian source term f(x,t)")
ax2.set_xlabel("x")
ax2.set_ylabel("t")
ax2.set_zlabel("f")

plt.tight_layout()
plt.show()
