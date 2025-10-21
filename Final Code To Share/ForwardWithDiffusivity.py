import numpy as np
import matplotlib.pyplot as plt

# ============================================
# Step 1: Initialize stuff
# ============================================

# grid dimensions
nx, nt = 2**8, 2**8               # number of pixels in length and time
L, T = 1.0, 1.0                   # length and time
x = np.linspace(0, L, nx)
t = np.linspace(0, T, nt)
dx = L / (nx - 1)
dt = T / (nt - 1)
X, Tm = np.meshgrid(x, t, indexing="ij")

# diffusivity parameter 
alpha = 0.01   

# initial condition
y0 = np.zeros(nx)

# ============================================
# Step 1a: Build spatial operator
# ============================================

def build_laplacian(nx, dx, alpha=1.0):
    N = nx - 2
    main = -2.0 * np.ones(N)
    off = np.ones(N - 1)
    Dxx = alpha * (np.diag(main) + np.diag(off, 1) + np.diag(off, -1)) / (dx ** 2)
    return Dxx

Dxx = build_laplacian(nx, dx, alpha)
I = np.eye(nx - 2)

# ============================================
# Step 2: Define forward solver
# ============================================

def forward_solve(u, y0, Dxx, dt, alpha=1.0):
    """Forward solve: y_t - alpha * y_xx = u."""
    y = np.zeros_like(u)
    y[:, 0] = y0
    M = I - dt * Dxx
    for n in range(nt - 1):
        rhs = y[1:-1, n] + dt * u[1:-1, n + 1]
        y[1:-1, n + 1] = np.linalg.solve(M, rhs)
    return y

# ============================================
# Step 3: Define a sample control u(x,t)
# ============================================

# Space–time Gaussian 
x0 = 0.5       # center in space
t0 = 0.5       # center in time
sigma_x = 0.05 # spatial width
sigma_t = 0.05 # temporal width
amplitude = 100

u = amplitude * np.exp(-((X - x0)**2) / (2 * sigma_x**2)) * np.exp(-((Tm - t0)**2) / (2 * sigma_t**2))

# ============================================
# Step 4: Solve forward system
# ============================================

y = forward_solve(u, y0, Dxx, dt, alpha)

# ============================================
# Step 4b: Save solution
# ============================================
np.save("y_desired_w_diff.npy", y)
print("Saved solution. OuterloopWithDiffusivity.py uses this solution")

# ============================================
# Step 5: Plot results (u and y)
# ============================================

fig = plt.figure(figsize=(12, 8))

# plot y(x,t)
ax1 = fig.add_subplot(121, projection="3d")
ax1.plot_surface(X, Tm, y, cmap="viridis", edgecolor="none", alpha=0.9)
ax1.set_title("State y(x,t)")
ax1.set_xlabel("x")
ax1.set_ylabel("t")
ax1.set_zlabel("y")

# plot u(x,t)
ax2 = fig.add_subplot(122, projection="3d")
ax2.plot_surface(X, Tm, u, cmap="plasma", edgecolor="none", alpha=0.9)
ax2.set_title("Control u(x,t)")
ax2.set_xlabel("x")
ax2.set_ylabel("t")
ax2.set_zlabel("u")

plt.suptitle(f"Heat Equation Solution (alpha = {alpha})")
plt.tight_layout()
plt.show()
