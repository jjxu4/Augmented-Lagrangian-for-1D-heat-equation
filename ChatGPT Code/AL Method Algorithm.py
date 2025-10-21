import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ==============================================================
# Forward System
# ==============================================================

alpha = 0.01
L = 1.0
T = 1.0
Nx = 60
Nt = 800
dx = L / Nx
dt = T / Nt

x = np.linspace(0, L, Nx + 1)
t = np.linspace(0, T, Nt + 1)
X, TT = np.meshgrid(x, t)

# --- True Gaussian source ---
A = 20.0
x0, t0 = 0.5, 0.3
sigma_x, sigma_t = 0.08, 0.05
f_true = np.zeros((Nt + 1, Nx + 1))
for n in range(Nt + 1):
    f_true[n, :] = A * np.exp(-((x - x0)**2)/(2*sigma_x**2)) * np.exp(-((t[n]-t0)**2)/(2*sigma_t**2))

# --- Forward solve to get true solution ---
y_true = np.zeros_like(f_true)
for n in range(Nt):
    for i in range(1, Nx):
        y_true[n + 1, i] = y_true[n, i] + dt * (
            alpha * (y_true[n, i + 1] - 2 * y_true[n, i] + y_true[n, i - 1]) / dx**2 + f_true[n, i]
        )

y_desired = y_true.copy()

# ==============================================================
# Augmented Lagrange Setup
# ==============================================================

y_bar = np.ones_like(y_desired)  # inactive state bound

gamma = 1e-3
rho0 = 10.0
rho_growth = 1.2
inner_iters = 25
lr = 0.08

def forward(u):
    y = np.zeros_like(u)
    for n in range(Nt):
        y[n+1,1:-1] = y[n,1:-1] + dt*(alpha*(y[n,2:]-2*y[n,1:-1]+y[n,:-2])/dx**2 + u[n,1:-1])
    return y

def adjoint(y, mu, rho):
    p = np.zeros_like(y)
    for n in range(Nt-1, -1, -1):
        viol = y[n,:] - y_bar[n,:] + mu[n,:]/rho
        rhs = (y[n,:] - y_desired[n,:]) + rho * np.maximum(0.0, viol)
        lap_p_next = np.zeros_like(p[n+1,:])
        lap_p_next[1:-1] = (p[n+1,2:] - 2*p[n+1,1:-1] + p[n+1,:-2]) / dx**2
        p[n,:] = p[n+1,:] + dt*(alpha * lap_p_next + rhs)
        p[n,0] = 0.0; p[n,-1] = 0.0
    return p

def al_objective(y, u, mu, rho):
    z = y - y_bar + mu / rho
    z_pos = np.maximum(0.0, z)
    penalty = 0.5 * rho * np.sum(z_pos**2) - (0.5/rho) * np.sum(mu**2)
    tracking = 0.5 * np.sum((y - y_desired)**2)
    reg = 0.5 * gamma * np.sum(u**2)
    return (tracking + reg + penalty) * dx * dt

# ==============================================================
# Initialize
# ==============================================================

u = np.zeros_like(f_true)   # start from zero control
mu = np.zeros_like(u)
rho = rho0

# ==============================================================
# AL Loop
# ==============================================================

snapshots = {}
prev_violation_norm = None
k = 0

print("Running indefinitely. Press Ctrl + C to stop.\n")

try:
    while True:
        k += 1

        # ---- Inner minimization (gradient descent) ----
        for _ in range(inner_iters):
            y = forward(u)
            p = adjoint(y, mu, rho)
            grad_u = gamma * u + p
            u -= lr * grad_u
            u[:,0] = 0.0; u[:,-1] = 0.0

        # ---- Outer updates ----
        y = forward(u)
        violation = y - y_bar
        mu = np.maximum(0.0, mu + rho * violation)

        viol_norm = np.sqrt(np.sum(np.maximum(0.0, violation)**2) * dx * dt)
        if prev_violation_norm is not None and viol_norm > 0.9 * prev_violation_norm and viol_norm > 1e-6:
            rho *= rho_growth
        prev_violation_norm = viol_norm

        # ---- Plot every 5 iterations ----
        if k % 1 == 0:
            print(f"Iteration {k}: rho={rho:.2f}, ||viol||={viol_norm:.2e}")
            ds_t, ds_x = 4, 2
            Xd, TTd = X[::ds_t, ::ds_x], TT[::ds_t, ::ds_x]
            ykd, ukd = y[::ds_t, ::ds_x], u[::ds_t, ::ds_x]
            y_bard = y_bar[::ds_t, ::ds_x]

            fig = plt.figure(figsize=(12,8))
            ax1 = fig.add_subplot(1,2,1, projection='3d')
            ax1.plot_surface(Xd, TTd, ykd, cmap='viridis')
            ax1.plot_wireframe(Xd, TTd, y_bard, rstride=6, cstride=6, alpha=0.3)
            ax1.set_title(f"Iter {k}: State y(x,t)")
            ax1.set_xlabel("x"); ax1.set_ylabel("t"); ax1.set_zlabel("y")

            ax2 = fig.add_subplot(1,2,2, projection='3d')
            ax2.plot_surface(Xd, TTd, ukd, cmap='plasma')
            ax2.set_title(f"Iter {k}: Control u(x,t)")
            ax2.set_xlabel("x"); ax2.set_ylabel("t"); ax2.set_zlabel("u")

            plt.tight_layout()
            plt.show()

except KeyboardInterrupt:
    print("\nOptimization stopped by user.\n")

    # ==============================================================
    # Comparisons
    # ==============================================================

    ds_t, ds_x = 4, 2
    Xd, TTd = X[::ds_t, ::ds_x], TT[::ds_t, ::ds_x]

    # -- Source comparison
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(121, projection='3d')
    ax.plot_surface(Xd, TTd, f_true[::ds_t,::ds_x], cmap='viridis', alpha=0.6)
    ax.plot_wireframe(Xd, TTd, u[::ds_t,::ds_x], color='r', alpha=0.8)
    ax.set_title("True source (surface) vs Recovered source (wireframe)")
    ax.set_xlabel("x"); ax.set_ylabel("t"); ax.set_zlabel("f/u")

    # -- Solution comparison
    ax2 = fig.add_subplot(122, projection='3d')
    ax2.plot_surface(Xd, TTd, y_desired[::ds_t,::ds_x], cmap='viridis', alpha=0.6)
    ax2.plot_wireframe(Xd, TTd, y[::ds_t,::ds_x], color='r', alpha=0.8)
    ax2.set_title("True solution (surface) vs Recovered solution (wireframe)")
    ax2.set_xlabel("x"); ax2.set_ylabel("t"); ax2.set_zlabel("y")

    plt.tight_layout()
    plt.show()
