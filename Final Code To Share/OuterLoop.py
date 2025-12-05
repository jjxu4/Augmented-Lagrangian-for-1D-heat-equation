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

# problem data
y0 = np.zeros(nx)                 # initial guess
print("Setting desired temperature to y_desired.npy")
y_d = np.load("y_desired.npy")

# ============================================
# Define Constraint Here
# ============================================
# Centered in the middle of the domain (x=0.5, t=0.5)
x0, t0 = 0.5, 0.5

# Widths (broad Gaussian)
sigma_x = 0.2
sigma_t = 0.3

# Peak height (somewhat below the true unconstrained maximum)
amplitude = 1.2

# Construct space–time Gaussian
y_C = amplitude * np.exp(-((X - x0)**2) / (2 * sigma_x**2)) \
                 * np.exp(-((Tm - t0)**2) / (2 * sigma_t**2))
# ============================================

beta = 0                          # beta thing

# penalty and guess initial multiplier (mu)
mu = np.ones((nx, nt))
rho = 1000.0

# ============================================
# Step 1a: Build spatial operator (Dirichlet Laplacian)
# ============================================

def build_laplacian(nx, dx):
    N = nx - 2
    main = -2.0 * np.ones(N)
    off = np.ones(N - 1)
    Dxx = (np.diag(main) + np.diag(off, 1) + np.diag(off, -1)) / (dx ** 2)
    return Dxx

Dxx = build_laplacian(nx, dx)
I = np.eye(nx - 2)

# ============================================
# Step 1b: Define Augmented Lagrangian
# ============================================

def augmented_lagrangian(y, u, mu, rho, y_d, y_C, beta):
    diff = 0.5 * np.sum((y - y_d)**2)
    ctrl = 0.5 * beta * np.sum(u**2)
    z = y - y_C
    rock = (1.0 / (2.0 * rho)) * np.sum(np.maximum(0, mu + rho * z)**2 - mu**2)
    return diff + ctrl + rock

# ============================================
# Step 2: PDE Solvers
# ============================================

def forward_solve(u, y0, Dxx, dt):
    y = np.zeros_like(u)
    y[:, 0] = y0
    M = I - dt * Dxx
    for n in range(nt - 1):
        rhs = y[1:-1, n] + dt * u[1:-1, n + 1]
        y[1:-1, n + 1] = np.linalg.solve(M, rhs)
    return y

def adjoint_solve(y, y_d, mu_hat, Dxx, dt):
    p = np.zeros_like(y)
    M = I - dt * Dxx
    for n in range(nt - 2, -1, -1):
        rhs_field = (y[:, n] - y_d[:, n]) + mu_hat[:, n]
        rhs = p[1:-1, n + 1] + dt * rhs_field[1:-1]
        p[1:-1, n] = np.linalg.solve(M, rhs)
    return p

# ============================================
# Step 2b: Helpers
# ============================================

def compute_mu_hat(mu, rho, y, yC):
    return np.maximum(0.0, mu + rho * (y - yC))

def update_mu(mu, y, y_C, rho):
    return np.maximum(0, mu + rho * (y - y_C))

# ============================================
# Step 2c: Inner Loop (Minimization of L_A)
# ============================================

def minimize_L_A(mu, rho, y_d, y_C, beta, max_iter=50000, alpha=0.1, eps=0.001, u_hotguess=np.zeros((nx,nt))):
    """Gradient descent to minimize L_A."""
    u = u_hotguess
    y = forward_solve(u, y0, Dxx, dt)
    J_prev = augmented_lagrangian(y, u, mu, rho, y_d, y_C, beta)
    passed = False

    for k in range(max_iter):
        mu_hat = compute_mu_hat(mu, rho, y, y_C)
        p = adjoint_solve(y, y_d, mu_hat, Dxx, dt)
        g = beta * u + p
        u -= alpha * g
        y = forward_solve(u, y0, Dxx, dt)
        J_curr = augmented_lagrangian(y, u, mu, rho, y_d, y_C, beta)
        
        # Progress
        diff = abs(J_curr - J_prev)
        print(f"Iteration {k+1}/{max_iter}  |  Difference: {diff:.6e}", end="\r", flush=True)
        
        if diff < eps:
            print("\nStopping criterion achieved")
            passed = True
            break
        J_prev = J_curr
    if passed == False:
        print(f"\nStopping criterion wasn't achieved. {max_iter} iterations reached first")
    return y, u


# ============================================
# Step 3: Outer Loop (runs indefinitely, plots each iteration)
# ============================================

plt.ion()
fig = plt.figure(figsize=(12, 8))
ax1 = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(122, projection='3d')

try:
    max_iterations = 50
    iteration = 0
    u_new=np.zeros((nx,nt))
    while True:
        if iteration == max_iterations:
            break
        print("")
        print(f"=== Outer loop iteration {iteration} ===")
        print("Minimizing Augmented Lagrangian")
        y_new, u_new = minimize_L_A(mu, rho, y_d, y_C, beta, u_hotguess=u_new)
        
        # Print violation with state constraint
        violation = np.maximum(0, y_new - y_C)
        L2_viol = np.sum(violation)
        max_viol = np.max(violation)
        print(f"||viol||_L2 = {L2_viol:.6e}, ||viol||_max = {max_viol:.6e}")
        
        print("Updating multipliers")
        mu = update_mu(mu, y_new, y_C, rho)
        # r_K = calculate_R_K()
        # if 

        # --- Plot y(x,t), y_C, and u(x,t) ---
        fig = plt.figure(figsize=(12, 8))
        ax1 = fig.add_subplot(121, projection='3d')
        ax2 = fig.add_subplot(122, projection='3d')

        # State plot
        ax1.plot_wireframe(X, Tm, y_C, color='skyblue', alpha=0.5)
        ax1.plot_surface(X, Tm, y_new, cmap='viridis', edgecolor='none', alpha=0.9)
        ax1.set_title(f"Iter {iteration}: State y(x,t)")
        ax1.set_xlabel("x")
        ax1.set_ylabel("t")
        ax1.set_zlabel("y")

        # Control plot
        ax2.plot_surface(X, Tm, u_new, cmap='plasma', edgecolor='none', alpha=0.9)
        ax2.set_title(f"Iter {iteration}: Control u(x,t)")
        ax2.set_xlabel("x")
        ax2.set_ylabel("t")
        ax2.set_zlabel("u")
        
        # # Plot where control is violated
        # mask = np.abs(y_new - y_C) < 0.1
        # ax1.scatter(
        #     X[mask], Tm[mask], y_new[mask],
        #     color='red', s=5, alpha=0.8, label='Constraint active region'
        # )

        plt.suptitle(f"Iteration {iteration}")
        plt.tight_layout()
        plt.show(block=False)
        plt.pause(0.5)
        plt.close(fig)


        iteration += 1
    
# ============================================
# Step 4: Comparison plot after manually stopping
# ============================================
except KeyboardInterrupt:
    print("")
    print("\nInterrupted. Plotting final results.")
    fig_final = plt.figure(figsize=(12, 8))

    # --- Left: Control (u) ---
    ax1 = fig_final.add_subplot(121, projection='3d')
    u_desired = np.load("u_desired.npy")
    ax1.plot_surface(X, Tm, u_desired, cmap='viridis', alpha=0.8)
    ax1.plot_wireframe(X, Tm, u_new, color='r', alpha=0.6)
    ax1.set_title("True source (surface) vs Recovered source (wireframe)")
    ax1.set_xlabel("x")
    ax1.set_ylabel("t")
    ax1.set_zlabel("f/u")

    # --- Right: Solution (y) ---
    ax2 = fig_final.add_subplot(122, projection='3d')
    y_desired = np.load("y_desired.npy")
    ax2.plot_surface(X, Tm, y_desired, cmap='viridis', alpha=0.8)
    ax2.plot_wireframe(X, Tm, y_new, color='r', alpha=0.6)
    ax2.set_title("True solution (surface) vs Recovered solution (wireframe)")
    ax2.set_xlabel("x")
    ax2.set_ylabel("t")
    ax2.set_zlabel("y")

    plt.tight_layout()
    plt.show()

# ============================================
# Step 4: Comparison plot after stopping
# ============================================ 
print("")
print("Plotting final results.")
fig_final = plt.figure(figsize=(12, 8))

# --- Left: Control (u) ---
ax1 = fig_final.add_subplot(121, projection='3d')
u_desired = np.load("u_desired.npy")
ax1.plot_surface(X, Tm, u_desired, cmap='viridis', alpha=0.8)
ax1.plot_wireframe(X, Tm, u_new, color='r', alpha=0.6)
ax1.set_title("True source (surface) vs Recovered source (wireframe)")
ax1.set_xlabel("x")
ax1.set_ylabel("t")
ax1.set_zlabel("f/u")

# --- Right: Solution (y) ---
ax2 = fig_final.add_subplot(122, projection='3d')
y_desired = np.load("y_desired.npy")
ax2.plot_surface(X, Tm, y_desired, cmap='viridis', alpha=0.8)
ax2.plot_wireframe(X, Tm, y_new, color='r', alpha=0.6)
ax2.set_title("True solution (surface) vs Recovered solution (wireframe)")
ax2.set_xlabel("x")
ax2.set_ylabel("t")
ax2.set_zlabel("y")

plt.tight_layout()
plt.show()