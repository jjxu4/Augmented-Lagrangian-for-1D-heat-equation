% main_new.m – AL version using ARMA461 source like run_opt_new

% Important notes:
% 1. mu{1}, ..., mu{4} must be column vectors.

clear variables
close all

[params, endpoints] = create_params_optimize_new();

% gridding parameters
params.nx = 24;   % you can change to 96 like run_opt_new if desired
params.nt = 24;
nx = params.nx;
nt = params.nt;

params.deltat = (endpoints.tend - endpoints.tstart)/params.nt;
params.deltax = (endpoints.xend - endpoints.xstart)/params.nx;

t = linspace(endpoints.tstart, endpoints.tend, params.nt+1);   % 1 x (nt+1)
x = linspace(endpoints.xstart, endpoints.xend, params.nx+1);  % 1 x (nx+1)

t_interior = t(2:end);   % nt time points for F, p, etc.
[X_full, T_full] = meshgrid(x, t);         % (nt+1) x (nx+1)
[X, T]           = meshgrid(x, t_interior);% nt x (nx+1)

% -------------------------------------------------------------------------
% Use ARMA461 in the SAME WAY as run_opt_new
% -------------------------------------------------------------------------
% NOTE: run_opt_new uses:
% [F, source.S, source.g, source.psi, ud_fun, pd_fun] = ARMA461(...)
% so ARMA461 has 6 outputs (not 7).
[F_fun, S_fun, g_fun, psi_fun, ud_fun, pd_fun] = ARMA461(params, endpoints.xend);

% Build q_actual as a grid (nx+1) x nt
% F_fun expects (y,t); we evaluate at our space-time grid for t(2:end)
q_actual = F_fun(X, T)';  % (nx+1) x nt

% Define source dictionary
source.S   = S_fun;          % function handle S(y,t)
source.g   = g_fun;          % boundary flux/traction on u
source.psi = psi_fun;        % boundary data for p

source.bcu = @(y) zeros(size(y));  % u boundary condition (can refine if needed)
source.bcp = @(y) zeros(size(y));  % p boundary condition
source.IC  = @(y) zeros(size(y));  % initial condition (matches exactu at t=0)

% Set the control to q_actual for the "truth" run
source.F = q_actual;

setup = 2;  % linear PVE

% -------------------------------------------------------------------------
% Compute u_tilde, p_tilde as in run_opt_new
% -------------------------------------------------------------------------
[u_tilde, p_tilde] = approxPVEsol(params, source, endpoints, setup);
% u_tilde: (nx+1) x (nt+1)
% p_tilde: (nx+1) x nt

% -------------------------------------------------------------------------
% Build desired states u_desired, p_desired the SAME WAY as run_opt_new
% -------------------------------------------------------------------------
% In run_opt_new:
%   ud = min(u_tilde, 0.0046);
%   pd = max(p_tilde, -0.0016);
u_desired = min(u_tilde, 0.0046);
p_desired = max(p_tilde, -0.0016);

% If you want to also have the analytic ud_fun/pd_fun (like run_opt_new's ud_actual, pd_actual),
% you could do:
% ud_actual = ud_fun(X_full, T_full)';   % (nx+1) x (nt+1)
% pd_actual = pd_fun(X, T)';             % (nx+1) x nt
% but for the AL cost, we use u_desired, p_desired above.

% -------------------------------------------------------------------------
% Now we initialize our OuterLoop call
% -------------------------------------------------------------------------
params.u_d = u_desired(:);  % vector length (nx+1)*(nt+1)
params.p_d = p_desired(:);  % vector length (nx+1)*nt

Uscale = max(1, max(abs(u_desired(:))));
Pscale = max(1, max(abs(p_desired(:))));

% state constraints big enough to be irrelevant for now
u_max_grid =  7 * ones(size(u_desired));
u_min_grid = -7 * ones(size(u_desired));
p_max_grid =  7 * ones(size(p_desired));
p_min_grid = -7 * ones(size(p_desired));

params.u_max = u_max_grid(:);
params.u_min = u_min_grid(:);
params.p_max = p_max_grid(:);
params.p_min = p_min_grid(:);

params.lambda = 1e-5;    % regularization
params.beta   = params.lambda;
params.W      = params.deltax * params.deltat;   % simple scalar weight

% AL multipliers and penalties
rho1 = 1.0; rho2 = 1.0; rho3 = 1.0; rho4 = 1.0;
params.rho_initial = {rho1, rho2, rho3, rho4};

mu1 = zeros((nx+1)*(nt+1), 1);
mu2 = zeros((nx+1)*(nt+1), 1);
mu3 = zeros((nx+1)*nt, 1);
mu4 = zeros((nx+1)*nt, 1);
params.mu_initial = {mu1, mu2, mu3, mu4};

params.gamma   = 2;
params.tau     = 0.9;
params.epsilon = 1e-3;

% Initial residual scale. same as the one used in 2024 paper
params.R = 1e3; 

% Run OuterLoop ===========================================================
% Initial guess for q
q_initial = zeros(nx+1, nt);

q_optimal = OuterLoop(q_initial, params, source, endpoints);

% Compare q_optimal vs q_actual
diff_q = q_optimal - q_actual;
rel_err_q = norm(diff_q(:)) / norm(q_actual(:));
fprintf('Relative error in q: %e\n', rel_err_q);

%% Everything after this is for plotting. TODO: make good looking plots
% Recover optimal states with q_optimal
source.F = q_optimal;
[u_optimal, p_optimal] = approxPVEsol(params, source, endpoints, setup);

% Prepare plotting grids ==================================================
X_full_plot = X_full';   
T_full_plot = T_full';   
X_plot = X';             
T_plot = T';             

% -------------------------------------------------------------------------
% 1. Plot q_actual
% -------------------------------------------------------------------------
figure(1);
surf(X_plot, T_plot, q_actual, 'EdgeColor', 'none');
title('q\_actual');
xlabel('x'); ylabel('t'); zlabel('q');
pbaspect([1 1 1])
view(135, 30);
colorbar;
rotate3d on;

% -------------------------------------------------------------------------
% 2A. Plot u_desired with bounds (SEPARATE FIGURE)
% -------------------------------------------------------------------------
figure(2);
surf(X_full_plot, T_full_plot, u_desired, 'EdgeColor', 'none');  
hold on;
surf(X_full_plot, T_full_plot, u_max_grid, 'FaceColor', 'none', 'EdgeColor', 'k');
surf(X_full_plot, T_full_plot, u_min_grid, 'FaceColor', 'none', 'EdgeColor', 'k');
title('u\_desired with bounds');
xlabel('x'); ylabel('t'); zlabel('u');
legend({'u\_desired','u\_max','u\_min'}, 'Location', 'best');
pbaspect([1 1 1])
view(135, 30);
colorbar;
rotate3d on;

% -------------------------------------------------------------------------
% 2B. Plot p_desired with bounds (SEPARATE FIGURE)
% -------------------------------------------------------------------------
figure(3);
surf(X_plot, T_plot, p_desired, 'EdgeColor', 'none');  
hold on;
surf(X_plot, T_plot, p_max_grid, 'FaceColor', 'none', 'EdgeColor', 'k');
surf(X_plot, T_plot, p_min_grid, 'FaceColor', 'none', 'EdgeColor', 'k');
title('p\_desired with bounds');
xlabel('x'); ylabel('t'); zlabel('p');
legend({'p\_desired','p\_max','p\_min'}, 'Location', 'best');
pbaspect([1 1 1])
view(135, 30);
colorbar;
rotate3d on;

% -------------------------------------------------------------------------
% 3. Plot q_optimal
% -------------------------------------------------------------------------
figure(4);
surf(X_plot, T_plot, q_optimal, 'EdgeColor', 'none');
title('q\_optimal');
xlabel('x'); ylabel('t'); zlabel('q');
pbaspect([1 1 1])
view(135, 30);
colorbar;
rotate3d on;

% -------------------------------------------------------------------------
% 4A. Plot u_optimal with bounds (SEPARATE FIGURE)
% -------------------------------------------------------------------------
figure(5);
surf(X_full_plot, T_full_plot, u_optimal, 'EdgeColor', 'none');  
hold on;
surf(X_full_plot, T_full_plot, u_max_grid, 'FaceColor', 'none', 'EdgeColor', 'k');
surf(X_full_plot, T_full_plot, u_min_grid, 'FaceColor', 'none', 'EdgeColor', 'k');
title('u\_optimal with bounds');
xlabel('x'); ylabel('t'); zlabel('u');
legend({'u\_optimal','u\_max','u\_min'}, 'Location', 'best');
pbaspect([1 1 1])
view(135, 30);
colorbar;
rotate3d on;

% -------------------------------------------------------------------------
% 4B. Plot p_optimal with bounds (SEPARATE FIGURE)
% -------------------------------------------------------------------------
figure(6);
surf(X_plot, T_plot, p_optimal, 'EdgeColor', 'none');  
hold on;
surf(X_plot, T_plot, p_max_grid, 'FaceColor', 'none', 'EdgeColor', 'k');
surf(X_plot, T_plot, p_min_grid, 'FaceColor', 'none', 'EdgeColor', 'k');
title('p\_optimal with bounds');
xlabel('x'); ylabel('t'); zlabel('p');
legend({'p\_optimal','p\_max','p\_min'}, 'Location', 'best');
pbaspect([1 1 1])
view(135, 30);
colorbar;
rotate3d on;
