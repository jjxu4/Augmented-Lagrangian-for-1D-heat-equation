clear variables
close all

% Parameters and gridding
[params, endpoints] = create_params_optimize_new();

params.nx = 24;
params.nt = 24;
nx = params.nx;
nt = params.nt;

params.deltat = (endpoints.tend - endpoints.tstart) / params.nt;
params.deltax = (endpoints.xend  - endpoints.xstart) / params.nx;

t = linspace(endpoints.tstart, endpoints.tend, nt+1);   
x = linspace(endpoints.xstart,  endpoints.xend,  nx+1); 

t_interior = t(2:end);  

[X_full, T_full] = meshgrid(x, t);       
[X, T] = meshgrid(x, t_interior);

X_full_plot = X_full';   
T_full_plot = T_full';
X_plot = X';       
T_plot = T';

setup = 2;  

% Build desired states u_d, p_d from ARMA461_new

[F_fun, S_fun, g_fun, psi_fun, exactu_fun, exactp_fun, ~] = ARMA461_new(params, endpoints.xend);

% Desired states 
u_desired = exactu_fun(X_full, T_full)';   % (nx+1) x (nt+1)
p_desired = exactp_fun(X, T)';   % (nx+1) x nt

params.u_d = u_desired(:);
params.p_d = p_desired(:);

% State constraints

% p constraints
p_max_x = 1.1 * ones(size(x));          
p_max_x(abs(x) > 0.4) = 0.2;         

p_min_x = -0.15 * ones(size(x));        

p_max_grid = repmat(p_max_x.', 1, nt);   % (nx+1) x nt
p_min_grid = repmat(p_min_x.', 1, nt);   % (nx+1) x nt

% u constraints
u_max_x = 1.1 * ones(size(x));
u_max_x(abs(x) > 0.5) = 0.1;

u_min_x = -0.1 * ones(size(x));

u_max_grid = repmat(u_max_x.', 1, nt+1); % (nx+1) x (nt+1)
u_min_grid = repmat(u_min_x.', 1, nt+1);

params.u_max = u_max_grid(:);
params.u_min = u_min_grid(:);
params.p_max = p_max_grid(:);
params.p_min = p_min_grid(:);

% AL and cost parameters
params.lambda = 1e-5;      
params.beta = params.lambda; 
params.W = params.deltax * params.deltat;

% Initial multipliers
mu1 = zeros((nx+1)*(nt+1), 1);
mu2 = zeros((nx+1)*(nt+1), 1);
mu3 = zeros((nx+1)*nt,     1);
mu4 = zeros((nx+1)*nt,     1);

params.mu_initial = {mu1, mu2, mu3, mu4};

% Initial penalty parameters
rho1 = 1.0; rho2 = 1.0; rho3 = 1.0; rho4 = 1.0;
params.rho_initial = {rho1, rho2, rho3, rho4};

% Algorithm 2 parameters
params.gamma = 2;     
params.tau = 0.9;   
params.epsilon = 1e-3;   
params.R = 1e3;    

% Define sources for constrained problem
source.S   = @(y,t) zeros(size(y));
source.psi = @(t) zeros(size(t));
source.g   = @(t) zeros(size(t));
source.bcu = @(y) zeros(size(y));
source.bcp = @(y) zeros(size(y));
source.IC  = @(y) zeros(size(y));
source.z3  = @(t) zeros(size(t));
source.z4  = @(t) zeros(size(t));

% Run OuterLoop
q_initial = zeros(nx+1, nt); % initial guess for control F
q_optimal = OuterLoop(q_initial, params, source, endpoints);

% Recover optimal states
source.F = q_optimal;
[u_opt, p_opt] = approxPVEsol(params, source, endpoints, setup);


%% EVERYTHING AFTER THIS IS PLOTTING
% Plot optimal control q_optimal
spacing = floor(nx / 12);

figure(301)
surf(X_plot, T_plot, q_optimal, 'EdgeColor', 'none');
hold on
[~, nCols] = size(X_plot);
for i = 1:spacing:nCols
    plot3(X_plot(:,i), T_plot(:,i), q_optimal(:,i), '-k');
    plot3(X_plot(i,:), T_plot(i,:), q_optimal(i,:), '-k');
end
hold off
title('q_{opt} (Constrained)');
xlabel('x'); ylabel('t'); zlabel('F');
colorbar; rotate3d on; view(135,30);

% Plot u_desired (Arma_unconstrained will plot this overlaid with
% constraints)

figure(302)
surf(X_full_plot, T_full_plot, u_desired, 'EdgeColor', 'none');
title('u_des');
xlabel('x'); ylabel('t'); zlabel('u_d');
colorbar; rotate3d on; view(135,30);

% Plot p_desired (Arma_unconstrained will plot this overlaid with
% constraints)

figure(303)
surf(X_plot, T_plot, p_desired, 'EdgeColor', 'none');
title('p_des');
xlabel('x'); ylabel('t'); zlabel('p_d');
colorbar; rotate3d on; view(135,30);

% Plot u_opt overlaid with constraints (wireframe)

figure(304)
surf(X_full_plot, T_full_plot, u_opt, 'EdgeColor', 'none');
hold on
surf(X_full_plot, T_full_plot, u_max_grid, ...
    'FaceColor', 'none', 'EdgeColor', 'k');
surf(X_full_plot, T_full_plot, u_min_grid, ...
    'FaceColor', 'none', 'EdgeColor', 'k');
hold off
title('u_{opt} (Constrained)');
xlabel('x'); ylabel('t'); zlabel('u');
colorbar; rotate3d on; view(135,30);

% Plot p_opt overlaiad with constraints (wireframe)

figure(305)
surf(X_plot, T_plot, p_opt, 'EdgeColor', 'none');
hold on
surf(X_plot, T_plot, p_max_grid, ...
    'FaceColor', 'none', 'EdgeColor', 'k');
surf(X_plot, T_plot, p_min_grid, ...
    'FaceColor', 'none', 'EdgeColor', 'k');
hold off
title('p_{opt} (Constrained)');
xlabel('x'); ylabel('t'); zlabel('p');
colorbar; rotate3d on; view(135,30);
