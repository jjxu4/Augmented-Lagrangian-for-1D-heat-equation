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

t = linspace(endpoints.tstart, endpoints.tend, nt+1);   % 1 x (nt+1)
x = linspace(endpoints.xstart,  endpoints.xend,  nx+1); % 1 x (nx+1)

t_interior = t(2:end);  
[X_full, T_full] = meshgrid(x, t);      
[X, T] = meshgrid(x, t_interior);

setup = 2;  

% Get ARMA461_new test case and build desired states u_d, p_d
[F_fun, S_fun, g_fun, psi_fun, exactu_fun, exactp_fun, ~] = ARMA461_new(params, endpoints.xend);

% Desired states
u_d_actual = exactu_fun(X_full, T_full)';   % (nx+1) x (nt+1)
p_d_actual = exactp_fun(X, T)';   % (nx+1) x nt

% These are the targets the unconstrained optimization will try to reach
ud = u_d_actual;
pd = p_d_actual;

% Define PDE sources for optimization problem 
source.S = @(y,t) zeros(size(y));  
source.psi = @(t) zeros(size(t)); 
source.g = @(t) zeros(size(t)); 
source.bcu = @(y) zeros(size(y)); 
source.bcp = @(y) zeros(size(y)); 
source.z3 = @(t) zeros(size(t)); 
source.z4 = @(t) zeros(size(t));
source.IC = @(y) zeros(size(y)); 

% Initial guess for control
F0 = zeros(nx+1, nt);
source.F = F0;  

% Control regularization parameter
params.lambda = 1e-5;

scaleu = 1;
scalep = 1;
d = params.deltax * ones(nx+1,1);
d(1) = d(1)/2;
d(end) = d(end)/2;
d = repmat(d, nt, 1);     
d = params.deltat * d;   
d2 = [scaleu * d; scalep * d];

Trap_doub = spdiags(d2, 0, nt*2*(nx+1), nt*2*(nx+1));
Trap = spdiags(d,  0, nt*(nx+1),   nt*(nx+1));

og = @(F) eval_objgrad_t(F, source, params, endpoints, ud, pd, Trap_doub, Trap);
oo = @(F) eval_obj_var(F, source, params, endpoints, ud, pd, Trap_doub, Trap); 

% Unconstrained optimization with fminunc 
options = optimoptions('fminunc', ...
    'Algorithm',              'quasi-newton', ...
    'SpecifyObjectiveGradient', true, ...
    'CheckGradients',         false, ...
    'Display',               'iter', ...
    'OptimalityTolerance',   1e-6 / params.nx);

fun = og;
[q_opt, fval, exitflag, output, grad] = fminunc(fun, source.F, options);

% fprintf('fminunc exitflag = %d\n', exitflag);
% disp(output);

% forward solve to get optimal states
source.F = q_opt;
[u_opt, p_opt] = approxPVEsol(params, source, endpoints, setup, params.nx);

% Prepare grids for plotting
X_full_plot = X_full';
T_full_plot = T_full';
X_plot = X';
T_plot = T';

u_d_plot = u_d_actual;
p_d_plot = p_d_actual;

%% 8.1 overlay results with state constraints that WILL BE used in the constrained version
% Define state constraints

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

% Plot results

% u_d
figure(201)
surf(X_full_plot, T_full_plot, u_d_plot, 'EdgeColor', 'none');
title('u_des');
xlabel('x'); ylabel('t'); zlabel('u_d');
colorbar; rotate3d on; view(135,30);

% p_d
figure(202)
surf(X_plot, T_plot, p_d_plot, 'EdgeColor', 'none');
title('p_des');
xlabel('x'); ylabel('t'); zlabel('p_d');
colorbar; rotate3d on; view(135,30);

% q_opt
figure(203)
surf(X_plot, T_plot, q_opt, 'EdgeColor', 'none');
title('q_{opt} (Unconstrained)');
xlabel('x'); ylabel('t'); zlabel('q_{opt}');
colorbar; rotate3d on; view(135,30);

% u_opt overlaid with constraints
figure(204)
surf(X_full_plot, T_full_plot, u_opt, 'EdgeColor', 'none');
hold on
surf(X_full_plot, T_full_plot, u_max_grid, ...
    'FaceColor', 'none', 'EdgeColor', 'k');
surf(X_full_plot, T_full_plot, u_min_grid, ...
    'FaceColor', 'none', 'EdgeColor', 'k');
hold off
title('u_{opt} (Unconstrained)');
xlabel('x'); ylabel('t'); zlabel('u');
colorbar; rotate3d on; view(135,30);

% p_opt overlaid with constraints
figure(205)
surf(X_plot, T_plot, p_opt, 'EdgeColor', 'none');
hold on
surf(X_plot, T_plot, p_max_grid, ...
    'FaceColor', 'none', 'EdgeColor', 'k');
surf(X_plot, T_plot, p_min_grid, ...
    'FaceColor', 'none', 'EdgeColor', 'k');
hold off
title('p_{opt} (Unconstrained)');
xlabel('x'); ylabel('t'); zlabel('p');
colorbar; rotate3d on; view(135,30);
