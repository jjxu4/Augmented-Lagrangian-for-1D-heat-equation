% Define params
clear variables
close all
[params, endpoints] = create_params_optimize();

%If you want k to be a function of time define params.k. If not, then
%simply defining params.kref is enough.
params.kref=1;
%params.k = @(t)4*t+1;
%params.k = @(t)t;

params.nx = 24;
params.nt = 24;
total_u = (params.nx+1) * (params.nt+1);
total_p = (params.nx+1) * (params.nt);
params.u_max=ones(1, total_u)*1e3;
params.p_max=ones(1, total_p)*1e3;
params.gamma=2;
params.tau=0.9;
params.R=1e5;
params.epsilon=1e-4;

% placeholder
params.u_d=ones(1, total_u);
params.p_d=ones(1, total_p);
params.mu_initial_u=zeros(1, total_u);
params.mu_initial_p=zeros(1, total_p);
params.rho_initial=1;

%if you want k as a function of x and t you have to put in a matrix
%params.k=5*rand(params.nx, params.nt);


% Parameters to scale a constant in front of u-ud term in cost function and
% in front of p-pd term in cost function
scaleu=1;
scalep=1;

%Set parameter to describe case. 1=constantlinear, 2=linear, 3=constantnonlinear, 4=nonlinear.
setup=2;

params.deltat = (endpoints.tend-endpoints.tstart)/params.nt;
t = linspace(endpoints.tstart,endpoints.tend,params.nt+1);
x = linspace(endpoints.xstart,endpoints.xend,params.nx+1);
params.deltax = (endpoints.xend-endpoints.xstart)/params.nx;
params.W = params.deltax*params.deltat;
params.beta=1;
params.mu={params.mu_initial_u, params.mu_initial_u, params.mu_initial_p, params.mu_initial_p};
params.mu_initial={params.mu_initial_u, params.mu_initial_u, params.mu_initial_p, params.mu_initial_p};
params.rho={params.rho_initial, params.rho_initial, params.rho_initial, params.rho_initial};
params.lambda = 1e-5;

%Find u tilde and p tilde which represent the solution with the control set
%to zero adn all other conditions set as desired.
[F,source.S,source.g,source.psi,ud_fun,pd_fun] = tc2source(params, endpoints.xend);

q0 = ones(1, total_p);
source.F=q0;

tend = t(2:end);
[X_full,T_full]=meshgrid(x,t);
[X,T] = meshgrid(x, tend);
%source.F = F(X, T)';
            
source.bcu=@(y)zeros(size(y));
source.bcp=@(y)zeros(size(y));
source.IC=@(y)zeros(size(y));

% Call OuterLoop
q_opt = OuterLoop(q0, params, source, endpoints);