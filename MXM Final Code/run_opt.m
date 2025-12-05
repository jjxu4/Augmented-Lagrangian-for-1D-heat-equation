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


%Find u tilde and p tilde which represent the solution with the control set
%to zero adn all other conditions set as desired.
[F,source.S,source.g,source.psi,ud_fun,pd_fun] = ...
    tc2source(params, endpoints.xend);
source.F=zeros(params.nx,params.nt);
source.F=@(x,t)0;


tend = t(2:end);
[X_full,T_full]=meshgrid(x,t);
[X,T] = meshgrid(x, tend);
%source.F = F(X, T)';
            
source.bcu=@(y)zeros(size(y));
source.bcp=@(y)zeros(size(y));
source.IC=@(y)zeros(size(y));

[u_tilde,p_tilde]=approxPVEsol(params, source, endpoints, setup);


ud_actual = ud_fun(X_full,T_full)';%ones(params.nx+1,params.nt+1);
pd_actual = pd_fun(X,T)';%ones(params.nx+1,params.nt);

ud=ud_actual-u_tilde;
pd=pd_actual-p_tilde;

source.S=@(y,t)zeros(size(y));
source.psi=@(t)zeros(size(t));
source.g=@(t)zeros(size(t));
source.bcu=@(y)zeros(size(y));
source.bcp=@(y)zeros(size(y));
source.z3=@(t)zeros(size(t));
source.z4=@(t)zeros(size(t));
source.F = source.S(X, T)';
params.lambda = 1e-5;

%source.ICm=@(y)zeros(size(y));
%source.ICh=@(y)zeros(size(y));

%Need to make these bigger. We need Trap_doub to have 2*n d's and Trap to
%have n ds.
d=params.deltax*ones(params.nx+1,1);
d(1,1)=d(1,1)/2;
d(params.nx+1,1)=d(params.nx+1,1)/2;

d=repmat(d,params.nt,1);
d=params.deltat*d;

d2=[scaleu*d;scalep*d];
%The Trap_doub and Trap compute trapezoidal integration with respect to
%space and right endpoint integration with respect to time.
Trap_doub=spdiags(d2,0,params.nt*2*(params.nx+1),params.nt*2*(params.nx+1));
Trap=spdiags(d,0,params.nt*(params.nx+1),params.nt*(params.nx+1));


og = @(F)eval_objgrad_t(F, source, params, endpoints, ud, pd,Trap_doub,Trap);
oo = @(F)eval_obj_var(F, source, params, endpoints, ud, pd,Trap_doub, Trap);

% disp('Run with numerical gradient information.')
% options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',false,'CheckGradients',false,'Display','iter','OptimalityTolerance',1e-6 / params.nx);
% fun = oo;
% [x_fd,~,~,~,~] = fminunc(fun,source.F,options);

disp('Run with adjoint gradient information.')
options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true,'CheckGradients',false,'Display','iter','OptimalityTolerance',1e-6 / params.nx);
fun = og;
[x_adj,fval,exitflag,output,grad] = fminunc(fun,source.F,options);

source.F=x_adj;
[u,p] = approxPVEsol(params, source, endpoints, setup, params.nx);

% figure(2)
% 
% surf(X,T, x_fd')
% hold on
% title('Optimal F with Numerical Gradient')


X=X';
T=T';
X_full=X_full';
T_full=T_full';

p_full=p+p_tilde;
u_full=u+u_tilde;



    
spacing =floor(params.nx/12);


figure(3)
surf(X,T, x_adj,'EdgeColor','none')
hold on 
for i = 1 : spacing : params.nx
    plot3(X(:,i), T(:,i), x_adj(:,i),'-k');
    plot3(X(i,:), T(i,:), x_adj(i,:),'-k');
end
title('Optimal F')
xlabel('x')
ylabel('t')
zlabel('F')

figure(4)
surf(X_full,T_full, u_full,'EdgeColor','none')
hold on
for i = 1 : spacing : length(X(:,1))
    plot3(X_full(:,i), T_full(:,i), u_full(:,i),'-k');
    plot3(X_full(i,:), T_full(i,:), u_full(i,:),'-k');
end
title('$u_{\bar{q}}+\tilde{u}$', 'Interpreter','latex')
xlabel('x')
ylabel('t')
zlabel('u')

figure(5)
surf(X,T,p_full, 'EdgeColor','none')
hold on
for i = 1 : spacing : params.nx
    plot3(X(:,i), T(:,i), p_full(:,i),'-k');
    plot3(X(i,:), T(i,:), p_full(i,:),'-k');
end
title('$p_{\bar{q}}+\tilde{p}$', 'Interpreter', 'latex')
xlabel('x')
ylabel('t')
zlabel('p')

figure(6)
surf(X_full,T_full,ud_actual, 'EdgeColor','none')
hold on
for i = 1 : spacing : params.nx
    plot3(X_full(:,i), T_full(:,i), ud_actual(:,i),'-k');
    plot3(X_full(i,:), T_full(i,:), ud_actual(i,:),'-k');
end
title('$u_d$', 'Interpreter', 'latex')
xlabel('x')
ylabel('t')
zlabel('u')

figure(7)
surf(X,T,pd_actual, 'EdgeColor','none')
hold on
for i = 1 : spacing : params.nx
    plot3(X(:,i), T(:,i), pd_actual(:,i),'-k');
    plot3(X(i,:), T(i,:), pd_actual(i,:),'-k');
end
title('$p_d$', 'Interpreter', 'latex')
xlabel('x')
ylabel('t')
zlabel('p')

figure(8)
F=F(X,T);
surf(X,T,F, 'EdgeColor','none')
hold on
for i = 1 : spacing : params.nx
    plot3(X(:,i), T(:,i), F(:,i),'-k');
    plot3(X(i,:), T(i,:), F(i,:),'-k');
end
title('Actual F')
xlabel('x')
ylabel('t')
zlabel('F')

