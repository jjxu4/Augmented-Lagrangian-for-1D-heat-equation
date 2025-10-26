clear all
[params, endpoints] = create_params_optimize();

params.nx = 128;
params.nt = 128;

params.deltat = (endpoints.tend-endpoints.tstart)/params.nt;
t = linspace(endpoints.tstart,endpoints.tend,params.nt+1);
x = linspace(endpoints.xstart,endpoints.xend,params.nx+1);
params.deltax = (endpoints.xend-endpoints.xstart)/params.nx;

%Set parameter to describe case. 1=constantlinear, 2=linear, 3=constantnonlinear, 4=nonlinear.
setup=2;

%To do Test Case 2
[source.F,source.S,g,source.psi,ud_fun,pd_fun] = ...
    tc2source(params, endpoints.xend);
source.g=zeros(params.nt,1);
%%%%%%%%%%%%%


tend = t(2:end);
[X_full,T_full]=meshgrid(x,t);
[X,T] = meshgrid(x, tend);
            
source.bcu=@(y)zeros(size(y));
source.bcp=@(y)zeros(size(y));
source.IC=@(y)zeros(size(y));

%%%%%%%%%%%%%%%%%%%%
[u_tilde,p_tilde]=approxPVEsol(params, source, endpoints, setup, params.nx);

%u_tilde=zeros(params.nx+1,params.nt+1);
p_tilde=zeros(params.nx+1,params.nt);
ud_actual = ud_fun(X_full,T_full)';
pd_actual = pd_fun(X,T)';

ud=ud_actual-u_tilde;
pd=pd_actual-p_tilde;
%%%%%%%%%%%%%%%%%%%


source.S=@(y,t)zeros(size(y));
source.psi=@(t)zeros(size(t));
source.bcu=@(y)zeros(size(y));
source.bcp=@(y)zeros(size(y));
source.z3=@(t)zeros(size(t));
source.z4=@(t)zeros(size(t));
source.F =@(y,t)zeros(size(y));
params.lambda = 1e-1;

source.g=g(tend)';
%ud = ones(params.nx+1);
%pd = ones(params.nx+1,params.nx);


%Need to make these bigger. We need Trap_doub to have 2*n d's and Trap to
%have n ds.
d=params.deltax*ones(params.nx+1,1);
d(1,1)=d(1,1)/2;
d(params.nx+1,1)=d(params.nx+1,1)/2;

d=repmat(d,params.nt,1);
d=params.deltat*d;

e=params.deltat*ones(params.nt,1);

%The Trap_doub and Trap compute trapezoidal integration with respect to
%space and right endpoint integration with respect to time.
Trap_doub=diag([d;d]);
Trap=diag(e);


og = @(g)eval_objgrad_g(g, source, params, endpoints, ud, pd,Trap_doub,Trap);
oo = @(g)eval_obj_var_g(g, source, params, endpoints, ud, pd,Trap_doub, Trap);

% disp('Run with numerical gradient information.')
% options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',false,'CheckGradients',false,'Display','iter','OptimalityTolerance',1e-6 / params.nx);
% fun = oo;
% [x_fd,fval,exitflag,output,grad] = fminunc(fun,source.g,options);

disp('Run with adjoint gradient information.')
options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true,'CheckGradients',false,'Display','iter','OptimalityTolerance',1e-6 / params.nx);
fun = og;
[x_adj,fval,exitflag,output,grad] = fminunc(fun,source.g,options);


source.g=x_adj;
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
plot(tend, x_adj)
hold on 
title('Optimal g with Adjoint Gradient')
xlabel('t')
ylabel('g')

figure(4)
surf(X_full,T_full, u_full,'EdgeColor','none')
hold on
for i = 1 : spacing : length(X(:,1))
    plot3(X_full(:,i), T_full(:,i), u_full(:,i),'-k');
    plot3(X_full(i,:), T_full(i,:), u_full(i,:),'-k');
end
title('$u+\tilde{u}$', 'Interpreter','latex')
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
title('$p+\tilde{p}$', 'Interpreter', 'latex')
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
