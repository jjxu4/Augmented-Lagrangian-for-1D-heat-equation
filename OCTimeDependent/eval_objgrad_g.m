function [f,g] = eval_objgrad_g(g, source, params, endpoints, ud, pd,Trap_doub,Trap)

source.g = g;

setup = 2;
[u,p,M,~,matrices] = approxPVEsol(params, source, endpoints, setup, params.nx);
source.z1 = u(:,2:end) - ud(:,2:end);
source.z2 = p - pd;
v=Trap_doub(1:2*params.nx+2,1:2*params.nx+2)*[source.z1;source.z2];
control_mat=zeros(2*params.nx+2,1);
control_mat(params.nx+1)=1;
adj_var=solve_adj(params, v,M,control_mat,matrices.Bcal);
adj_var=adj_var(end,:);

f = eval_obj_g(source, u, p, ud, pd, params,Trap_doub,Trap);

lam_trap_g=params.lambda*Trap*source.g;
g = adj_var'+lam_trap_g;
end