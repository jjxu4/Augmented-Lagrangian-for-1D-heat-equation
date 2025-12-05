function [f,g] = eval_objgrad_t(F, source, params, endpoints, ud, pd,Trap_doub,Trap)

source.F = F;

setup = 2;
[u,p,M,~,matrices] = approxPVEsol(params, source, endpoints, setup, params.nx);
source.z1 = u(:,2:end) - ud(:,2:end);
source.z2 = p - pd;
v=Trap_doub(1:2*params.nx+2,1:2*params.nx+2)*[source.z1;source.z2];
if isfield(matrices,'storeBcal')==1
    adj_var=solve_adj(params, v,matrices.storeM,matrices.BcalF,matrices.storeBcal);
else
    adj_var=solve_adj(params, v,M,matrices.BcalF,matrices.Bcal);
end


f = eval_obj(source, u, p, ud, pd, params,Trap_doub,Trap);

c=(params.nx+1)*params.nt;
F=reshape(source.F,c,1);
lam_trap_F=params.lambda*Trap*F;
g2=reshape(lam_trap_F,params.nx+1,params.nt);
g = adj_var+g2;
end