function [obj] = eval_obj_var(F, source, params, endpoints, ud, pd,Trap_doub,Trap)

source.F = F;

setup = 2;
[u,p] = approxPVEsol(params, source, endpoints, setup, params.nx);

obj=eval_obj(source, u, p, ud, pd,params,Trap_doub,Trap);
end