function [obj] = eval_obj_var_g(g, source, params, endpoints, ud, pd,Trap_doub,Trap)

source.g = g;

setup = 2;
[u,p] = approxPVEsol(params, source, endpoints, setup, params.nx);

obj=eval_obj_g(source, u, p, ud, pd,params,Trap_doub,Trap);
end