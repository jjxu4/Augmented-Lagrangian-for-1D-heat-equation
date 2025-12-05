function [obj] = eval_obj(source, u, p, ud, pd,params,Trap_doub,Trap)
source.z1 = u(:,2:end) - ud(:,2:end);
source.z2 = p - pd;
v=[source.z1;source.z2];
v_vec=reshape(v,2*(params.nx+1)*params.nt,1);

c=(params.nx+1)*(params.nt);
F=reshape(source.F,c,1);

obj=1/2*v_vec'*Trap_doub*v_vec+params.lambda/2*F'*Trap*F;
end