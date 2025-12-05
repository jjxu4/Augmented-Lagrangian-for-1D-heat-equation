function adj_var=solve_adj(params,v,M,Betauf,B)
%Solving the adjoint for the time dependent case 
%   Detailed explanation goes here
Bdim=ndims(B);
sol=zeros(2*params.nx+2,params.nt);
vec=v(:,end);
if Bdim==3
    sol(:,end)=M(:,:,params.nt)'\vec;
    for i=1:params.nt-1
        vec=v(:,end-i)+B(:,:,params.nt-i+1)'*sol(:,end-i+1);
        sol(:,end-i)=M(:,:,params.nt-i)'\vec;
    end
else
    sol(:,end)=M'\vec;
    for i=1:params.nt-1
        vec=v(:,end-i)+B'*sol(:,end-i+1);
        sol(:,end-i)=M'\vec;
    end
end

adj_var=Betauf'*sol;
    
end

