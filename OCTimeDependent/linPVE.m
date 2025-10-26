function [uhatm,phatm,M,b,matrices]=linPVE(source,params,endpoints,a)
x=linspace(endpoints.xstart,endpoints.xend,params.n+1);
xm=x(1:end-1)/2+x(2:end)/2;

    
if nargin <4
    [matrices,params]=algmatrices(params,2);
    if isfield(params,'k')==0
        matrices=algkmatrices(params,2,matrices,params.kref);
        kfun=0;
    else
        matrices.storeBcal=zeros(2*params.nx+2,2*params.nx+2,params.nt);
        matrices.storeBS=zeros(2*params.nx+2,params.nx+1,params.nt);
        kfun=1;
        matrices.storeM=zeros(2*params.nx+2, 2*params.nx+2, params.nt);
    end
    t=endpoints.tstart+params.deltat;
    uhatm=zeros(params.nx+1,params.n+1);
    phatm=zeros(params.n+1,params.n);
    uhatm(:,1)=source.IC(params.x);
    for k=1:params.nt
        if kfun==1
            if isa(params.k, 'function_handle')==1
                params.kref=params.k(t);
            else
                params.kref=params.k(:,k);
            end
            matrices=algkmatrices(params,2,matrices,params.kref);
        end
        matrices=bmatrices(params,matrices);

        if isa(source.g,'function_handle')==1
            BC.g=source.g(t);
        else
            BC.g=source.g(k);
        end
        
        if isa(source.psi, 'function_handle')==1
            BC.psi=source.psi(t);
        else
            BC.psi=source.psi(k);
        end
        BC.bcu=source.bcu(t);
        BC.bcp=source.bcp(t);
        [b]=derive_b(matrices,source,uhatm(:,k),x,t,2,endpoints,params,BC);
        [uhat,phat,M]=solveuhatphat(matrices.Luu, matrices.Lup,...
        matrices.Lpu, matrices.Lpp,b,params.n);
        uhatm(:,k+1)=uhat;
        phatm(:,k)=phat;
        if kfun==1
            matrices.storeBcal(:,:,k)=matrices.Bcal;
            matrices.storeBS(:,:,k)=matrices.BcalS;
            matrices.storeM(:,:,k)=M;
        end
        %check=M*[uhat;phat]-matrices.Bcal*[uhatm(:,k);zeros(size(phat))]-matrices.BcalF*source.F(:,k)
        t=t+params.deltat;
    end
else
    [matrices,params]=algmatrices(params,2,a);
    matrices=algkmatrices(params,2,matrices,params.kref,a);
    firstuhat=source.ICm(params.x');
    pi=zeros(params.n,1);
    for k=1:params.n
        pi(k)=integral(source.ICh,params.x(k),params.x(k+1))*params.n;
    end
    Pi=solvePi(firstuhat,params.lambdae,matrices.G);
    t=endpoints.tend-params.deltat;
    u=zeros(params.n,params.n);
    p=zeros(params.n,params.n);
    uhatm=zeros(params.n+1,params.n);
    phatm=zeros(params.n+1,params.n);
    u(:,end)=source.ICm(xm);
    p(:,end)=source.ICh(xm);
    for k=1:params.n-1
        matrices=algPivectors(params,2,matrices, source,endpoints,t,Pi,a,pi');
        [bu,bp,matrices]=adjbvectors(params,matrices,params.kref,2);
        BC.g=source.g(t);
        BC.psi=source.psi(t);
        BC.bcu=source.bcu(t);
        BC.bcp=source.bcp(t);
        [uhat,phat,M,b]=solveuhatphat(matrices.Luu, matrices.Lup,...
        matrices.Lhm, matrices.Lpp,bu,bp,BC,params.n);
        Pi=solvePi(uhat,params.lambdae,matrices.G);
        [u(:,end-k),p(:,end-k)]=solveup(uhat,phat,matrices);
        uhatm(:,end-k)=uhat;
        phatm(:,end-k)=phat;
        
        t=t-params.deltat;
    end
end


end