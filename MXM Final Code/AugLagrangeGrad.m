function [gJ] = AugLagrangeGrad(q, source, params, endpoints, ud, pd, u_max, p_max, u_min, p_min)
    % gJ is returned as (nx + 1) X (nt)
    
    % ud, pd, u_max, u_min, p_max, p_min are expected to be stacked vectors
    
    
    
        nx=params.nx;
        nt=params.nt;
        scaleu=1;
        scalep=1;
        u_max=reshape(u_max, nx+1, nt+1);
        p_max=reshape(p_max, nx+1, nt);
        u_min=reshape(u_min, nx+1, nt+1);
        p_min=reshape(p_min, nx+1, nt);
        ud=reshape(ud, nx+1, nt+1);
        pd=reshape(pd, nx+1, nt);
    
        %%% Trap, Trap_doub
        % Need to make these bigger. We need Trap_doub to have 2*n d's and Trap to have n ds.
        d=params.deltax*ones(params.nx+1,1);
        d(1,1)=d(1,1)/2;
        d(params.nx+1,1)=d(params.nx+1,1)/2;
        d=repmat(d,nt,1);
        d=params.deltat*d;
        d2=[scaleu*d;scalep*d];
    
        % The Trap_doub and Trap compute trapezoidal integration with respect to space and right endpoint integration with respect to time.
        Trap_doub=spdiags(d2,0,nt*2*(nx+1),nt*2*(nx+1));
        Trap=spdiags(d,0,nt*(nx+1),nt*(nx+1));
    
        %%% unpacking
        source.F=q;
        setup=2;
        mu1=reshape(params.mu{1}, nx+1, nt+1);   % vector
        mu2=reshape(params.mu{2}, nx+1, nt+1);
        mu3=reshape(params.mu{3}, nx+1, nt);
        mu4=reshape(params.mu{4}, nx+1, nt);
        rho1=params.rho{1}; % scalar
        rho2=params.rho{2};
        rho3=params.rho{3};
        rho4=params.rho{4};
        [u,p,M,~,matrices]=approxPVEsol(params,source,endpoints,setup,params.nx);
        W=Trap_doub(1:2*nx+2,1:2*nx+2);
    
        %%% g1: GTv, v = W(y-yd)
        source.z1=u(:,2:end)-ud(:,2:end);
        source.z2=p-pd;
        v=W*[source.z1;source.z2];
        g1=solve_adj(params,v,M,matrices.BcalF,matrices.Bcal);
    
        %%% g2: betaWq term
        c=(nx+1)*nt;
        F=reshape(source.F,c,1);
        lam_trap_F=params.lambda*Trap*F;
        g2=reshape(lam_trap_F,nx+1,nt);
    
        %%% gAL: rho1, rho3
        d1=(mu1+rho1.*(u-u_max));
        d1_use=d1(:,2:end);
        d3=mu3+rho3.*(p-p_max);
        d_diff=[d1_use;d3];
        d13_penalized=max(d_diff,0);
        d13=W*d13_penalized;
    
        %%% gAL: rho2, rho4
        d2=(mu2+rho2.*(u_min-u));
        d4=(mu4+rho4.*(p_min-p));
        d2_use=d2(:,2:end);
        d_weighted=[d2_use;d4];
        d24_penalized=-max(d_weighted,0);
        d24=W*d24_penalized;
    
        dAL=d13+d24;
        gAL=solve_adj(params,dAL,M,matrices.BcalF,matrices.Bcal);
    
        gJ = g1 + g2 + gAL;
    end