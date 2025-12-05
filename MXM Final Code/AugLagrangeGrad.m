function [gJ] = AugLagrangeGrad(q, source, params, endpoints, ud, pd, u_max, p_max)
    
    nx=params.nx;
    nt=params.nt;
    scaleu=1;
    scalep=1;

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
    mu1=params.mu{1};   % vector
    mu2=params.mu{2};
    mu3=params.mu{3};
    mu4=params.mu{4};
    rho1=params.rho{1}; % scalar
    rho2=params.rho{2};
    rho3=params.rho{3};
    rho4=params.rho{4};
    [u,p,M,~,matrices]=approxPVEsol(params,source,endpoints,setup,params.nx);
    W=Trap_doub(1:2*nx+2,1:2*nx+2);

    u = u(:).';
    p = p(:).';
    
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

    %%% gAL: ρ1, ρ3
    d1=(mu1+rho1.*(u-u_max));
    d3=mu3+rho3.*(p-p_max);
    d1_use=d1(:,2:end);
    d_diff=[d1_use;d3];
    d13_penalized=max(d_diff,0);
    d13=W*d13_penalized;

    %%% gAL: ρ2, ρ4
    d2=(mu2+rho2.*(-u));
    d4=(mu4+rho4.*(-p));
    d2_use=d2(:,2:end);
    d_weighted=[d2_use;d4];
    d24_penalized=-max(d_weighted,0);
    d24=W*d24_penalized;

    dAL=d13+d24;
    gAL=solve_adj(params,dAL,M,matrices.BcalF,matrices.Bcal);

    gJ = g1 + g2 + gAL;