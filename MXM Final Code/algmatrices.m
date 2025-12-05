function [matrices,params]=algmatrices(params,setup,a)
%Defines Ai, B, C, D, G, MM,scriptB,scriptBuinv, Q,QQ,Qvec,Lup
    phi1=@(z)1-z/params.h;
    phi2=@(z)z/params.h;

    %product of basis functions
    phi11=@(z)phi1(z).^2;
    phi12=@(z)phi1(z).*phi2(z);
    phi22=@(z)phi2(z).^2;

    aa=integral(phi11,0,params.h); %integral of product of first basis function on first interval
    ab=integral(phi12,0,params.h); %integral of product of first two basis functions on first interval
    bb=integral(phi22,0,params.h); %integral of product of second basis function on first interval
    
    params.mu = 2. * params.mue ...
    + 2.*params.delta*params.muv/params.deltat;
    params.mp = 1+params.delta...
        * params.lambdav / (params.deltat*params.lambdae);

%   Matrices not specific to case
    matrices.Ai=1/(aa*bb-ab^2)*[bb,-1*ab;-1*ab,aa];
    matrices.B=[-1,1];
    matrices.C=[-1,0;0,1];
    matrices.D=[params.h/2;params.h/2];
    matrices.G=[-1/params.h,1/params.h];
    matrices.MM=matrices.Ai*(params.mu*matrices.C...
        +params.lambdae*params.mp*matrices.D*matrices.G);
    matrices.scriptB=matrices.B*matrices.Ai*matrices.B.';
    
 

    matrices.scriptBuinv=1/(params.mu*matrices.scriptB);
    matrices.Q=1/matrices.scriptB*(matrices.B*matrices.Ai*matrices.C);
    matrices.Qvec=reshape(matrices.Q',1,numel(matrices.Q));
    DQ=matrices.D*matrices.Qvec;
    
    if nargin>2
        if (-1)^setup==1
            matrices.QQ=matrices.scriptBuinv/params.deltat*matrices.B*matrices.Ai*matrices.D*matrices.Q';
        else
            matrices.QQ=zeros(2,1);
        end
    else
        matrices.QQ=-matrices.scriptBuinv*matrices.B*matrices.Ai*matrices.D*matrices.Q';
    end
        %Vectorize Lup
    QQvec=reshape(matrices.QQ,1,numel(matrices.QQ));
    matrices.Lup=-matrices.Ai*params.mu*matrices.B'*QQvec;
    if nargin>2
        if (-1)^setup==1
            matrices.Lup=matrices.Lup+1/params.deltat*matrices.Ai*DQ;
        end
          matrices.Lhm=-matrices.Ai*matrices.B'*1/matrices.scriptB*params.h*matrices.G;
          matrices.Lhm=repmat(matrices.Lhm,1,params.nx);
    else
        matrices.Lup=matrices.Lup-matrices.Ai*(DQ);
    end
    matrices.Lup=repmat(matrices.Lup,1,params.nx);
end