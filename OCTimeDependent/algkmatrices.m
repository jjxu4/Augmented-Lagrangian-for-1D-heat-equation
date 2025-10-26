function [matrices]=algkmatrices(params,setup,matrices,k,a)
%Defines scriptBpinv, R, RR, Luu, AiBR,Lpu, Lpp, 
matrices.scriptBpinv=(1./(k.*matrices.scriptB));
    if (setup==1 && nargin<5)
        matrices.R=[0,0];
    elseif (setup==3 && nargin<5)
        matrices.R=zeros(params.n,2);
    else
        matrices.R=-params.h*matrices.scriptBpinv*matrices.G;
        if nargin<5
            matrices.R=matrices.R./params.deltat;
%         else
%             matrices.R=-matrices.R;
        end
    end
%Create RR
AiDR1=matrices.Ai*matrices.D*matrices.R(1:end,1)';
AiDR2=matrices.Ai*matrices.D*matrices.R(1:end,2)';
bigAiDR = reshape([AiDR1;AiDR2], size(AiDR1,1), []);
[~,c]=size(bigAiDR);
matrices.bigMM=repmat(matrices.MM,1,c/2);
matrices.RR=matrices.scriptBuinv*matrices.B*matrices.bigMM;
if nargin<5
    matrices.RR=matrices.RR-matrices.scriptBuinv*matrices.B*bigAiDR;
elseif (-1)^setup==1
    matrices.RR=matrices.RR-1/params.deltat*matrices.scriptBuinv*matrices.B*bigAiDR;
end


%Vectorize Luu and Lpu for PVE model
matrices.Rvec=reshape(matrices.R',1,numel(matrices.R));
DR=matrices.D*matrices.Rvec;
Luu1=matrices.bigMM-matrices.Ai*(params.mu*matrices.B'*matrices.RR);

kjvec=reshape([k';k'],size(k',1),[]);
if nargin<5
    matrices.Luu=Luu1-matrices.Ai*DR;
    matrices.AiBR=matrices.Ai*matrices.B.'*matrices.Rvec;
    matrices.Lpu=kjvec.*matrices.AiBR;
    matrices.Lpu=repmat(matrices.Lpu,1,2*params.nx/c);    
else
%     matrices.Lhm=matrices.B'*matrices.R;
%     matrices.Lhm=repmat(matrices.Lhm,1,params.n);
    matrices.Luu=Luu1+1/params.deltat*matrices.Ai*DR;
end
[~,c]=size(DR);
matrices.Luu=repmat(matrices.Luu,1,2*params.nx/c);

%Vectorize Lpp
Lppnok=matrices.Ai*matrices.C-matrices.Ai*matrices.B'*matrices.Qvec;
BigLppnok=repmat(Lppnok,1,length(k));


matrices.Lpp=-kjvec.*BigLppnok;
matrices.Lpp=repmat(matrices.Lpp,1,2*params.nx/c);
end