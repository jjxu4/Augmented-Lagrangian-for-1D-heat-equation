function matrices=bmatrices(params,matrices)
    matrices.ri_no_u=-params.delta*params.Hv/params.deltat*matrices.Ai*matrices.D*matrices.G;
    matrices.lu=params.h/params.deltat*matrices.G;
    if size(params.kref)==[1,1]
        matrices.fu=matrices.B*matrices.ri_no_u-matrices.B*matrices.Ai*matrices.D*matrices.lu*matrices.scriptBpinv;
    else
        a=matrices.B*matrices.ri_no_u;
        b=-matrices.B*matrices.Ai*matrices.D*matrices.lu;
        a=repmat(a,12,1);
        b=repmat(b,12,1);
        doubBpinv=repmat(matrices.scriptBpinv,1,2);
        b=b.*doubBpinv;
        matrices.fu=a+b;
    end
    
    halfh=[params.h/2 params.h/2];
    
    matrices.fqS=-matrices.B*matrices.Ai*matrices.D*matrices.scriptBpinv*halfh;
    if size(params.kref)==[1,1]
        matrices.buu=matrices.ri_no_u-params.mu*matrices.Ai*matrices.B'*matrices.scriptBuinv*matrices.fu...
        -matrices.Ai*matrices.D*matrices.scriptBpinv*matrices.lu;
    else
        fu=[matrices.fu,matrices.fu];
        fu=reshape(fu,2*params.nx,2);
        multBpinv=repmat(matrices.scriptBpinv,1,4);
        multBpinv=reshape(multBpinv,2*params.nx,2);
        a=matrices.ri_no_u;
        b=-params.mu*matrices.Ai*matrices.B'*matrices.scriptBuinv;
        c=-matrices.Ai*matrices.D*matrices.lu;
        b=repmat(b,12,2);
        A=repmat(a,12,1);
        c=repmat(c,12,1);
        B=b.*fu;
        C=multBpinv.*c;
       
        matrices.buu=A+B+C;
    end
        matrices.bpu=1/matrices.scriptB*matrices.Ai*matrices.B'*matrices.lu;

%     matrices.betauS=params.mu*matrices.Ai*matrices.B'*matrices.scriptBuinv*matrices.B*matrices.Ai*matrices.D*matrices.scriptBpinv*[params.h/2 params.h/2]...
%         -matrices.Ai*matrices.D*matrices.scriptBpinv*[params.h/2 params.h/2];
%     matrices.betapS=params.kref*matrices.Ai*matrices.B'*matrices.scriptBpinv*[params.h/2 params.h/2];
    if size(params.kref)==[1,1]
        matrices.betauS=-params.mu*matrices.Ai*matrices.B'*matrices.scriptBuinv*matrices.fqS-matrices.Ai*matrices.D*matrices.scriptBpinv*halfh;
    else
        doubfqS=repmat(matrices.fqS,1,2);
        doubfqS=reshape(doubfqS,2*params.nx,2);
        a=-params.mu*matrices.Ai*matrices.B'*matrices.scriptBuinv;
        a=repmat(a,params.nx,2);
        b=-matrices.Ai*matrices.D*halfh;
        b=repmat(b,params.nx,1);
        matrices.betauS=a.*doubfqS+b.*multBpinv;
    end
    matrices.betapS=matrices.Ai*matrices.B'*1/matrices.scriptB*halfh;
    matrices.betauF=-params.mu*matrices.Ai*matrices.B'*matrices.scriptBuinv*[params.h/2 params.h/2];
    
    matrices.Bcal=zeros(2*params.nx+2);
    matrices.BcalS=zeros(2*params.nx+2,params.nx+1);
    matrices.BcalF=zeros(2*params.nx+2,params.nx+1);
    
    offset = params.nx + 1;
    for ix = 1:ceil(params.nx / 2)
        if 2*ix <= params.nx + 1
            if size(params.kref)==[1,1]
                matrices.Bcal(2*ix - 1:2*ix, 2*ix - 1:2*ix) = (-1)^params.nx*matrices.buu;
                matrices.BcalS(2*ix - 1:2*ix, 2*ix - 1:2*ix) = (-1)^params.nx*matrices.betauS;
            else
                matrices.Bcal(2*ix - 1:2*ix, 2*ix - 1:2*ix) = (-1)^params.nx*matrices.buu(4*ix-3:4*ix-2,1:2);
                matrices.BcalS(2*ix - 1:2*ix, 2*ix - 1:2*ix) = (-1)^params.nx*matrices.betauS(4*ix-3:4*ix-2,1:2);
            end
            matrices.Bcal(offset + 2*ix - 1:offset + 2*ix, 2*ix - 1:2*ix) = (-1)^params.nx*matrices.bpu;
            

            matrices.BcalS(offset + 2*ix - 1:offset + 2*ix, 2*ix - 1:2*ix) = (-1)^params.nx*matrices.betapS;
            
            matrices.BcalF(2*ix - 1:2*ix, 2*ix - 1:2*ix) = (-1)^params.nx*matrices.betauF;
        end
    end
    
    for ix = 1:ceil(params.nx / 2)
        if (2*ix + 1 <= params.nx + 1)
            if size(params.kref)==[1,1]
                matrices.Bcal(2*ix:2*ix + 1, 2*ix:2*ix + 1) =matrices.Bcal(2*ix:2*ix + 1, 2*ix:2*ix + 1)- (-1)^params.nx*matrices.buu(1:2,1:2);
                matrices.BcalS(2*ix:2*ix + 1, 2*ix:2*ix + 1) =matrices.BcalS(2*ix:2*ix + 1, 2*ix:2*ix + 1)- (-1)^params.nx*matrices.betauS(1:2,1:2);
            else
                matrices.Bcal(2*ix:2*ix + 1, 2*ix:2*ix + 1) =matrices.Bcal(2*ix:2*ix + 1, 2*ix:2*ix + 1)- (-1)^params.nx*matrices.buu(4*ix-1:4*ix,1:2);
                matrices.BcalS(2*ix:2*ix + 1, 2*ix:2*ix + 1) =matrices.BcalS(2*ix:2*ix + 1, 2*ix:2*ix + 1)- (-1)^params.nx*matrices.betauS(4*ix-1:4*ix,1:2);
            end
            matrices.Bcal(offset + 2*ix:offset + 2*ix + 1, 2*ix:2*ix + 1) = matrices.Bcal(offset + 2*ix:offset + 2*ix + 1, 2*ix:2*ix + 1)-(-1)^params.nx*matrices.bpu;           
  
            matrices.BcalS(offset + 2*ix:offset + 2*ix + 1, 2*ix:2*ix + 1) = matrices.BcalS(offset + 2*ix:offset + 2*ix + 1, 2*ix:2*ix + 1)-(-1)^params.nx*matrices.betapS(1:2,1:2);
            
            matrices.BcalF(2*ix:2*ix + 1, 2*ix:2*ix + 1) =matrices.BcalF(2*ix:2*ix + 1, 2*ix:2*ix + 1)- (-1)^params.nx*matrices.betauF(1:2,1:2);

        end
    end
    matrices.Bcal(1,:)=zeros(1,2*params.nx+2);
    matrices.Bcal(params.nx+2,:)=zeros(1,2*params.nx+2);
    matrices.BcalS(1,:)=zeros(1,params.nx+1);
    matrices.BcalS(params.nx+2,:)=zeros(1,params.nx+1);
    matrices.BcalF(1,:)=zeros(1,params.nx+1);
    matrices.BcalF(params.nx+2,:)=zeros(1,params.nx+1);
end