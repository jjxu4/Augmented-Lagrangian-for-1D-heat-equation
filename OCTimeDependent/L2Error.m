function [erroru,errorp] = L2Error(Exactu,Exactp,Approxu, Approxp,endpoints)
%Gives L2 Error 5 arguments does it for 2 functions. 3 arguments does it just
%for 1 function
%   Approx u and Approx p are step functions
if nargin==3
    endpoints=Approxu;
    Approxu=Exactp;
end
    [xnodes,tnodes]=size(Approxu);
    erroru=0;
    errorp=0;

    if tnodes==1
        factor=1;
    else 
        factor=(endpoints.tend-endpoints.tstart)/tnodes;
    end

    x=linspace(endpoints.xstart,endpoints.xend,xnodes);
    for k=1:tnodes
        for j=1:xnodes-1
            if tnodes==1
                difference_u=@(y)(Exactu(y)-Approxu(j+1,k)).^2;
                if nargin==5
                    difference_p=@(y)(Exactp(y)-Approxp(j,k)).^2;
                end
            else
                difference_u=@(y)(Exactu(y,endpoints.tstart+factor*k)-Approxu(j+1,k)).^2;
                if nargin==5
                    difference_p=@(y)(Exactp(y,endpoints.tstart+factor*k)-Approxp(j,k)).^2;
                end
            end
            erroruj=integral(difference_u,x(j),x(j+1));
            erroru=erroru+erroruj*factor;
            if nargin==5
                errorpj=integral(difference_p,x(j),x(j+1));
                errorp=errorp+errorpj*factor;
            end
        end
    end
    erroru=sqrt(erroru);
    if nargin==5
        errorp=sqrt(errorp);
    end
end

