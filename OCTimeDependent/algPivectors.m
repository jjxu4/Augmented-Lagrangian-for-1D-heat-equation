function [matrices]=algPivectors(params,setup,matrices,source,endpoints,t,Pi,a,hi)
% Defines biplus1, liplus1, ri, fiplus1nok
    x_input=.5*params.x(1:end-1)+.5*params.x(2:end);

if (-1)^setup==1 
    if isa(source.F,'function_handle')==1
        Fvec=source.F(x_input,t);
    else        
        %If F is a matrix the first column should include the initial
        %condition.
        %if nargin>7
            Fvec=source.F(:,round((t-endpoints.tstart)/params.deltat));
            [a,b]=size(Fvec);
            if a==params.n+1 || b==params.n+1
                Fvec=.5*Fvec(1:end-1)+.5*Fvec(2:end);
            end
        %else
        %    Fvec=source.F(:,round((t-endpoints.tstart)/params.deltat));
        %end
    end
    
    if isa(source.S,'function_handle')==1
        Svec1=source.S(x_input,t);
    else
        Svec1=source.S(:,round((t-endpoints.tstart)/params.deltat)+1);
            [a,b]=size(Svec1);
            if a==params.n+1 || b==params.n+1
                Svec1=.5*Svec1(1:end-1)+.5*Svec1(2:end);
            end
    end
    Fvec=reshape(Fvec,params.n,1);
    
    Svec=reshape(Svec1,params.n,1);
    matrices.biplus1=(params.h*Fvec);
    liplus1=(params.h*Svec)';
    if nargin==7
        matrices.liplus1=liplus1-(params.h/(params.lambdae*params.deltat)*Pi);
    elseif nargin~=4 && nargin~=6
        matrices.liplus1=-liplus1;
    end
    matrices.ri=params.delta*params.Hv/(params.deltat*params.lambdae)*matrices.Ai*matrices.D*Pi;
    
    if nargin>7
        hi=reshape(hi,1,[]);
        matrices.ri=-matrices.ri+1/params.deltat*matrices.Ai*matrices.D*hi;
    end

    fiplus1nok=matrices.biplus1;
    if nargin<8
        matrices.fiplus1nok=fiplus1nok+transpose(matrices.B*matrices.ri);
    else
        matrices.fiplus1nok=fiplus1nok-transpose(matrices.B*matrices.ri);
    end
else
    if isa(source.F,'function_handle')==1
        Fvec=source.F(x_input);
    else
        Fvec=source.F(1:end-1)/2+source.F(2:end)/2;
    end
    
    if isa(source.S,'function_handle')==1
        Svec=source.S(x_input);
    else
        Svec=source.S(1:end-1)/2+source.S(2:end)/2;
    end

    l=numel(Fvec);
    Fvec=reshape(Fvec,l,1);
    
    Svec=reshape(Svec,l,1);
    matrices.biplus1=(params.h*Fvec);
    matrices.liplus1=(params.h*Svec)';
    if (setup==1 && nargin>5)
        matrices.liplus1=-matrices.liplus1;
    end
    matrices.fiplus1nok=matrices.biplus1;
    matrices.ri=zeros(2,length(matrices.fiplus1nok));
end
end