function b=derive_b(matrices,source,uhat,x,t,setup,endpoints,params,BC)

if (-1)^setup==1 
    if isa(source.F,'function_handle')==1
        Fvec=source.F(x,t);
        Fvec=Fvec';
    else        
        %If F is a matrix the first column should include the initial
        %condition since uhat imcludes initial condition
        Fvec=source.F(:,round((t-endpoints.tstart)/params.deltat));
    end
    
    if isa(source.S,'function_handle')==1
        Svec=source.S(x,t);
        Svec=Svec';
    else
        Svec=source.S(:,round((t-endpoints.tstart)/params.deltat)+1);
    end

else
    if isa(source.F,'function_handle')==1
        Fvec=source.F(x);
        Fvec=Fvec';
    else
        Fvec=source.F;
    end
    
    if isa(source.S,'function_handle')==1
        Svec=source.S(x);
        Svec=Svec';
    else
        Svec=source.S;
    end
end

y=zeros(params.n+1,1);
y=[uhat;y];

b=matrices.Bcal*y+matrices.BcalS*Svec+matrices.BcalF*Fvec;

%Add in g and psi to match 197 and 199.

if (-1)^setup==1 
    b(params.n+1,1)=b(params.n+1,1)+BC.g;
    b(2*params.n+2,1)=b(2*params.n+2,1)+BC.psi;
else
    b(params.n+1,1)=b(params.n+1,1)+source.g;
    b(2*params.n+2,1)=b(2*params.n+2,1)+source.psi;
end

b(1,1)=source.bcu(0);
b(params.n+2,1)=source.bcp(0);


end