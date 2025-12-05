function [uhat,phat,M]=solveuhatphat(Luu, Lup, Lpu, Lpp, b,n)
        M=zeros(2*n+2);
        %b=zeros(2*n+2,1);
        for k=1:n
    
            M(k:k+1,k:k+1) = M(k:k+1,k:k+1)+ (-1)^(n-k) * Luu(1:2,2*k-1:2*k);
            M(k:k+1,k+n+1:k+n+2) = M(k:k+1,k+n+1:k+n+2)+(-1)^(n-k) * Lup(1:2,2*k-1:2*k);
            M(k+n+1:k+n+2,k:k+1) = M(k+n+1:k+n+2,k:k+1)+(-1)^(n-k) * Lpu(1:2,2*k-1:2*k);
            M(k+n+1:k+n+2,k+n+1:k+n+2) = ...
                M(k+n+1:k+n+2,k+n+1:k+n+2)+(-1)^(n-k) * Lpp(1:2,2*k-1:2*k);
            %b(k:k+1) = b(k:k+1)+(-1)^(n-k+1)*bu(1:2,k);
            %b(k+n+1:k+n+2,1)= b(k+n+1:k+n+2,1)+(-1)^(n-k+1)*bp(1:2,k);
        end
        
       
        %Make sure the first u is 0 and the first p is zero
        M(1,:)=0;
        M(1,1)=1;
        M(n+2,:)=0;
        M(n+2,n+2)=1;
        nv=M\b;
        uhat=nv(1:n+1,1); 
        phat=nv(n+2:2*n+2,1);
end