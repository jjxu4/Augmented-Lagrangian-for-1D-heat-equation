function [mhat,hhat,M,b,matrices]=constlinPVE(source,params,a)


if nargin<3
    [matrices,params]=algmatrices(params,1);
    matrices=algkmatrices(params,1,matrices,params.kref);
    matrices=algPivectors(params,1,matrices, source);
    matrices=bmatrices(params,matrices);
    [mhat, hhat, M, b]=solveuhatphat(matrices.Luu, matrices.Lup,...
        matrices.Lpu, matrices.Lpp, bu,bp, source, params.n,matrices,params);
else
    [matrices,params]=algmatrices(params,1,a);
    matrices=algkmatrices(params,1,matrices,params.kref,a);
    matrices=algPivectors(params,1,matrices, source,a);
    [bm,bh,matrices]=adjbvectors(params,matrices,params.kref,1);
    [mhat, hhat, M, b]=solveuhatphat(matrices.Luu, matrices.Lup,...
        -matrices.Lhm, -matrices.Lpp, bm,bh, source, params.n,matrices,params);
end
[u,p]=solveup(mhat, hhat, matrices);
end