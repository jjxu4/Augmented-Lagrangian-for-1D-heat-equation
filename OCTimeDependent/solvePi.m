function [Pi]=solvePi(uhat,lambdae,G)
    uhatasmatrix=[uhat(1:length(uhat)-1,1)';uhat(2:length(uhat),1)'];
    Pi=-lambdae*G*uhatasmatrix;
end