function [u,p]=solveup(uhat, phat,matrices)
    phatasmatrix=[phat(1:length(phat)-1,1)';phat(2:length(phat),1)'];
    uhatasmatrix=[uhat(1:length(uhat)-1,1)';uhat(2:length(uhat),1)'];
    reshapeuhat=reshape(uhatasmatrix,1,numel(uhatasmatrix));
    [~,c]=size(matrices.RR);
    BigRR=repmat(matrices.RR,1,numel(uhatasmatrix)/(c));
    RRuhatfull=BigRR.*reshapeuhat;
    RRuhateven=RRuhatfull;
    RRuhateven(1:2:end) = [];
    RRuhatodd=RRuhatfull;
    RRuhatodd(2:2:end)=[];
    RRuhat=RRuhateven+RRuhatodd;
    QQphat=matrices.QQ.*phatasmatrix;
    QQphat=QQphat(1,:)+QQphat(2,:);   
    p=matrices.Q*phatasmatrix+sum(matrices.R'.*uhatasmatrix)...
        +(matrices.scriptBpinv'.*matrices.liplus1);
    p=p';
    u=RRuhat'+QQphat'+matrices.scriptBuinv*matrices.fiplus1';

    
end