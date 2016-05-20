function [C,T]=plspred(X,P,Q,W,b,N,amean,ascal,cmean, cscal)
% function [C,t]=plspred(X,P,Q,W,b,N,amean,ascal,cmean, cscal)
%X: data matrix for PLS predction
%P,Q,W,b have the same meaning as in PLS
%N: the number of components for prediction
%amean (optional): spectral matrix average of training set
%ascale (optional): scale coefficent of spectral matrix in training set,e.g., stand deviation
%cmean (optional): average of concentraion matrix in training set
%cscale (optional): scale coefficent of concentration matrix in training set,e.g., stand deviation
%C: predicted concentrations
%T: projected pls scores of X

if nargin<5
    help plspred;
    return;
elseif nargin==5
    N=size(P,2);
elseif nargin>=8
    [m,n]=size(X);
    X=X-ones(m,1)*amean;
    X=X./(ones(m,1)*ascal);
end
[m,n]=size(X);
c=zeros(m,size(Q,1));
for i=1:N
    T(:,i)=X*W(:,i);
    c=c+b(i)*T(:,i)*Q(:,i)';
    X=X-T(:,i)*P(:,i)';
end
C=c;
if nargin==10
    
    C=C.*(ones(m,1)*cscal);
    C=C+ones(m,1)*cmean;
end
    
    