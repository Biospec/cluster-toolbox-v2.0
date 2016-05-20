function [C,BV]=plspred2(X,P,Q,W,b,N,amean,cmean,ascal,cscal)
%function [C,BV]=plspred2(X,P,Q,W,b,n,amean, cmean, ascal, cscal);
%Same to plspred. But quick prediction procedure is used.
%please refer to plspred
%BV: regression coefficents

[m,n]=size(X);
if nargin<5
    help plspred;
    return;
elseif nargin==5
    N=size(P,2);
elseif nargin==8
%     [m,n]=size(X);
    X=X-ones(m,1)*amean;
elseif nargin>8
    X=X-ones(m,1)*amean;
    X=X./(ones(m,1)*ascal);
end
W=W(:,1:min(N, size(W,2)));
P=P(:,1:min(N, size(W,2)));
b=b(1:min(N, size(W,2)));
Q=Q(:,1:min(N, size(W,2)));
R=Q*diag(b);

B=W*inv(P'*W)*R';
C=X*B;
BV=B;


%  [m,n]=size(X);

if nargin==8
    C=C+ones(m,1)*cmean;
end   
if nargin==10
    C=C+ones(m,1)*cmean;
    C=C.*(ones(m,1)*cscal);
end
end