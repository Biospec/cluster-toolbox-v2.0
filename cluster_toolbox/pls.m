function [T,P,Q,W,b,X,Y]=pls(A,C,maxrank)
% function [T,P,Q,W,b,X,Y]=pls(A,C,maxrank);
%A: spectral matrix in training set
%C: concentration matrix in training set
%maxrank: the number of components to be extracted.
%T: scores of spectral matrix
%P: loadings of spectral matrix
%Q: loadings of concentration matrix
%W: weight matrix
%b: coefficient
%X: Residual spectral matrix
%Y: Residual concentration matrix
%This program is for general use, eithe PLS1 or PLS2.


[m,n]=size(A);
maxrank=min([m n maxrank]);

X=A;
% X=[A; diag(ones(n,1)*10)];
Y=C;
% Y=[Y; zeros(n,size(C,2))];
for i=1:maxrank
    XY=X'*Y;
    [u,s,v]=svd(XY,0);
    q=v(:,1);
    w=u(:,1);
    t=X*w;
    p=t'*X/(t'*t);p=p';
    u=Y*q;
    b(i,1)=1/(t'*t)*t'*u;
    X=X-t*p';
    Y=Y-b(i)*t*q';
    
    T(:,i)=t;
    W(:,i)=w;
    Q(:,i)=q;
    P(:,i)=p;
    
end
end