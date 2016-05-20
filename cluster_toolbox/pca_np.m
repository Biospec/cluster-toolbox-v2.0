function [tt,pp,pr] = PCA(X,comp)
% [tt,pp,pr]=PCA(X,comp);
% comp is the number of principal components
% Here the NIPALS algorithm is used
% pr = percentage explained variance 

X0 = X;
m = mean(X);
[a,b] = size(X);
X =  X - ones(a,1)*m;
% Centering the matrix
% the first approximation to t is the first
% column in X

stop = 3;
for i = 1:comp
i
av = std(X);

[av2,indx]=sort(av);
    t0 = X(:,indx(b));
       while stop > 1,
          p0 = t0'*X;
          p1 = p0/norm(p0);
          t1 = X*p1';
          if norm(t1-t0) < 0.00005,
             tt(:,i) = t1;
             pp(i,:) = p1;
             stop = 0;
          end;
          t0 = t1;
       end;
     stop = 3;
     X = X - t1*p1;
  end;


pr = (explv(X0,tt,pp))';

%eig = diag(tt'*tt);
%ssum = cumsum(eig);
%pr   = 100*ssum/ssum(comp);









