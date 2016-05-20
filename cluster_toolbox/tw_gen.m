function [T,W] = tw_gen(X,group)
% [T,W] = TW_gen(X,group)
% Generates the T and W matrices used in DISCRIMINANT FUNCTION ANALYSIS
% X contains m groups 
% group is a vector containing a number corresponding
% to a group for every row in X. If you have m groups
% there will be numbers in the range 1:m in this vector.
%
% Copyright, B.K. Alsberg, 1996
%

un = unique(group);
no_groups = length(un);
% The use of unique is to make the routine more general.
% In this case we may have different names on classes that are
% not necessarily in successive order. With the old program
% we could not handle class vector having members like [1 5 21]
% it had to be [1 2 3].


[n,m]=size(X);

mx = mean(X);

% Making T and W

for j = 1:no_groups
  idx = find(group==un(j));
  L = length(idx);
  K = X(idx,:);
  [nn,mm]=size(K);
  if nn > 1,
     zz = mean(K);
  else
     % In the special case of one single member of a class the
     % mean function in MATLAB will return a single number for a vector
     zz = K;
  end;
  A = K - ones(L,1)*mx;
  C = K - ones(L,1)*zz;
  %C = K - ones(L,1)*mean(K);
  if j > 1
    T = T + A'*A;
    W = W + C'*C;
  elseif j==1
    T = A'*A;
    W = C'*C;
  end;
end;

%xm = mean(X);
%T =zeros(m,m);
%W =T;

%for i = 1:no_groups
%  idxi = find(group==i);
%  Xi = X(idxi,:);
%  xmi = mean(Xi);
%  for j = 1:no_groups
%    idxj = find(group==j);
%    Xj = X(idxj,:);
%    xmj = mean(Xj);
%    T1 = Xj-ones(size(Xj))*xm(j);
%    T2 = Xi-ones(size(Xi))*xm(i);
%    T(i,j) = sum(sum(T1.*T2));
%    W1 = Xj-ones(size(Xj))*xmj(j);
%    W2 = Xi-ones(size(Xi))*xmi(i);
%    W(i,j) = sum(sum(W1.*W2));
%  end;
%end;


