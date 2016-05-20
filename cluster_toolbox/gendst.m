function [D] = gendst(X,type,param)
% [D] = gendst(X,type,[param])
% This is a general distance calculation procedure
% type = {Mink, Cheb}.
% param = 1,2,3...
% Note that {Mink,1}   = Manhattan (City block) distance
%           {Mink,2}   = ordinary Euclidean distance
%           {Cheb,*}   = Minkowski for param = infinity
%           {Maha,*} = Mahalanobis distance using the total
%                             covariance matrix of X
%           {Canb,*}    = "Canberra" metric
%           {Angu,*}     = Angular separation
%            * = means there are no parameters 
%
% Copyright (c) 1996, B.K. Alsberg
%

[n,m] = size(X);

% This speeds up the calculation
D = zeros(n,n);

% A more extensive list of distance measures can be found in:
% Gower, J.C. (1985) Measures of similarity, dissimilarity and
% distance. In Encyclopedia of Statistical Sciences, Volume 5
% (S. Kotz, N.L. Johnson and C.B. Read, eds.) John Wiley & Sons,
% New York


if type == 'Cheb'

% Calculating the Chebychev distance matrix

 for i = 1:n
  for j = 1:n
   td = abs(X(i,:) - X(j,:));
   D(i,j) = max(td);
  end
 end

elseif type == 'Mink'

% Calculating the Minkowski distance matrix
 for i = 1:n
  for j = 1:n
   p = (abs(X(i,:)-X(j,:))).^param;
   D(i,j) = (sum(p)).^(1/param);
  end
 end

elseif type == 'Maha'

% We must have a check here that
% uses SVD for singular matrices and
% uses the k first eigenvectors

 C = X'*X;
 S = inv(C);

 for i = 1:n
  for j = 1:n
   v = X(i,:)-X(j,:);
   D(i,j) = v*S*v';
  end
 end

elseif type == 'Canb'

 for i = 1:n
  for j = 1:n
   td = abs(X(i,:) - X(j,:))./(X(i,:) + X(j,:));
   D(i,j) = sum(td);
  end
 end

elseif type == 'Angu'

 for i = 1:n
  for j = 1:n
   a = X(i,:).*X(j,:);
   b = X(i,:).^2;
   c = X(j,:).^2;
   D(i,j) = sum(a)/sqrt(sum(b)*sum(c));
  end
 end

else
 error('Not a valid distance measure!');
end;
