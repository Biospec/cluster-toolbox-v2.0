function U = explv(X,T,P)
% U = explv(X,T,P)
% X = raw data matrix (before subtraction of mean)
% T = score matrix
% P = transposed loadings matrix.

[n,m] = size(X);
[n,comp] = size(T);

E0 = X - ones(n,1)*mean(X);
s0 = sum(sum(E0.^2));

for i = 1:comp
 E = E0-T(:,1:i)*P(1:i,:);
 s(i) = sum(sum(E.^2)); 
end;

U = (1-(s/s0))*100;

