function K = polyker(X, Y, p, b)
% K = ploy(X, Y, p, b)
% compute polynomial kernel function
% X, Y are the data matrices, if only X is given then Y = X
% p: the oder of polynomial, set to 1, i.e. linear kernel if ignored
% b: the bias term, set to 0 if ignored


if nargin<2 || isempty(Y)
    Y = X;
end
if nargin<3 || isempty(p)
    % linear kernel if p is ignored
    p = 1;
end
if nargin<4
    % no bias term if b is ignored
    b = 0;
end

K = (Y*X'+b).^p;
