function [U,V,eigenvals] = dfa(X,group,maxfac)
%[U,V,eigenvals] = DFA(X,group,maxfac)
% Performs DISCRIMINANT FUNCTION ANALYSIS
%
% INPUT VARIABLES
%
% X              = data matrix that contains m groups
%                  Dim(X) = [N x M]. All columns must be independent.
% group          = a vector containing a number corresponding
%                  to a group for every row in X. If you have 
%                  m groups there will be numbers in the range 
%                  1:m in this vector.
% maxfac         = the maximum number of DFA factors extracted
%
% OUTPUT VARIABLES
%
% U              = DFA scores matrix (Dim(U) = [N x maxfac])
%                  the eigenvalues are multiplied with each column
%                  in this matrix.
% V              = DFA loadings matrix, Dim(V) = [M x maxfac]
% eigenvals      = a vector of DFA eigenvalues
%
%
% Copyright, B.K. Alsberg, 1996
% Updated, Yun Xu, 2016

no_groups = length(unique(group));
% Sanity check on maxfac
if maxfac > no_groups - 1
    warning('Too many DFs to be extracted, set it to group-1')
    maxfac = no_groups - 1;
end
[T,W]=tw_gen(X,group);

B = T-W;

invW = inv(W);
P = invW*B;
% replace eig with eigs for automated a selected few eigenvalues extraction
% and no longer need to sort after decomposition, by Yun
% [vec1,val]=eig(P);
[vec,val]=eigs(P, maxfac);
% d=(diag(val))';
% eigenvals = d(1:maxfac);
eigenvals = diag(val);
[eigenvals, sort_idx] = sort(eigenvals, 'descend');
vec=vec(:,sort_idx);
% sanity check of the solution
if any(~isreal(eigenvals))
    error('complex eigenvalues detected, maxfac is probably too high')
end
% Here we sort the eigenvectors w.r.t. the eigenvalues:
% [dummy,idx]=sort(-eigenvals);
% vec = vec1(:,idx);
% eigenvals = eigenvals(idx);

%% V is the matrix of canonical variates directions %%%
V = vec;

%% U is the matrix of scores %%%
U = X*V;

% new line to multiply eigenvalues to DFA directions:
% This is not right, U should already carry the variance, no need to
% multiply eigenvals again. by Yun
% U = U*diag(eigenvals);
% we need a reference to see how this is done properly


% U = real(U);

