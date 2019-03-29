function [pval, D_ob, D_null] = protest(mat1, mat2, n)
% [pval, D_ob, D_null] = protest(mat1, mat2)
% This function assess the statistical significance of the similarity
% between two matrices measured by Procrustes Analysis.
% mat1 is PCA scores matrix for data set 1. 
% mat2 is the scores matrix for data set 2.
% pval is the p-value of the statistical significance of Procrustes
% distance estimated by permutation test
% D_ob is the observed Procrustes distance
% D_null is the NULL distribution of Procrustes distances of permuted data.

[m1, n1] = size(mat1);
[m2, n2] = size(mat2);

if ~exist('n', 'var')
    n = 1000;
end

if m1 ~= m2
    error('The number of rows of the two matrices to be compared must be the same!')
end

if n1 < n2
    error('The number of variables in mat2 cannot be greater than those in mat1')
end

D_ob = procrustes(mat1, mat2); 
D_null = zeros(n, 1);
for i = 1:n % 1000 permutations, i.e. n = 1000
    D_null(i) = procrustes(mat1, mat2(randperm(size(mat2,1)), :));
end

pval = numel(find(D_null < D_ob)) / numel(D_null);