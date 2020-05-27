function [normal_tot] = normtot(M)
%[normal_tot] = normtot(M)
%normalizes to total 
%takes matrix M and makes normal_tot
%
% Copyright (c) 1996, Royston Goodacre
% Updated 2020, Yun Xu

% rowsum=sum(M');
% normal_tot=(M'*inv(diag(rowsum)))';

% New code
rowsum = sum(M, 2);
normal_tot = M ./ repmat(rowsum, 1, size(M, 2));

