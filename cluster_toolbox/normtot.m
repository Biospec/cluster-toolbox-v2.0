function [normal_tot] = normtot(M)
%[normal_tot] = normtot(M)
%normalizes to total 
%takes matrix M and makes normal_tot
%
% Copyright (c) 1996, Royston Goodacre
%

rowsum=sum(M');
normal_tot=(M'*inv(diag(rowsum)))';

