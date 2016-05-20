function [normal_high] = normhigh(M)
%[normal_high] = normhigh(M)
%normalizes to highest 
%takes matrix M and makes normal_high
%
% Copyright (c) 1996, Royston Goodacre
%

[maxheight,where]=max(M');
[rows,cols]=size(M);
normal_high=zeros(rows,cols);
for i=1:rows
   normal_high(i,:)=M(i,:)/maxheight(i);
end
