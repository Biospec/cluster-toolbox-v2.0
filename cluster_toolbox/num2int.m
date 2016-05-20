function [ints]=num2int(M)
% [ints]=num2int(M)
%converts numerical matrix M into integer matrix
%
% Copyright (c) 1996, Royston Goodacre
%

[Mrows,Mcols] = size(M);

ints=[];

for i = 1:Mrows
 hh=sprintf('%-5.0f',M(i,1));
 ints=[ints;hh];
end
