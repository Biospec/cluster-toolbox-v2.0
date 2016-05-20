function [basl] = detrendm(M)
%[basl] = detrendM(M)
%detrends from first and last bin
%takes matrix M and makes basl
%
% Copyright (c) 1997, Royston Goodacre
%


[rows,cols]=size(M);
basl=zeros(rows,cols);

for i=1:rows
   first=M(i,1);
   last=M(i,cols);
   range=M(i,cols)-M(i,1);
   slide = lintrans(1:cols,1,cols,first,last);
   basl(i,:)=M(i,:)-slide;
end

