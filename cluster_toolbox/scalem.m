function [basl] = scaleM(M)
%[basl] = scaleM(M)
%scales so min and max are between 0 and 1
%takes matrix M and makes basl
%
% Copyright (c) 1997, Royston Goodacre
%


[rows,cols]=size(M);
basl=zeros(rows,cols);

for i=1:rows
   Mmin=min(M(i,:));
   Mmax=max(M(i,:));
   for j=1:cols
      basl(i,j)=(M(i,j)-Mmin)/(Mmax-Mmin);
   end
end


