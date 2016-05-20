 function [Mvecnorm] = vecnorm(M);
%[Mvecnorm] = vecnorm(M);
%
%computer vector normalisation a la Bruker s/w
%well almost
%
% Copyright (c) 2000, Royston Goodacre
%

%size matrix first
   [n,m]=size(M);

%average of y-value is first calculated:
   mean_y=mean(M');

%This y-value is then subtracted from the 
%spectrum from which it was calculated in
%such that the spectrum is centered about y=0
   for i = 1:n
     Z(i,:) = M(i,:)-mean_y(i);
   end;


%the sum of the squared y-values are then calculated
Z2=Z.*Z;
Z2sum=sum(Z2');

%the square root of this sum is calculated
Z2sum2half=(Z2sum.^0.5)';

%the spectrum is then divided by the square root of this sum
   for i = 1:n
     Mvecnorm(i,:) = M(i,:)/Z2sum2half(i,:);
   end;

%the Euclidean vector norm of the new spectrum is 1
%I have bidged this bit...
Mvnorm=scalem(Mvecnorm);

