function [MF, MS]=mos(x)

% [MF, MS]=mos(x)
% MF: Morphological factor
% MS: Morphological Score
% x: a vector of signals

diff_x=diff(x);
sign_x=sign(diff_x);
zcp=length(find(abs(diff(sign_x))==2));
MF=norm(x)/(norm(diff_x)*zcp);
x_m=x-mean(x);
MS=norm(x_m)/(norm(diff(x_m)));