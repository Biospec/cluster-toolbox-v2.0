function [ut] = lintrans(x,x1,x2,y1,y2)
% Makes a linear transformation from
% (x1,y1) ---> (x2,y2). E.g. if x = 1:4 is
% to be mapped into xx = -5:5 we have to following usage:
% ut = lintrans(x,1,4,-5,5);
a = (y1-y2)/(x1-x2);
b = y1 - a*x1;
ut = a*x + b;


