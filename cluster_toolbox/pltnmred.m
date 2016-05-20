function [] = plotnm2(X,col1,col2,nam)
%function [] = plotnm2(X,col1,col2,nam)
% numplot plots a 2D plot with numbers from 
% 1 to number of rows
% nam is a matrix [rows x q] where q is the max allowed
% numbers of characters for each name
%
% Copyright (c) 1997, Royston Goodacre
%

if nargin == 4
 cc = 'r';
end;
 
[n,m] = size(X);
 plot(X(:,col1),X(:,col2),'.w');
for i = 1:n
 h=text(X(i,col1),X(i,col2),nam(i,:));
 set(h,'color',cc);
% set(h,'color',cc);
% if nargin==6
% end;

%marksiz=8;
%   set(h,'fontsize',marksiz);


end;

