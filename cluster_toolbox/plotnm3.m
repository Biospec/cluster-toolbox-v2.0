function [] = plotnm3(X,col1,col2,col3,nam)
%function [] = plotnm3(X,col1,col2,col3,nam)
% numplot plots a 3D plot with numbers from 
% 1 to number of rows
% nam is a matrix [rows x q] where q is the max allowed
% numbers of characters for each name
%
% Copyright (c) 1997, Royston Goodacre
%

if nargin == 4
 cc = 'w';
end;
 
[n,m] = size(X);
 plot3(X(:,col1),X(:,col2),X(:,col3),'.w');
for i = 1:n
 h=text(X(i,col1),X(i,col2),X(i,col3),nam(i,:));
% set(h,'color',cc);
 if nargin==6
   set(h,'fontsize',marksiz);
 end;
end;

