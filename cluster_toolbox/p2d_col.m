function [] = p2d_col(X,col1,col2,nam,cc)
%function [] = p2d_col(X,col1,col2,nam,cc)
% numplot plots a 2D plot with col1 against col2
% nam is a matrix [rows x q] where q is the max allowed
% numbers of characters for each name
% cc is 'colour' to plot, white (i.e. not shown) if not given
%
% Copyright (c) 1998, Royston Goodacre
% Updated 2016 by Yun Xu

if nargin == 4
 cc = 'w';
end;


[n,m] = size(X);

plot(X(:,col1),X(:,col2),'.','color',cc);

for i = 1:n
 h=text(X(i,col1),X(i,col2),nam(i,:));
 if ~strcmp(cc, 'w')
    set(h,'color',cc);
 end
end;

