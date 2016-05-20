function [] = plot_pca(tt,PC_x,PC_y,names)
%plot_pca(tt,PC_x,PC_y,names)
%Plot PCA with labels, title, etc
%where tt = PC scores matrix
%PC_x and PC_y are the two ordinates to be plotted
%names contains integer identifier file
%
% Copyright (c) 1996, Royston Goodacre
%

%asks for title
tit = input('Title please (in quotes)  ');

figure

%turns number to integers for labels
PC_X1=num2int(PC_x);
PC_Y1=num2int(PC_y);
PC_X2=['Principal component ',PC_X1];
PC_Y2=['Principal component ',PC_Y1];

%plot with axis labels
plotnm2(tt,PC_x,PC_y,names);
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',[14])
hold on
xlabel(PC_X2)
ylabel(PC_Y2)
title(tit,'FontSize',[16],'FontWeight','bold')
hold off

