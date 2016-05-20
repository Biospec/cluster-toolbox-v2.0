function [] = plot_dfa(U,DF_x,DF_y,names)
%plot_dfa(U,DF_x,DF_y,names)
%Plot DFA/CVA with labels, title, etc
%where U = DF scores matrix
%DF_x and DF_y are the two ordinates to be plotted
%names contains integer identifier file
%
% Copyright (c) 1997, Royston Goodacre
%

%asks for title
tit = input('Title please (in quotes)  ');

figure

%turns number to integers for labels
DF_X1=num2int(DF_x);
DF_Y1=num2int(DF_y);
DF_X2=['Discriminant function ',DF_X1];
DF_Y2=['Discriminant function ',DF_Y1];

%plot with axis labels
plotnm2(U,DF_x,DF_y,names);
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',[14])
%hold on
xlabel(DF_X2)
ylabel(DF_Y2)
title(tit,'FontSize',[16],'FontWeight','bold')
%hold off

