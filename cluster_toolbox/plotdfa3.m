function [] = plotdfa3(U,DF_x,DF_y,DF_z,names)
%plotdfa3(U,DF_x,DF_y,DF_z,names)
%Plot DFA/CVA with labels, title, etc
%where U = DF scores matrix
%DF_x, DF_y  and DF_z are the three ordinates to be plotted
%names contains integer identifier file
%
% Copyright (c) 1996, Royston Goodacre
%

%asks for title
tit = input('Title please (in quotes)  ');

figure

%turns number to integers for labels
DF_X1=num2int(DF_x);
DF_Y1=num2int(DF_y);
DF_Z1=num2int(DF_z);
DF_X2=['DF ',DF_X1];
DF_Y2=['DF ',DF_Y1];
DF_Z2=['DF ',DF_Z1];

%plot with axis labels
plotnm3(U,DF_x,DF_y,DF_z,names);
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',[14])
hold on
xlabel(DF_X2)
ylabel(DF_Y2)
zlabel(DF_Z2)
title(tit,'FontSize',[16],'FontWeight','bold')
hold off

