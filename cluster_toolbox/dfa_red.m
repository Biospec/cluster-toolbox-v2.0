function [] = dfa_red(tt,DF_x,DF_y,names)
%plot_pca(tt,PC_x,PC_y,names)
%Plot PCA with labels, title, etc
%where tt = PC scores matrix
%PC_x and PC_y are the two ordinates to be plotted
%names contains integer identifier file
%
% Copyright (c) 1996, Royston Goodacre
%

%asks for title
%tit = input('Title please (in quotes)  ');
%tit='DFA; red = train, blue = test';

figure

%turns number to integers for labels
DF_X1=num2int(DF_x);
DF_Y1=num2int(DF_y);
DF_X2=['Discriminant function ',DF_X1];
DF_Y2=['Discriminant function ',DF_Y1];

%plot with axis labels
pltnmred(tt,DF_x,DF_y,names);
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',[14])
hold on
xlabel(DF_X2)
ylabel(DF_Y2)
%title(tit,'FontSize',[16],'FontWeight','bold')
hold off

