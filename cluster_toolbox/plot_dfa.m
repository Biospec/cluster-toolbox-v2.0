function [] = plot_dfa(U,DF_x,DF_y,names)
%plot_dfa(U,DF_x,DF_y,names)
%Plot DFA/CVA with labels, title, etc
%where U = DF scores matrix
%DF_x and DF_y are the two ordinates to be plotted
%If there is only 1 DF, set DF_y to 0
%names contains integer identifier file
%
% Copyright (c) 1997, Royston Goodacre
% Updated 2020, Yun Xu

%asks for title
tit = input('Title please (in quotes)  ');
if iscell(names)
    names = cell2mat(names);
end
if isnumeric(names)
    names = num2str(names(:));
end

if any([DF_x == 0, DF_y == 0])
    figure
    plot(1:numel(U), U, 'w.')
    text(1:numel(U), U, names);
    current_YLim = get(gca, 'YLim');
    set(gca, 'YLim', [current_YLim(1)*.9, ...
        current_YLim(2) + current_YLim(1)*0.1]);
    h = xlabel('Sample id.');
    set(h, 'fontsize', 14);
    h = ylabel('Discriminant function scores');
    set(h, 'fontsize', 14);
    title(tit, 'FontSize', 16, 'FontWeight', 'bold')
else
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
end
