function [] = plot_pca(tt,PC_x,PC_y,names, pr)
%plot_pca(tt,PC_x,PC_y,names, pr)
%Plot PCA with labels, title, etc
%where tt = PC scores matrix
%PC_x and PC_y are the two ordinates to be plotted
%names contains integer identifier file
%pr, total explained variance calculated by pca_np, if pr is provided, it
%will be shown in the captions of X and Y axis.
% Copyright (c) 1996, Royston Goodacre
% Updated May 2019, Yun Xu

%asks for title
tit = input('Title please (in quotes)  ');
if isnumeric(names)
    names = num2str(names(:));
end
if iscell(names)
    for i = 1:length(names)
        if isnumeric(names{i})
            names{i} = num2str(names{i});
        end
    end
end

figure

%turns number to integers for labels
PC_X1=num2int(PC_x);
PC_Y1=num2int(PC_y);
if ~exist('pr', 'var')
    PC_X2=['Principal component ',PC_X1];
    PC_Y2=['Principal component ',PC_Y1];
else
    if PC_x == 1
        PC_X2 = ['Principal component ', PC_X1, ', TEV = ', ...
            num2str(round(pr(PC_x)*100)/100), '%'];
    else
        PC_X2 = ['Principal component ', PC_X1, ', TEV = ', ...
            num2str(round((pr(PC_x) - pr(PC_x - 1))*100)/100), '%'];
    end
    if PC_y == 1
        PC_Y2 = ['Principal component ', PC_Y1, ', TEV = ', ...
            numstr(round(pr(PC_y)*100)/100), '%'];
    else
        PC_Y2 = ['Principal component ', PC_Y1, ', TEV = ', ...
            num2str(round((pr(PC_y) - pr(PC_y - 1))*100)/100), '%'];
    end
end

%plot with axis labels
plotnm2(tt,PC_x,PC_y,names);
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',[14])
hold on
xlabel(PC_X2)
ylabel(PC_Y2)
title(tit,'FontSize',[16],'FontWeight','bold')
hold off

