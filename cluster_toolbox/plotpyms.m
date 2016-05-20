function [] = plotpyms(M, ut)
%plotpyms(M, ut)
%Plots a single PyMS spectrum as % total
% M is the mass intensity and ut is the corresponding m/z
% Copyright (c) 1997, Royston Goodacre
%

%creates matrix of correct dimensions for PyMS
if nargin<2
    ut = 51:200;
end
%asks for title
tit = input('Title please (in quotes)  ');

%new figure
figure

%finds max and min of matrix M
maxvalM=max(max(M'));
minvalM=min(min(M'));
[rowM,colM]=size(M);

%turns into % ion count
totalions=sum(M);
percM=M*(100/totalions);

%does plot
bar(ut,percM(1,:),0.001,'k')

%sets Y limit so have + 5% on top
axis([51 200 0 ((maxvalM*100)+(((maxvalM*100)-(minvalM*100))/20))])

%labelling etc
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',[20])
set(gca,'TickDir','out')
hold on
xlabel('Mass (m/z)')
ylabel('Percentage total ion count')
title(tit,'FontSize',[26],'FontWeight','bold')

