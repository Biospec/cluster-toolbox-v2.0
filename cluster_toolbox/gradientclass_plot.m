function gradientclass_plot(u, label1,symbol)

%gradientclass_plot(u, label1, symbol)
%u: the data matrix 
%label1: the class label which shall be separated by colour gradient
%symbol: the symbol to be used in the plot, can be one of the following:
%   'o';'x';'+';'*';'s';'d';'v';'^';'<';'>';'p';'h';'.'
%   'o' will be used if ignored.


if iscell(label1)
    unique_color=unique(label1);
    no_colors=length(unique_color);
elseif ischar(label1)
    unique_color=cellstr(unique(label1,'rows'));
    no_colors=length(unique_color);
    label1=cellstr(label1);
else
    unique_color=cellstr(num2str(unique(label1)));
    no_colors=length(unique_color);
    label1=cellstr(num2str(label1));
end

color_map=jet(no_colors);

if nargin<3
    symbol='o';
end


figure
hold on

for i=1:no_colors
    sample_idx=find(ismember(label1, unique_color{i}));
    h=plot(u(sample_idx,1), u(sample_idx,2), symbol);
    set(h,'color', color_map(i,:))
    label_legend{i}=[unique_color{i}];
end

legend(label_legend);
% h=xlabel('Discriminant function 1'); set(h,'fontsize',14)
% h=ylabel('Discriminant function 2'); set(h,'fontisze',14)

                
