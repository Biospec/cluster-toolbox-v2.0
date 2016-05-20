function gradientclass_plot2(u, label1, label2)

%gradientclass_plot(u, label1, label2)
%u: the data matrix 
%label1: the class label which shall be separated by colour gradient
%label2: the class label which shall be separated by different symbols


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
end

if isnumeric(label2)
    label2=cellstr(num2str(label2));
elseif ischar(label2)
    label2=cellstr(label2);
end

color_map=jet(no_colors);

symbol=['o';'x';'+';'*';'s';'d';'v';'^';'<';'>';'p';'h';'.'];

figure
hold on
count=1;
for i=1:no_colors
    label2_now=label2(ismember(label1, unique_color{i}));
    symbols_now=unique(label2_now);
    no_symbols=length(symbols_now);
    for ii=1:no_symbols
        if ii>length(symbol)
            symbol_idx=rem(ii, length(symbol));
        else
            symbol_idx=ii;
        end
        sample_idx=find(ismember(label1, unique_color{i}) & ismember(label2, symbols_now{ii}));
        h=plot(u(sample_idx,1), u(sample_idx,2), symbol(symbol_idx));
        set(h,'color', color_map(i,:))
        label_legend{count}=[unique_color{i} ' ' symbols_now{ii}];
        count=count+1;
    end
end

legend(label_legend);
% h=xlabel('Discriminant function 1'); set(h,'fontsize',14)
% h=ylabel('Discriminant function 2'); set(h,'fontisze',14)

                
