function gradientclass_plot2(u, label1, label2)

%gradientclass_plot(u, label1, label2)
%u: the data matrix with no more than 3 columns. 
%label1: the class label which shall be separated by colour gradient
%label2: the class label which shall be separated by different symbols


if iscell(label1)
    unique_label1 = unique(label1); 
    num_flag = true;
    for i = 1:numel(unique_label1)
        if isempty(str2num(unique_label1{i}))
            num_flag = false;
        end
    end
    if ~num_flag
        unique_color=unique(label1);        
    else        
        for i = 1:size(label1, 1)
            label1_tmp(i,1) = str2num(label1{i});
        end
        unique_label1 = unique(label1_tmp);
        for i = 1:numel(unique_label1)
            unique_color{i} = num2str(unique_label1(i));
        end
        clear label1_tmp
    end
    no_colors=length(unique_color);
elseif ischar(label1)
    num_flag = true;
    unique_label1 = unique(label1, 'rows');
    for i = 1:size(unique_label1, 1)
        if isempty(str2num(unique_label1(i, :)))
            num_flag = false;
        end
    end
    if ~num_flag
        unique_color = cellstr(unique(label1,'rows'));
        no_colors = length(unique_color);
        label1 = cellstr(label1);
    else
        for i = 1:size(label1, 1)
            label1_tmp(i, 1) = str2num(label1(i, :));
        end
        unique_label1 = unique(label1_tmp);
        for i = 1:numel(unique_label1)
            unique_color{i} = num2str(unique_label1(i));
        end
        no_colors = length(unique_color);
        label1 = cellstr(label1);
    end
else
    unique_color=cellstr(num2str(unique(label1)));
    no_colors=length(unique_color);
    label1 = cellstr(num2str(label1));
end

if isnumeric(label2)
    label2=cellstr(num2str(label2));
elseif ischar(label2)
    label2=cellstr(label2);
end

color_map=parula(no_colors);

symbol=['o';'s';'d';'v';'^';'<';'>';'p';'h';'x';'+';'*';'.'];

n = size(u, 2);
if n == 2
    figure % 2-D plot
    hold on
    set(gcf, 'colormap', parula(256));
    count=1;
    for i=1:no_colors
        label1_now=label1(ismember(label1, unique_color{i}));
        symbols_now=unique(label2);
        no_symbols=length(symbols_now);
        for ii=1:no_symbols
            if ii>length(symbol)
                symbol_idx=rem(ii, length(symbol));
            else
                symbol_idx=ii;
            end
            sample_idx=find(ismember(label1, unique_color{i}) & ismember(label2, symbols_now{ii}));
            h=plot(u(sample_idx,1), u(sample_idx,2), symbol(symbol_idx));
            set(h,'markersize', 8, 'color', color_map(i,:), 'markerfacecolor', ...
                color_map(i, :))
            if i == 1
                label_legend{count}=symbols_now{ii};
                count=count+1;
            end
        end
    end
    h = colorbar;
    set(h, 'Ytick', linspace(0, 1, length(unique_color)), 'YtickLabel', unique_color);
    legend(label_legend);
    % h=xlabel('Discriminant function 1'); set(h,'fontsize',14)
    % h=ylabel('Discriminant function 2'); set(h,'fontisze',14)
else
    if n > 3
        warning('There are more than 3 columns in the data matrix, only first 3 are used to generate the 3-D plot!')
    end
    figure % 3-D plot
    hold on
    set(gcf, 'colormap', parula(256));
    count=1;
    for i=1:no_colors
        label1_now=label1(ismember(label1, unique_color{i}));
        symbols_now=unique(label2);
        no_symbols=length(symbols_now);
        for ii=1:no_symbols
            if ii>length(symbol)
                symbol_idx=rem(ii, length(symbol));
            else
                symbol_idx=ii;
            end
            sample_idx=find(ismember(label1, unique_color{i}) & ismember(label2, symbols_now{ii}));
            h=plot3(u(sample_idx,1), u(sample_idx,2), u(sample_idx, 3), symbol(symbol_idx));
            set(h,'markersize', 8, 'color', color_map(i,:), 'markerfacecolor', ...
                color_map(i, :))
            if i == 1
                label_legend{count}=symbols_now{ii};
                count=count+1;
            end
        end
    end
    h = colorbar;
    set(h, 'Ytick', linspace(0, 1, length(unique_color)), 'YtickLabel', unique_color);
    legend(label_legend);
    % h=xlabel('Discriminant function 1'); set(h,'fontsize',14)
    % h=ylabel('Discriminant function 2'); set(h,'fontisze',14)
end
                
