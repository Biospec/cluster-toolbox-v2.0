function plot_idx=multiclass_plot(u, label, plot_idx)
%plot_idx=multiclass_plot(u, label)
%u: the data matrix
%label: the class label
%plot_idx: a numeric index which records the combination 
%   of the color, symbol and fill/not fill used in the figure. 
%   Will be randomly generated if not provided as an input, 
%   the output plot_idx can be used as the third input to 
%   produce more figures using the same combination.
%note: the supported maximum number of different classes  
%   is 155.

color=['r';'g';'b';'c';'m';'k';'y'];
symbol=['o';'x';'+';'*';'s';'d';'v';'^';'<';'>';'p';'h';'.'];
to_fill=[0;1];
if ischar(label)
    unique_label=unique(label,'rows');
else
    unique_label=unique(label);
end
no_cls=length(unique_label);
count=1;
for i=1:length(color)
    for ii=1:length(symbol)
        for iii=1:length(to_fill)
            symbol_to_plot.color(count)=color(i);
            symbol_to_plot.symbol(count)=symbol(ii);
            symbol_to_plot.fill(count)=to_fill(iii);
            if any(ismember({'x','+','.','*'}, symbol(ii))) && to_fill(iii)==1
                continue
            else
                count=count+1;
            end
        end
    end
end

no_combos=count-1;
if nargin<3
    plot_idx=randperm(no_combos);
    plot_idx=plot_idx(1:no_cls);
end

figure
hold on
for i=1:no_cls
    if iscell(label)
        h=plot(u(ismember(label,unique_label(i)),1),u(ismember(label,unique_label(i)),2),symbol_to_plot.symbol(plot_idx(i)));
        set(h,'color',symbol_to_plot.color(plot_idx(i)));
        if symbol_to_plot.fill(plot_idx(i))==1
            set(h,'markerfacecolor',symbol_to_plot.color(plot_idx(i)));
        end
    elseif ischar(label)
        h=plot(u(ismember(label,unique_label(i,:),'rows'),1),u(ismember(label,unique_label(i,:),'rows'),2),symbol_to_plot.symbol(plot_idx(i)));
        set(h,'color',symbol_to_plot.color(plot_idx(i)));
        if symbol_to_plot.fill(plot_idx(i))==1
            set(h,'markerfacecolor',symbol_to_plot.color(plot_idx(i)));
        end
    else
        h=plot(u(label==unique_label(i),1),u(label==unique_label(i),2),symbol_to_plot.symbol(plot_idx(i)));
        set(h,'color',symbol_to_plot.color(plot_idx(i)));
        if symbol_to_plot.fill(plot_idx(i))==1
            set(h,'markerfacecolor',symbol_to_plot.color(plot_idx(i)));
        end
    end

end
if isnumeric(unique_label)
    unique_label=num2str(unique_label);
end
legend(unique_label);
% h=xlabel('Discriminant function 1'); set(h,'fontsize',14)
% h=ylabel('Discriminant function 2'); set(h,'fontisze',14)
plot_idx=plot_idx;
                
