function gradientclass_plot(u, label1, symbol, xaxis)

%gradientclass_plot(u, label1, symbol, xaxis)
%u: the data matrix, if there are 3 or less columns in u, 2D or 3D plot 
%   will be generated, if there are more than 3 columns, the whole row 
%   will be ploted (spectra mode)
%label1: the class label which shall be separated by colour gradient, it
%   can be numbers, char array or cell array
%symbol: the symbol to be used in the plot, can be one of the following:
%   'o';'x';'+';'*';'s';'d';'v';'^';'<';'>';'p';'h';'.'
%   for spectra model, '-' will be automatically added ('solid line'), use 
%   [] as a space holder if you want to add xaxis to the figure (see below)
%   and keep using solid line 
%   if ignored 'o' will be used for 2D and 3D mode and '-' for spectra
%   mode.
%xaxis, the X axis unit, useful if plot u in spectra mode, will be ignored
%   if plot a 2D or 3D figure.
%By Yun Xu
%updated 23/05/2016: 3D plot and spectra mode added

[m,n]=size(u);
if n<2
    error('Need at least 2 columns in u!')
end
if nargin==4
    if length(xaxis)~=n
        error('The length of xaxis must be the same as the number of columns in u!')
    end
end
% parse the label
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

% Use Jet for gradient colouring
color_map=jet(no_colors);

if nargin<3
    symbol='o';
end


figure
hold on

for i=1:no_colors
    sample_idx=find(ismember(label1, unique_color{i}));
    switch n
        %2-D plot
        case 2
            h=plot(u(sample_idx,1), u(sample_idx,2), symbol);
            set(h,'color', color_map(i,:))
            label_legend{i}=[unique_color{i}];
        %3-D plot
        case 3
            h=plot3(u(sample_idx,1), u(sample_idx,2), u(sample_idx,3), symbol);
            set(h,'color', color_map(i,:))
            label_legend{i}=[unique_color{i}];
        otherwise
        %Spectra mode
            if nargin<4
                xaxis=1:n;
            end
            if strcmp(symbol,'o')
                symbol='-';
            elseif isempty(symbol) ||(length(symbol)==1 && ~strcmp(symbol, '-'))
                symbol=[symbol '-'];
            end
            h=plot(xaxis, u(sample_idx,:),symbol);
            set(h,'color',color_map(i,:));
            h_(i)=h(1);
            label_legend{i}=[unique_color{i}];
    end
end
axis('tight')
if n<=3
    legend(label_legend);
else
    legend(h_, label_legend);
end


                
