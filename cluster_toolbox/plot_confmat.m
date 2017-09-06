function handles = plot_confmat(conf_mat)
% handles = plot_confmat(conf_mat)
%   plot a confusion matrix with numbers on top of it
%   conf_mat: the confusion matrix in %.
[m,n] = size(conf_mat);
if max(max(conf_mat))<1
    conf_mat = conf_mat*100;
end
conf_mat = round(conf_mat*100)/100;
figure
h=image(conf_mat);
set(h,'CDatamapping','scaled');
for i=1:m
    for ii=1:n
        h(i,ii) = text(ii-.1, i, [num2str(conf_mat(i,ii)) '%']);
        if i~=ii
            set(h(i,ii),'color','w');
        end
    end
end
set(gca,'XTick', 1:m);
set(gca,'YTick', 1:n);
handles = h;
h=title('Confusion matrix');
set(h,'fontsize',14)

end

