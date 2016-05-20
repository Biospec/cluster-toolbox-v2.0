function [peak_areas, data_int, baseline]=band_area2(data, start_pt, end_pt)

data_int=data(:,start_pt:end_pt);
for i=1:size(data,1)
    data_int2(i,:)=data_int(i,:);
    baseline(i,:)=linspace(data_int2(i,1), data_int2(i,end), length(data_int2(i,:)));
    data_int2(i,:)=data_int2(i,:)-baseline(i,:);
    peak_areas(i)=trapz(data_int2(i,:));
end
peak_areas(peak_areas<0)=NaN;