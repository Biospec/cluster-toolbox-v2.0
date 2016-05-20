function [peak_areas, data_int, baseline]=band_area(data)

figure
plot(data');
disp('Click the start and end point of the band and then press enter to continue');
pause;
[x,y]=ginput;
x=round(x);
data_int=data(:,x(1):x(2));
for i=1:size(data,1)
    data_int2(i,:)=data_int(i,:);
    baseline(i,:)=linspace(data_int2(i,1), data_int2(i,end), length(data_int2(i,:)));
    data_int2(i,:)=data_int2(i,:)-baseline(i,:);
    peak_areas(i)=trapz(data_int2(i,:));
end
peak_areas(peak_areas<0)=NaN;