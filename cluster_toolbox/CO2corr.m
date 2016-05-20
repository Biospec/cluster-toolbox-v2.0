function [co2basl] = CO2corr(M,type)
%[co2basl] = CO2corr(M,type)
%takes matrix M and corrects for CO2 interference
%CO2 bands which are treated are:
%2403.21 to 2272.06 (i.e. bins closest to 2400-2275)
%682.77 to 655.77 (i.e. bins closest to 680-660)
%when:
%      type = 1, bands are filled with a trend
%      type = 2, bands are replaced with zeros
%      type = 3, bands are chopped out
%
% Copyright (c) 1997, Royston Goodacre
%

%to work out matrix dimensions
[rows,cols]=size(M);


if type == 1

%to set a trend in the two CO2 bands
co2basl1=zeros(rows,35);
for i=1:rows
   co2begin1=M(i,415);
   co2end1=M(i,449);
   slide1 = lintrans(1:35,1,35,co2begin1,co2end1);
   co2basl1(i,:)=slide1;
end

co2basl2=zeros(rows,8);
for i=1:rows
   co2begin2=M(i,861);
   co2end2=M(i,868);
   slide2 = lintrans(1:8,1,8,co2begin2,co2end2);
   co2basl2(i,:)=slide2;
end

for i=1:rows
   co2basl(i,:)=[M(i,1:414)';co2basl1(i,:)';M(i,450:860)';co2basl2(i,:)';M(i,869:882)']';
end


elseif type == 2

%to set a two CO2 bands to zeros
co2basl1=zeros(rows,35);
co2basl2=zeros(rows,8);
for i=1:rows
   co2basl(i,:)=[M(i,1:414)';co2basl1(i,:)';M(i,450:860)';co2basl2(i,:)';M(i,869:882)']';
end


elseif type == 3

%to remove the two CO2 bands
for i=1:rows
   co2basl(i,:)=[M(i,1:414)';M(i,450:860)';M(i,869:882)']';
end


end
