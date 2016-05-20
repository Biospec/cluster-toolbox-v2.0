function [] = plotftir(M)
%plotftir(M)
%Plot axis correctly for FT-IR in colour
%where M = IR spectrum or spectra
%
%colours=['rbgcm']';
%
% Copyright (c) 1997, Royston Goodacre
%

%asks for title
tit = input('Title please (in quotes)  ');

%new figure
figure

%finds max and min of matrix M
maxvalM=max(max(M'));
minvalM=min(min(M'));
[rowM,colM]=size(M);

%set up colour index ('y' is missing but could be added)
%set up size so if number of spectra >5 loop idx matrix
%if you want reps used replace colidx with the correct order
colours=['rbgcm']';
numcols=floor((rowM/5))+1;
for i=1:numcols
  colidx(((i*5)-4):(i*5),:)=colours;
end
colidx=colidx';

%creates matrix of correct dimensions for FTIR
ut = lintrans(1:colM,1,colM,4000,600);

if rowM==1

plot(ut,M(1,:),'r')

else

plot(ut,M(1,:),'r')
hold on

%plots multiple pictures
for i=2:rowM
 plot(ut,M(i,:),colidx(1,i))
end

end

%sets X and Y limits so have +/- 5% on absorbances
axis([500 4100 (minvalM-((maxvalM-minvalM)/20)) (maxvalM+((maxvalM-minvalM)/20))])

%reverse X axis
set(gca,'xdir','reverse');
set(gca,'FontName','Times New Roman')
set(gca,'FontSize',[14])
hold on
xlabel('Wavenumber (cm^-^1)')
ylabel('Absorbance (arbitrary)')
title(tit,'FontSize',[16],'FontWeight','bold')




