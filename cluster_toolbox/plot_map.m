function [] = plot_map(Surface,type)
%plot_map(Surface,type)
%Plots a map file of size according to matrix
%where Surface is coordinate
%type can be:   1 = colour 2-D contour map
%               2 = b&w 3-D surface map
%               3 = combination of 1 and 2
%
% Copyright (c) 1996, Royston Goodacre
%

%asks for title
tit = input('Title please (in quotes)  ');

figure

if type == 1
   pcolor(Surface)
   shading interp
   colormap(jet)
   colorbar('vert')

elseif type == 2
   surfl(Surface)
   view(276,60)
   shading interp
   colormap(gray)

elseif type == 3
   surfc(Surface)
   shading interp

end

hold on
title(tit,'FontSize',[16],'FontWeight','bold')
hold off

