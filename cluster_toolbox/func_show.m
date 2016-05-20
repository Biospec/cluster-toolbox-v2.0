function func_show

global fig_h axe_superimpose axe_bc button_show button_accept edt_lamda edt_p d;
set(button_show,'Enable','off');
set(button_accept,'Enable','off');
prog=get(fig_h,'UserData');
chroms=prog.chroms; point=prog.point;
lamda=str2num(get(edt_lamda,'String'));
p=str2num(get(edt_p,'String'));
temp_tic=asysm(chroms(point,:)',lamda,p,d);
temp_tic=temp_tic';
plot(axe_superimpose,chroms(point,:));
hold(axe_superimpose,'on');
plot(axe_superimpose,temp_tic,'r');
hold(axe_superimpose,'off');
title_superimpose=title(axe_superimpose,'Chromatogram before baseline correction', ...
    'FontSize', 14);
plot(axe_bc,chroms(point,:)-temp_tic);
title_bc=title(axe_bc,'Chromatogram after baseline correction', 'FontSize', 14);
prog.temp_tic=[];prog.temp_tic=temp_tic;
set(fig_h,'UserData',prog);
set(button_show,'Enable','on');
set(button_accept,'Enable','on');