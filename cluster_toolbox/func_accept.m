function func_accept

global fig_h title_text axe_superimpose axe_bc edt_lamda edt_p d;
prog=get(fig_h,'UserData');
point=prog.point;
chroms=prog.chroms;
no_chroms=size(chroms,1);
prog.chroms_bc(point,:)=chroms(point,:)-prog.temp_tic;
point=point+1;
if point>no_chroms
    msgbox('Baseline Correction Finished');
    assignin('base','chroms_bc',prog.chroms_bc);
    close(fig_h);
    return
end


prog.point=point;
set(fig_h,'UserData',prog);
disp_string=['No. ' num2str(point) ', ' num2str(no_chroms) ' chromatogram(s) in total']; 
title_text=uicontrol('Parent', fig_h,...
    'Units', 'normalized',...
    'Position',[.12 .96 .6 .028],...
    'String',disp_string,...
    'FontSize', 14,...
    'FontWeight', 'bold',...
    'Style','Text',...
    'Tag','tip_text');
feval('func_show');