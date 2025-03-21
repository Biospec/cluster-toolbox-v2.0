function output_txt = var_id_from_WS_datatip(obj,event_obj)
% Display data cursor position in a data tip
% obj          Currently not used
% event_obj    Handle to event object
% output_txt   Data tip text, returned as a character vector or a cell array of character vectors

pos = event_obj.Position;
Data_idx = event_obj.DataIndex;
var_id_in = evalin('base', 'var_id4;');
%********* Define the content of the data tip here *********%

output_txt = [num2str(var_id_in(Data_idx))];
%***********************************************************%



%***********************************************************%

function formattedValue = formatValue(value,event_obj)
% If you do not want TeX formatting in the data tip, uncomment the line below.
% event_obj.Interpreter = 'none';
if strcmpi(event_obj.Interpreter,'tex')
    valueFormat = ' \color[rgb]{0 0.5 1}\bf';
    removeValueFormat = '\color[rgb]{.15 .15 .15}\rm';
else
    valueFormat = ': ';
    removeValueFormat = '';
end
formattedValue = [valueFormat num2str(value,4) removeValueFormat];
