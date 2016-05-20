function output = gaussian_smooth(input, window_width)
% output = gaussian_smooth(input, window_width)
% Gaussian smooth function 
% input: a vector of signal need to be smoothed
% window_width: the length of the moving window
% output: the smoothed signal

% window_width=int16(window_width);
half_width=round(window_width/2);
gaussFilter = gausswin(double(window_width));
gaussFilter=gaussFilter/sum(gaussFilter);
output=conv(input, gaussFilter);
output=output(half_width:end-half_width);
