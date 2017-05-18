function signal_out = sg_process(signal_in,w,p,d)
% signal_out = sg_process(signal_in,w,p,d)
% Savitzky-Golay signal processing
% signal_in: spectra matrix
% w: window width
% p: order of polynomial
% d: d=0 for smoothing, d>1 for approximation of dth derivative

[~,g] = sgolay(p,w);
no_spc = size(signal_in,1);
dt = 1;
signal_out = zeros(size(signal_in));
for i=1:no_spc
    signal_out(i,:) = conv(signal_in(i,:)', factorial(d)/(-dt)^d*g(:,d+1), ...
        'same')';
end

