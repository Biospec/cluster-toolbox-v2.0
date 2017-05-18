function mat_out = normminmax(mat_in,min_val,max_val)
% mat_out = normminmax(mat_in,min_val,max_val)
% Normalise each row to a specific range, defined by min_val and max_val.
% If min_val, max_val ignored, it will normalise each row to a range of -1
% to 1.

[m,n] = size(mat_in);
if nargin == 1
    min_val = -1;
    max_val = 1;
end
mat_out = zeros(size(mat_in));
for i=1:m
    a = (max_val - min_val)/(max(mat_in(i,:))-min(mat_in(i,:)));
    b = max_val - a*max(mat_in(i,:));
    mat_out(i,:) = a.*mat_in(i,:) + b;
end

