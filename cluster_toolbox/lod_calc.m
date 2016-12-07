function lod = lod_calc(x,y,alpha)
% Limit of detection calculation based on known concentration x and
%   response y
% x = known/reference concentrations
% y = response
% alpha = specified false positve and negatvie error, set to 0.05 if
%    ignored.

if isempty(alpha)
    alpha = 0.05;
end
n = length(x);
unique_x = unique(x); m = length(unique_x);
x_bar = mean(x);
y_bar = mean(y);
Qxx = 0; Qxy = 0; r = zeros(m,1);
for i=1:m
    idx = find(x == unique_x(i));
    r(i) = length(idx);
    Qxx = Qxx + r(i)*(unique_x(i) - x_bar)^2;
    y_i = mean(y(x == unique_x(i)));
    Qxy = Qxy + r(i)*(unique_x(i) - x_bar)*(y_i - y_bar);
end
w0_sq = 1/median(r) +1/n +x_bar^2/Qxx;
b = pinv([ones(n,1) x])*y;
y_hat = [ones(n,1) x]*b;
% alpha_hat = b(1);
beta_hat = b(2);
sigma_sq = sum((y_hat-y).^2)/(n-2);
delta_pq = tinv(1-alpha,n-2)*2;
lod = delta_pq*sqrt(w0_sq)*sqrt(sigma_sq)/beta_hat;