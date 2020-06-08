function data_corr = co2corr(data_mat, wavenumber)
%A more adaptive CO2 correction routine which select CO2 regions (i.e.
% (1) between 655.77 and 682.77 cm^-1 and (2) between 2272.06 and 2403.21 cm^-1)
% automatically.
% data_corr = co2corr(data_mat, wavenumber)
% data_mat: a matrix of FT-IR spectra in which each row is a sample
% wavenumber: a vector providing wave numbers for the variables

var_seg1 = find(wavenumber >= 655.77 & wavenumber <= 682.77);
var_seg2 = find(wavenumber >= 2272.06 & wavenumber <= 2403.21);

if ~isempty(var_seg1)
    for i = 1:size(data_mat, 1)
        data_mat(i, var_seg1(1):var_seg1(end)) = linspace(data_mat(i, ...
            var_seg1(1)), data_mat(i, var_seg1(end)), length(var_seg1));
    end
end

if ~isempty(var_seg2)
    for i = 1:size(data_mat, 1)
        data_mat(i, var_seg2(1):var_seg2(end)) = linspace(data_mat(i, ...
            var_seg2(1)), data_mat(i, var_seg2(end)), length(var_seg2));
    end
end

data_corr = data_mat;
