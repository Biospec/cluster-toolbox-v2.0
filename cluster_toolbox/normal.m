function ans = normal(spec);
% ans = normal(spec)
% spec is single spectrum or multiple spectra of form spec(1:number of spectra, 1:number of bins)
% Normalise spectrum to top in user-specified zone. Bottom is lowest point in spec,

topscale = 1000; % Scaled in range 0 to this number

[numspecs, numbins] = size(spec);

%prompt = ['Start of range for spectrum top (1)'];
%range1 = input(prompt);

%if isempty(range1)
   range1 = 1;
%end%if

%prompt = ['End of range for spectrum top (', num2str(numbins), ')'];
%range2 = input(prompt);

%if isempty(range2)
   range2 = numbins;
%end%if

for i=1:numspecs
   top = max(spec(i,range1:range2));
   bottom = min(spec(i,1:numbins));
   top = top - bottom;
   spec(i,:) = (spec(i,:)-bottom)*topscale/top;
end;

ans = spec;
