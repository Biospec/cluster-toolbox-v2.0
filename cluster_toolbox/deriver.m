function[diff]=deriver(raw,punkt)
% function[diff]=deriver(raw,punkt)
% Savitzky & Golay-derivasjon langs n-retning, 5,7 eller 9 punkts
% antar:raw=n*m=(# retensjonstider)*

[n,m]=size(raw);
convolute5=[-1 8 0 -8 1];
convolute7=[-22 67 58 0 -58 -67 22];
convolute9=[-86 142 193 126 0 -126 -193 -142 86]; 

if (punkt == 5)
        convolute = convolute5; 
	number = 12;       
elseif (punkt == 7)
        convolute = convolute7;
	number = 252;
else
        convolute = convolute9;
	number = 1188;
end;
start=((punkt+1)/2);
slutt=(n-1+(punkt+1)/2);

for i=1:m
%        a(:)=raw(:,i);
        a=raw(:,i);
        b=conv(convolute,a);
        diff(:,i)=b(start:slutt)./number;
end;

