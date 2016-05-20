function [] = floatout(M,filename,type)
%[]=floatout(M,filename,type)
%This function outputs floating point data of 12.8
%of the matrix M and stores it in filename - remember I need 'file.txt'
%if type = 1 then tab separated
%          2 then comma separated
%
% Copyright (c) 1996, Royston Goodacre
%


if type==1

fid = fopen(filename,'w');

[Nocols,Norows]=size(M);
for i=1:Nocols;
   fprintf(fid, '%12.8f	', M(i,1:(Norows-1)));
   fprintf(fid, '%12.8f', M(i,Norows));
   fprintf(fid, '\n');
   i=i+1;
end;
fclose(fid);
clear Nocols Norows ans fid i



elseif type ==2

fid = fopen(filename,'w');

[Nocols,Norows]=size(M);
for i=1:Nocols;
   fprintf(fid, '%12.8f,', M(i,1:(Norows-1)));
   fprintf(fid, '%12.8f', M(i,Norows));
   fprintf(fid, '\n');
   i=i+1;
end;
fclose(fid);
clear Nocols Norows ans fid i

end
