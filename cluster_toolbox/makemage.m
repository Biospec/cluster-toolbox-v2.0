function [] = makemage(filename,U,name);
%makemage(filename,U,name)
%This function writes a file in kinemage format
%The data are automatically scaled
%
%   For the filename - remember I need 'file.txt'
%   U      -  is the co-ordinates for X, Y and Z
%   name   -  is an string matrix containing names
%
% Copyright (c) 1997, Royston Goodacre
%

[Nocols,Norows]=size(U);
[OverallMax]=max(max(abs(U)));
[Uscaled]=(U/OverallMax)*10;
[maxval]=max(max(Uscaled));
[minval]=min(min(Uscaled));


fid = fopen(filename,'w');
   fprintf(fid, '@kinemage 1 '); fprintf(fid, '\n');
   fprintf(fid, '@group {CVA plots}'); fprintf(fid, '\n');
   fprintf(fid, '@viewid {X vs Y}'); fprintf(fid, '\n');
   fprintf(fid, '@zoom 0.7'); fprintf(fid, '\n');
   fprintf(fid, '@zslab 800'); fprintf(fid, '\n');
   fprintf(fid, '@ztran 200'); fprintf(fid, '\n');
   fprintf(fid, '@center 0.0 0.0 0.0'); fprintf(fid, '\n');
   fprintf(fid, '@matrix 1 0 0 0 1 0 0 0 1'); fprintf(fid, '\n');
   fprintf(fid, '@2viewid {X vs Z}'); fprintf(fid, '\n');
   fprintf(fid, '@2zoom 0.7'); fprintf(fid, '\n');
   fprintf(fid, '@2zslab 800'); fprintf(fid, '\n');
   fprintf(fid, '@2ztran 200'); fprintf(fid, '\n');
   fprintf(fid, '@2center 0.0 0.0 0.0'); fprintf(fid, '\n');
   fprintf(fid, '@2matrix 1 0 0 0 0 -1 0 1 0'); fprintf(fid, '\n');
   fprintf(fid, '@3viewid {Y vs Z}'); fprintf(fid, '\n');
   fprintf(fid, '@3zoom 0.7'); fprintf(fid, '\n');
   fprintf(fid, '@3zslab 800'); fprintf(fid, '\n');
   fprintf(fid, '@3ztran 200'); fprintf(fid, '\n');
   fprintf(fid, '@3center 0.0 0.0 0.0'); fprintf(fid, '\n');
   fprintf(fid, '@3matrix 0 0 1 1 0 0 0 1 0'); fprintf(fid, '\n');
   fprintf(fid, '@labellist {samples} color=blue'); fprintf(fid, '\n');

      for i=1:Nocols;
         fprintf(fid, '{');
         fprintf(fid, name(i,:));
         fprintf(fid, '}	');
         fprintf(fid, '%12.8f	', Uscaled(i,1:(Norows-1)));
         fprintf(fid, '%12.8f', Uscaled(i,Norows));
         fprintf(fid, '\n');
         i=i+1;
      end;

   fprintf(fid, '@vectorlist {axis} color=black'); fprintf(fid, '\n');
      fprintf(fid, 'P 0 0 0 {X} '); fprintf(fid, '%f ', maxval); fprintf(fid, '0 0'); fprintf(fid, '\n');
      fprintf(fid, 'P 0 0 0 {Y} 0 '); fprintf(fid, '%f ', maxval); fprintf(fid, '0'); fprintf(fid, '\n');
      fprintf(fid, 'P 0 0 0 {Z} 0 0 '); fprintf(fid, '%f', maxval); fprintf(fid, '\n');

   fprintf(fid, '@labellist {CVA labels} color=black'); fprintf(fid, '\n');
      fprintf(fid, '{X} '); fprintf(fid, '%f ', maxval); fprintf(fid, '0 0'); fprintf(fid, '\n');
      fprintf(fid, '{Y} 0 '); fprintf(fid, '%f ', maxval); fprintf(fid, '0'); fprintf(fid, '\n');
      fprintf(fid, '{Z} 0 0 '); fprintf(fid, '%f', maxval); fprintf(fid, '\n');

   fprintf(fid, '@whitebkg'); fprintf(fid, '\n');
   fprintf(fid, '@keepthinline'); fprintf(fid, '\n');

fclose(fid);

runmage=['!d:\programs\mage\mage_4_5.exe ',filename];
eval(runmage);

