function [name] = names(samps,reps,A)
%[name] = names(samps,reps,A)
%returns a names (numerical integer) files for plotting in order
%  A=1      numerical ABC repeat (i.e., A,B,C, ...., A,B,C etc)
%  A=2      numerical AAA repeat (i.e., A,A,A, ...., Z,Z,Z)        
%samps = number of groups
%reps = number of replicates
%
% Copyright (c) 1996, Royston Goodacre
%

if A == 1

%returns a numerical matrix for names in order ABC repeat (i.e., A,B,C, ...., A,B,C etc)
hh = (1:samps)';
nam=[];

for i=1:reps
 nam = [nam;hh];
end

elseif A==2

%returns a numerical matrix for names in order AAA repeat (i.e., A,A,A, ...., Z,Z,Z)
nam=[];

for i=1:samps
 for j=1:reps
  nam = [nam;i];
 end
end

end


%Converts Matrix file nam into integer file

[Mrows,Mcols] = size(nam);
name=[];

for i = 1:Mrows
 hh=sprintf('%-5.0f',nam(i,1));
 name=[name;hh];
end
