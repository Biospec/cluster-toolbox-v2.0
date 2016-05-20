function [nameM] = filtnam(M,samps,reps,A)
%[nameM] = filtnam(M,samps,reps,A)
%this calculates means based on :
%  A=1      ABC repeat (i.e., A,B,C, ...., A,B,C etc)
%  A=2      AAA repeat (i.e., A,A,A, ...., Z,Z,Z)        
%M = matrix to mean
%samps = number of groups (excluding reps)
%reps = number of replicates
%
% Copyright (c) 1996, Royston Goodacre
%

[Mrows,Mcols] = size(M);
idx=1:Mrows;

%to choose AAA or ABC order
if A==1
 idxreshape=reshape(idx,samps,reps);
 elseif A==2
 idxreshape=reshape(idx,reps,samps)';
end

%this calculates means based on above choice
for i=1:samps
 nameM(i,:)=M(idxreshape(i,:));
end

