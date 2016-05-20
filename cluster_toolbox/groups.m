function [group] = groups(samps,reps,A)
%[group] = groups(samps,reps,A)
%returns a group files for dfa in order
%  A=1      ABC repeat (i.e., A,B,C, ...., A,B,C etc)
%  A=2      AAA repeat (i.e., A,A,A, ...., Z,Z,Z)        
%samps = number of groups
%reps = number of replicates
%
% Copyright (c) 1996, Royston Goodacre
%

if A == 1

%returns a group files for dfa in order ABC repeat (i.e., A,B,C, ...., A,B,C etc)
hh = (1:samps)';
group=[];

for i=1:reps
 group = [group;hh];
end

elseif A==2

%returns a group files for dfa in order AAA repeat (i.e., A,A,A, ...., Z,Z,Z)
group=[];

for i=1:samps
 for j=1:reps
  group = [group;i];
 end
end

end
