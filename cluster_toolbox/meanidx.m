function [MeanM] = meanidx(M,group)
%[MeanM] = meanidx(M,group)
%this calculates means of M based on index from a group file
%the sizes of the group file can be different
%use namesidx with this to get names correct
%
% Copyright (c) 2001, Royston Goodacre
%

%find the unique objects in group
classes=unique(group);

%find the number of different groups
len=length(classes);

%calculates means based on sub-index using findstr.m
for i=1:len
 MeanM(i,:)=mean(M([findstr(classes(i,:),group')],:));
end


