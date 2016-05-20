function [NewName] = namesidx(name,group)
%[NewName] = namesidx(name,group)
%this orders names in order of groups on first occurence
%so can be used with meanidx for resorting names
%
%use when you have used meanidx to take the mean of a matrix
%
% Copyright (c) 1997, Royston Goodacre
%

%find the unique objects in group
classes=unique(group);

%find the number of different groups
len=length(classes);

%this orders names according to groups
for i=1:len
 where=findstr(classes(i,:),group');
 NewName(i,:)=name(where(1,1),:);
end
