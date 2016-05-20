function [idx] = makeidx(groups,subset)
%[idx] = makeidx(groups,want)
%returns an index (idx) files base on find numbers you
%want in a groups file
%
%   idx = index file output in order of occurence
%   groups = replicate indexing
%   subset = subset that you want from groups
%          single row vector of type [a b c...]
%
% Copyright (c) 1996, Royston Goodacre
%

lengthsubset=length(subset);

%for first instance
whereis=find(groups==subset(1,1));
lengthwhereis=length(whereis);
idx((1):(lengthwhereis),1)=whereis((1:lengthwhereis),1);


for i=2:lengthsubset
   whereis=find(groups==subset(1,i));
   lengthwhereis=length(whereis);
   lengthidx=length(idx);
   idx((lengthidx+1):(lengthidx+lengthwhereis),1)=whereis((1:lengthwhereis),1);
end

