function [differ] = derivats(M,window)
%[differ] = derivats(M,window)
%takes Savitzky & Golay derivatives
%to take *1st* differential do [df1]=derivats(M,window);
%to take *1st* differential do [df2]=derivats(df1,window);
%M = FT-IR matrix 
%window = sliding window size (NORMALLY SET TO 5, accepts 5, 7 and 9)
%nb you will lose the ends of spectrum due to half window size
%
% Copyright (c) 1996, Royston Goodacre
% 

%sizing and calculating ends for M
[No_samps,No_waves] = size(M);

%returns an integer rounded down of half window
ends = floor((window/2));

%to take *1st* differential of Savitzky & Golay
df1=deriver(M',window);
differ=df1((1+ends):(No_waves-ends),:)';

%calculating ends for diff1
%[No_samps,No_diff1s] = size(diff1);

%to take *2nd* differential of Savitzky & Golay
%df2=deriver(diff1',window);
%diff2=df2((1+ends):(No_diff1s-ends),:)';




