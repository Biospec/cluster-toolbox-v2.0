function [tt,projXtt,pp,pr,U,projXU,V,eigenvals] = projpcdf(trainx,trainy,train_names,testx,test_names,n_pc,n_df)
%function [tt,projXtt,pp,pr,U,projXU,V,eigenvals] = projpcdf(trainx,trainy,train_names,testx,test_names,n_pc,n_df);
%
%       PC-DFA using training data (trainx)
%then projection of test data (testx) in PC-DFA space
%
%Inputs:
%  trainx :      training data
%  trainy :      list of classes from 1:n 
%  train_names : names for training data
%  testx :       test data to be projected
%  test_names :  names for test data
%  n_pc :        number of PCs to extract (min=3)
%  n_df :        number of DFs to extract (min=3, must be =< n_pc)
%
%Output as PCA and DFA plots - 1 vs 2; 1 vs 3; 2 vs 3.
%
% Copyright (c) 2002, Royston Goodacre
%

%PCA and PC-DFA
[tt,pp,pr]=pca_np(trainx,n_pc);
[U,V,eigenvals] = dfa(tt,trainy,n_df);

%PCA projection
m = mean(trainx);
[a,b] = size(testx);
Xmeansub =  testx - ones(a,1)*m;
projXtt=Xmeansub*pp';

%plotting PCA
pca_red(tt,1,2,train_names)
title('PCA - PC1 vs. PC2; red = train, blue = test');
hold on; p2d_col(projXtt,1,2,test_names,'b');

if n_pc>2
    pca_red(tt,1,3,train_names)
    title('PCA - PC1 vs. PC3; red = train, blue = test');
    hold on; p2d_col(projXtt,1,3,test_names,'b');

    pca_red(tt,2,3,train_names)
    title('PCA - PC2 vs. PC3; red = train, blue = test');
    hold on; p2d_col(projXtt,2,3,test_names,'b');
end
%PC-DFA projection
% projXU=projXtt*V*diag(eigenvals);
projXU=projXtt*V;

%plotting PC-DFA
dfa_red(U,1,2,train_names)
title('DFA - DF1 vs. DF2; red = train, blue = test');
hold on; p2d_col(projXU,1,2,test_names,'b')

if n_df>2
    dfa_red(U,1,3,train_names)
    title('DFA - DF1 vs. DF3; red = train, blue = test');
    hold on; p2d_col(projXU,1,3,test_names,'b')

    dfa_red(U,2,3,train_names)
    title('DFA - DF2 vs. DF3; red = train, blue = test');
    hold on; p2d_col(projXU,2,3,test_names,'b')
end

