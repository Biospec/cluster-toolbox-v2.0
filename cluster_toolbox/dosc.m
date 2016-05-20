function output = dosc(Xcal,Ycal,nOSC,Xtest)
% DOSC Direct Orthogonal Signal Correction
%
% output = dosc(Xcal,Ycal,nOSC)
%
% Input
% Xcal  training data matrix
% Ycal  training concentration matrix
% nOSC  number of DOSC components to calculate 
% Xtest test data matrix (optional)
%
% output.
%       .Xcal   DOSC corrected training data matrix 
%       .W      weights used to determine DOSC components
%       .P      loadings used to remove DOSC component from Xcal
%       .T      scores of DOSC components 
%       .Xtest  DOSC corrected test data matrix if provided
% 
% 
% Once the calibration is done, new (scaled) x-data can be corrected by 
% newx = x - x*W*P'; Or use mfile dosc_pred
%
% See Reference
% Westerhuis JA, de Jong S and Smilde AK, Direct orthogonal signal correction, 
% Chemometrics and Intelligent Laboratory Systems, 56, (2001), 13-25.

if nargin<4
    Xtest=[];
end
% project Y onto X (step 1)
Yhat = Xcal*(pinv(Xcal')'*Ycal);

% deflate X wrt Yhat (step 2)
AyX = Xcal-Yhat*(pinv(Yhat)*Xcal); 

% find major PCs of AyX
[Ta,D] = eigs(AyX*AyX',nOSC);

[U,S,V] = svd(Xcal,0);
S_ratio = diag(S)./sum(diag(S));
idx = find(cumsum(S_ratio)>.9999);
Xcorr = U(:,1:idx)*S(1:idx,1:idx)*V(:,1:idx)';
pinvX = pinv(Xcorr);

W = pinvX*Ta;
T = Xcal*W;

% Calculate loadings to remove DOSC component (step 7a)
% P = Xcal'*T*inv(T'*T);
P = (Xcal'*T)/(T'*T);
% deflate X wrt to DOSC components (step 6a)
Z = Xcal - T*P';

output.Xcal = Z;
output.W_orth = W;
output.T_orth = T;
output.P_orth = P;

% apply to new samples if provided
if ~isempty(Xtest)
    Ztest = Xtest - Xtest*W*P';
    output.Xtest = Ztest;
end


