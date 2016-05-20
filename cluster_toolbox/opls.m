function  output = opls(X,Y,nOSC, nPLS, Xtest, Ytest)
% output = opls(X,Y,nOSC, nPLS, Xtest, Ytest)
% Orthogonal Projection to Latent Structure (OPLS) for PLS 1 (Y is a vector)
% X: spectra matrix in training set
% Y: concentration vector in training set
% nOSC: the number of orthogonal components, set to 1 by default
% nPLS: the number of PLS components, set to 1 by default
% Xtest: the spectra matrix in test set (optional)
% Ytest: the concentration matrix in test set (optional)
% output is a structure array with following fields:
% output
%       .OSCmodel: the parameters of the OSC filter
%                .W_orth: orthogonal weight 
%                .T_orth: orthogonal scores
%                .P_orth: orthogonal loadings
%                .Xcal: the processed training data by OSC filter
%                .R2Xcal: variation explained by Xcal over X
%       .PLSmodel: the parameters of the PLS model
%                .T: PLS X scores
%                .P: PLS X loadings
%                .Q: PLS Y loadings
%                .W: PLS weight matrix
%                .B: PLS regression coefficient
%                .Yhat_cal: the predicted concentration of training set
%                .R2: squared correlation coefficients of training set
%                .RMSEC: root mean square error of calibration (training
%                   set)
%       .Testset: the results of test set if Xtest is provided
%                .Xtest: the processed test data by OSC filter
%                .R2test: variation explayed by processed Xtest over
%                   original Xtest
%                .Yhat_test: the predicted concentration of test set
%        the following fileds will be given if Ytest is also provided
%                .Q2p: squared correlation coefficients of test set
%                .RMSEP: root mean square error of prediction (test set)
% By Yun Xu, 2016

[m,n] = size(X);
[m2,c] = size(Y);

if m ~=m2
    error('X and Y need to have same number of rows!')
end

if c>1
    error('Y need to be a single column vector, if Y has more than one component, use opls2.m instead!')
end

if nargin<3
    nOSC=1;
    nPLS=1;
end

if nargin<4
    nPLS=1;
end

if nargin<5
    Xtest=[];
    Ytest=[];
end

if nargin<6
    Ytest=[];
end
% mean centre pre-processing
X_centre = mean(X);
Y_centre = mean(Y);
X_mc = X - repmat(X_centre, m,1);
Y_mc = Y - repmat(Y_centre, m,1);
output.TrainingSet.X_centre = X_centre;
output.TrainingSet.Y_centre = Y_centre;

% Orthogonal projection to latent structure
ssqXcal=sum(sum((X_mc.^2))); 
w=X_mc'*Y_mc;         % calculate weight vector
w=w/norm(w);    % normalization
for i=1:nOSC 
    t=X_mc*w;                     % calculate scores vector
    p=X_mc'*t/(t'*t);             % calculate loadings of X
    wosc=p-(w'*p)/(w'*w)*w;    % orthogonal weight
    wosc=wosc/norm(wosc);      % normalization
    tosc=X_mc*wosc;               % orthogonal components
    posc=X_mc'*tosc/(tosc'*tosc); % loadings 
    
    X_mc=X_mc-tosc*posc';            % remove orthogonal components
    W_orth(:,i)=wosc;          % record results
    T_orth(:,i)=tosc;          % record results
    P_orth(:,i)=posc;          % record results
    
end
% Output orthogonal projection model parameters
output.OSCmodel.W_orth=W_orth;
output.OSCmodel.T_orth=T_orth;
output.OSCmodel.P_orth=P_orth;
Zcal=X_mc;  % after subtracting all orthogonal components, X_mc is the new X for PLS 
output.OSCmodel.Xcal=Zcal;
R2Xcal=1-sum(sum(Zcal.^2))/ssqXcal; 
output.OSCmodel.R2Xcal=R2Xcal;

% PLS model
[T,P,Q,W,b] = pls(Zcal, Y_mc, nPLS);
[Yhat_cal,B] = plspred2(Zcal,P,Q,W,b,nPLS);
Yhat_cal = Yhat_cal + repmat(Y_centre,size(Y,1),1);
R2=1-sum((Yhat_cal-Y).^2)./sum((Y-repmat(mean(Y), size(Y,1),1)).^2);
RMSEC=sqrt(sum((Yhat_cal-Y).^2)./size(Y,1));
output.PLSmodel.T=T;
output.PLSmodel.P=P;
output.PLSmodel.Q=Q;
output.PLSmodel.W=W;
output.PLSmodel.B=B;
output.PLSmodel.Yhat_cal = Yhat_cal;
output.PLSmodel.R2Y=R2;
output.PLSmodel.RMSEC=RMSEC;

%Correcting new samples if given
if ~isempty(Xtest)
    [m2,n2]=size(Xtest);
    Xtest_mc = Xtest - repmat(X_centre, m2,1);
    ssqXtest=sum(sum((Xtest_mc.^2))); 
    for i=1:nOSC
     t=Xtest_mc*W_orth(:,i);
     Xtest_mc=Xtest_mc-t*P_orth(:,i)';
    end
    Ztest=Xtest_mc;
    R2Xtest=1-sum(sum(Ztest.^2))/ssqXtest;
    output.TestSet.Xtest=Ztest;
    output.TestSet.R2Xtest=R2Xtest;
    Yhat_test = plspred2(Ztest,P,Q,W,b,nPLS);
    Yhat_test = Yhat_test + repmat(Y_centre,m2,1); 
    output.TestSet.Yhat_test = Yhat_test;
end

if ~isempty(Ytest)
    Q2=1-sum((Yhat_test-Ytest).^2)./sum((Ytest-repmat(mean(Ytest),...
        size(Ytest,1),1)).^2);
    RMSEP=sqrt(sum((Yhat_test-Ytest).^2)./size(Ytest,1));    
    output.TestSet.Q2p=Q2;
    output.TestSet.RMSEP = RMSEP;
end



end
