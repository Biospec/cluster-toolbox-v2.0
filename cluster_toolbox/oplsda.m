function  output = oplsda(X,Y,nOSC, nPLS, Xtest, Ytest)
% output = opls(X,Y,nOSC, nPLS, Xtest, Ytest)
% Orthogonal Projection to Latent Structure (OPLS)
% X: spectra matrix in training set
% Y: class membership matrix in training set
% nOSC: the number of orthogonal components, set to 1 by default
% nPLS: the number of PLS components, set to 1 by default
% Xtest: the spectra matrix in test set (optional)
% Ytest: the concentration matrix in test set (optional)
% output is a structure array with following fields:
% output
%       .TrainingSet: the centres of training set
%                .X_centre: the centre of X matrix
%                .Y_centre: the centre of Y matrix
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
%                .Yhat_cal: the numeric prediction of training set
%                .C_cal: class membership prediction of training set
%                .ccr_auto: correct classification rate of training set
%               
%       .Testset: the results of test set if Xtest is provided
%                .Xtest: the processed test data by OSC filter
%                .R2test: variation explayed by processed Xtest over
%                   original Xtest
%                .Yhat_test: the predicted concentration of test set
%                .C_test: class membership prediction of test set
%                .ccr_test: correct classification rate of test set, only
%                given if Ytest is also provided
% By Yun Xu, 2016

[m,n] = size(X);
[m2,c]=size(Y);
if m~=m2
    error('The number of rows in X and Y must be the same!')
end
if c==1
    class_train = Y;
    unique_cls = unique(Y);
    c = length(unique_cls);
    Ymat=zeros(m, c);
    for i=1:length(unique_cls)
        Ymat(Y==unique_cls(i),i)=1;
    end
    Y=Ymat;
else
    [~, class_train]=max(Y,[],2);
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

% build the OSC filter
if isempty(Xtest)
    OSC_model = osc(X, Y, nOSC);
else
    [m2,n2] = size(Xtest);
    if n2 ~= n
        error('The number of columns in Xtest and X must be the same!')
    end
    Xtest_mc = Xtest - repmat(X_centre, m2,1);  
    OSC_model = osc(X, Y, nOSC, Xtest);
end

% Output OSC model parameters
output.OSCmodel.W_orth=OSC_model.W_orth;
output.OSCmodel.T_orth=OSC_model.T_orth;
output.OSCmodel.P_orth=OSC_model.P_orth;
output.OSCmodel.Xcal=OSC_model.Xcal;
output.OSCmodel.R2Xcal=OSC_model.R2Xcal;

% Build PLS model on OSC corrected data
Xcal=OSC_model.Xcal;
[T,P,Q,W,b] = pls(Xcal, Y, nPLS);
[Yhat_cal,B] = plspred2(Xcal,P,Q,W,b,nPLS);
Yhat_cal = Yhat_cal + repmat(Y_centre,m,1);
[~, C_cal] = max(Yhat_cal,[],2);
% R2Y=1-sum(sum((Yhat_cal-Y).^2))/sum(sum((Y-repmat(mean(Y), size(Y,1),1)).^2)); 
% RMSEC=sqrt(sum((Yhat_cal-Y).^2)./size(Y,1));
ccr_auto = length(find(C_cal == class_train))/length(class_train);
output.PLSmodel.T=T;
output.PLSmodel.P=P;
output.PLSmodel.Q=Q;
output.PLSmodel.W=W;
output.PLSmodel.B=B;
output.PLSmodel.Yhat_cal = Yhat_cal;
output.PLSmodel.C_cal = C_cal;
% output.PLSmodel.R2=R2Y;
% output.PLSmodel.RMSEC=RMSEC;
output.PLSmodel.ccr_auto = ccr_auto;

%Correcting new samples if given
if ~isempty(Xtest)
    output.TestSet.Xtest=OSC_model.Xtest;
    Ztest=OSC_model.Xtest;
    R2Xtest=OSC_model.R2Xtest;
    Yhat_test = plspred2(Ztest,P,Q,W,b,nPLS);
    Yhat_test = Yhat_test + repmat(Y_centre,m2,1); 
    [~, C_test] = max(Yhat_test,[],2);    
    output.TestSet.Xtest=Ztest;
    output.TestSet.R2Xtest=R2Xtest;
    output.TestSet.Yhat_test = Yhat_test;    
    output.TestSet.C_test = C_test;        
end

if ~isempty(Ytest)
    [n,c]=size(Ytest);
    if c==1
        class_test=Ytest;
        unique_cls = unique(Ytest);
        c = length(unique_cls);
        Ymat=zeros(n, c);
        for i=1:length(unique_cls)
            Ymat(Ytest==unique_cls(i),i)=1;
        end
        Ytest=Ymat;
    else
        [~, class_test]=max(Ytest,[],2);
    end
    ccr_test = length(find(C_test==class_test))/length(class_test);
%     Q2=1-sum(sum((Yhat_test-Ytest).^2))./sum(sum((Ytest-repmat(mean(Ytest),...
%          size(Ytest,1),1)).^2));
%     RMSEP=sqrt(sum((Yhat_test-Ytest).^2)./size(Ytest,1));
%     output.TestSet.Q2p=Q2;
%     output.TestSet.RMSEP = RMSEP;
    output.TestSet.ccr_test = ccr_test;
end



end

