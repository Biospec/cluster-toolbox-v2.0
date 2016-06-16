function [alpha,testY,varargout] = dualpls(K,Ktest,Y,T,varargin)

%function [alpha,testY,varargout] = dualpls(K,Y,Ytest,T,varargin)
%
% Performs dual (kernel) PLS regression or discrimination
%
%INPUTS
% K = the training kernel matrix
% Ktest = the kernel matrix (dimension ell x elltest)
% Y = the training label matrix (ell x m)
% T = the number of PLS components to take
% varargin = optional argument specifying the true test label matrix
%            of size elltest x m, if 
%
%OUTPUTs
% alpha = the dual vectors corresponding to the PLS classfier
% testY = the estimated label matrix for the test samples
% varargout = the test error, optional when varargin is specified (the
%             true test labels)
%
%Note: this function is mostly based on the one published on 
%      www.kernel-methods.net with modifications to prevent undesired
%      behaviours (infinit loop, rank deficiency) from happening.

% K is an ell x ell kernel matrix
% Y is ell x m containing the corresponding output vectors
% T gives the number of iterations to be performed
ell=size(K,1);
trainY=0;
KK = K; YY = Y;
for i=1:T
    YYK = YY*YY'*KK;
    beta(:,i) = YY(:,1)/norm(YY(:,1));
    if size(YY,2) > 1, % only loop if dimension greater than 1
        bold = beta(:,i) + 1;
        count=1; % a counter to prevent infinit loop
        while norm(beta(:,i) - bold) > 0.001,
            bold = beta(:,i);
            tbeta = YYK*beta(:,i);
            beta(:,i) = tbeta/norm(tbeta);
            count=count+1;
            if count>100
                disp('The algorithm failed to converge, skip...')
                break
            end
        end
    end
    tau(:,i) = KK*beta(:,i);
    val = tau(:,i)'*tau(:,i);
    c(:,i) = YY'*tau(:,i)/val;
    trainY = trainY + tau(:,i)*c(:,i)';
    trainerror = norm(Y - trainY,'fro')/sqrt(ell);
    w = KK*tau(:,i)/val;
    KK = KK - tau(:,i)*w' - w*tau(:,i)' + tau(:,i)*tau(:,i)'*(tau(:,i)'*w)/val;
    YY = YY - tau(:,i)*c(:,i)';
end

% Avoid extracting too many components

rnk = rank(tau'*K*beta);
if rnk < size(beta,2) 
    disp(['NULL components detected, T might be too high, reduce it to' ...
        num2str(rnk)]);
    beta = beta(:, rnk);
    tau = tau(:, rnk);
end
% Regression coefficients for new data
alpha = beta * ((tau'*K*beta)\tau')*Y;

%  Ktest gives new data inner products as rows, Ytest true outputs
elltest = size(Ktest,1);
testY = Ktest * alpha;
if ~isempty(varargin)
    Ytest = varargin{1};
    testerror = {norm(Ytest - testY,'fro')/sqrt(elltest)};
    varargout = testerror;
end
