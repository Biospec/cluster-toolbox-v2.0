function output=osc(Xcal,Ycal,nOSC, Xtest)
%Orthogonal Signal Correction using Fearn's method
%output=osc(Xcal,Ycal,nOSC, Xtest)
%Xcal, Ycal: the training data and concentration matrix
%nOSC: number of OSC components, default 1
%Xtest: the test data to be corrected
%reference:Tom Fearn, On orthogonal signal correction, 
%Chemometrics and Intelligent Laboratory Systems 50,2000. 47-52


if nargin<3
    nOSC=1;
end

if nargin<4
    Xtest=[];
end

[n,p]=size(Xcal);
D=eye(p);
SSQcal=sum(sum(Xcal.^2));
SSQtest=sum(sum(Xtest.^2));

% compute each orthogonal component and subtract them from X
for i=1:nOSC
    M=D-Xcal'*Ycal*pinv(Ycal'*Xcal*Xcal'*Ycal)*Ycal'*Xcal; 
    [ev,lambda]=eigs(Xcal*M*M'*Xcal',1); 
    w=1/sqrt(lambda)*M*Xcal'*ev; 
    t=Xcal*w;
    p=Xcal'*t/(t'*t);
    Xcal=Xcal-t*p';
    W(:,i)=w;
    T(:,i)=t;
    P(:,i)=p;   
end
Zcal=Xcal;


%ratio of explained variance of OSC components
R2Xcal=1-sum(sum(Zcal.^2))/SSQcal;
output.W_orth=W;
output.T_orth=T;
output.P_orth=P;
output.Xcal=Xcal;
output.R2Xcal=R2Xcal;

%Correcting new samples
if ~isempty(Xtest)
    for i=1:nOSC
     t=Xtest*W(:,i);
     Xtest=Xtest-t*P(:,i)';
    end
    Ztest=Xtest;
    output.Xtest=Ztest;
    R2Xtest=1-sum(sum(Ztest.^2))/SSQtest;
    output.R2Xtest=R2Xtest;
end

end