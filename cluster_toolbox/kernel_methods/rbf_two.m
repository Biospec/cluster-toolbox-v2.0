function K = rbf_two(X, Y, sig)
%function K = rbf_two(X,Y, sig)
%
% Computes an rbf kernel matrix from the two matrices
%
%INPUTS
% X =  sample matrix 1 containing all samples as rows
% Y =  sample matrix 2 containing all samples as rows
% sig = sigma, the kernel width; squared distances are divided by
%       squared sig in the exponent
%
%OUTPUTS
% K = the rbf kernel matrix 
%
%
% For more info, see www.kernel-methods.net

%
% New code by Yun (10/06/2016)

ssqX = sum(X.^2,2);
nX = size(X,1);
ssqY = sum(Y.^2,2);
nY = size(Y,1);
D = ones(nY, 1)*ssqX' + (ones(nX,1)*ssqY')' - 2*Y*X';
K = exp(-D/(2*sig^2));

%% Slow version 
% n1 = size(X,1);
% n2 = size(Y,1);
% K=zeros(n1,n2);
% for i=1:n1
%     for j=1:n2
%         K(i,j)=exp(-norm(X(i,:)-Y(j,:))^2/(2*sig^2));
%        %K(i,j)=exp(-norm(X(i,:)-Y(j,:))^2/(2*sig));
%     end
% end