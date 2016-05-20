function [Tb,Pb,Wt,Tt,ssq,Rbo,Rbv,Lbo,Lbv,Lto] = cpca(X,Xin,nPC,tol)

% [Tb,Pb,Wt,Tt,ssq,Rbo,Rbv,Lbo,Lbv,Lto] = cpca(X,Xin,nPC,tol)
% X (objects x all variables) single, augmented data-block
% Xin (number of blocks x 2) = begin to end variable index for this block; index for X-block
% nPC (number of blocks x 1) number of principal components extracted per block
% tol (1 x 1) tolerance for convergence in NIPALS(default 1e-12)
%
% out :
% Tb (objects x number of blocks.max(nPC)) block scores, [t1-block-1 t1-block-2 ...t2-block-1...]
% Pb (all variables x max(nPC)) block loadings
% Wt (number of blocks x max(nPC)) super weights
% Tt (objects x max(nPC)) super scores
% ssq (max(nPC) x 1 + number of blocks) cumulative sum of squares for whole and individual blocks
% Rbo (objects x max(nPC)) object residuals
% Rbv (all variables x max(nPC)) block variable residuals
% Lbo (objects x max(nPC)) object leverages
% Lbv (all variables x max(nPC)) variable leverages
% Lto (objects x max(nPC)) super score object leverages

if nargin == 0
   help cpca
   return
elseif nargin == 3
   tol = 1e-12;
end

maxiter = 5000;
[n,m] = size(X);
nb = size(Xin,1);
X=X-ones(size(X,1),1)*mean(X);
[maxPC,maxb] = max(nPC);
% handling missing values, i.e. NaNs
MV = sum(sum(isnan(X)));
if MV
   Xmv = sparse(isnan(X));
   X(Xmv) = 0;
end
Tb = zeros(n,maxPC*nb);
Pb = zeros(m,maxPC);
Wt = zeros(nb,maxPC);
Tt = zeros(n,maxPC);
ssq = zeros(maxPC,1+nb);
Rbo = zeros(n,maxPC*nb);
Rbv = zeros(m,maxPC);
Lbo = zeros(n,maxPC*nb);
Lbv = zeros(m,maxPC);
Lto = zeros(n,maxPC);
ssqX(1) = sum(sum(X.^2));
for aa=1:nb
   rowi = Xin(aa,1):Xin(aa,2);
   ssqX(aa+1) = sum(sum(X(:,rowi).^2));
end
echo off;
opts.issym=1;
opts.isreal=1;
opts.tol=eps(1);
opts.maxit=300;
opts.p=size(X,1);
opts.v0=rand(size(X,1),1);
opts.disp=0;
opts.cholB=false;
opts.permB=1:size(X,1);
[v,d]=eigs(X*X',maxPC,'BE', opts);
for a=1:maxPC
   iter = 0;
%    Tt(:,a) = X(:,Xin(maxb,1)); % Old initialization method which could
%    results in differnt solutions when rearrange the order of samples.
   Tt(:,a)=v(:,a);   
   t_old = Tt(:,a)*100;
   if MV
      while (sum((t_old - Tt(:,a)).^2) > tol) & (iter < maxiter)
         iter = iter + 1;
         t_old = Tt(:,a);
         for aa=1:nb
            if nPC(aa) >= a
               rowi = Xin(aa,1):Xin(aa,2);
               coli = (a-1)*nb+aa;
               Pb(rowi,a) = X(:,rowi)'*Tt(:,a);
               % ---- Ignore missing values ----
               for aaa=rowi
                  Pb(aaa,a) = Pb(aaa,a)/(Tt(~Xmv(:,aaa),a)'*Tt(~Xmv(:,aaa),a));
               end
               % ----
               null_var=isnan(Pb(rowi,a));
               Pb(rowi,a) = Pb(rowi,a)/norm(Pb(rowi(~null_var),a));
               Tb(:,coli) = X(:,rowi(~null_var))*Pb(rowi(~null_var),a);
               % ---- Ignore missing values ----
               for aaa=1:n
                  Tb(aaa,coli) = Tb(aaa,coli)/(Pb(rowi(~null_var),a)'*Pb(rowi(~null_var),a));
               end
               % ----
            end
         end
         index = (a-1)*nb+1:a*nb;
         Wt(:,a) = Tb(:,index)'*Tt(:,a)/(Tt(:,a)'*Tt(:,a));
         Wt(:,a) = Wt(:,a)/norm(Wt(:,a));
         Tt(:,a) = Tb(:,index)*Wt(:,a)/(Wt(:,a)'*Wt(:,a));
      end
   else
      while (sum((t_old - Tt(:,a)).^2) > tol) & (iter < maxiter)
         iter = iter + 1;
         t_old = Tt(:,a);
         for aa=1:nb
            if nPC(aa) >= a
               rowi = Xin(aa,1):Xin(aa,2);
               coli = (a-1)*nb+aa;
               Pb(rowi,a) = X(:,rowi)'*Tt(:,a)/(Tt(:,a)'*Tt(:,a));
               Pb(rowi,a) = Pb(rowi,a)/norm(Pb(rowi,a));
               Tb(:,coli) = X(:,rowi)*Pb(rowi,a)/(Pb(rowi,a)'*Pb(rowi,a));
            end
         end
         index = (a-1)*nb+1:a*nb;
         Wt(:,a) = Tb(:,index)'*Tt(:,a)/(Tt(:,a)'*Tt(:,a));
         Wt(:,a) = Wt(:,a)/norm(Wt(:,a));
         Tt(:,a) = Tb(:,index)*Wt(:,a)/(Wt(:,a)'*Wt(:,a));
      end
   end
   if iter == maxiter
      s = ['WARNING: maximum number of iterations (' num2str(maxiter) ') reached before convergence'];
      disp(s)
   end
 
   for aa=1:nb
      if nPC(aa) >= a
         rowi = Xin(aa,1):Xin(aa,2);
         if MV
            Pb(rowi,a) = X(:,rowi)'*Tt(:,a);
            for aaa=rowi
               Pb(aaa,a) = Pb(aaa,a)/(Tt(~Xmv(:,aaa),a)'*Tt(~Xmv(:,aaa),a));
            end
         else
            Pb(rowi,a) = X(:,rowi)'*Tt(:,a)/(Tt(:,a)'*Tt(:,a));
         end
         X(:,rowi) = X(:,rowi) - Tt(:,a)*Pb(rowi,a)';
      end
   end
   
   if MV
      X(Xmv) = 0;
   end
   ssq(a,1) = (ssqX(1) - sum(sum(X.^2)))/ssqX(1);
   for aa=1:nb
      rowi = Xin(aa,1):Xin(aa,2);
      coli = (a-1)*nb+aa;
      ssq(a,aa+1) = (ssqX(aa+1) - sum(sum(X(:,rowi).^2)))/ssqX(aa+1);
      Rbo(:,coli) = sqrt(sum(X(:,rowi).^2,2));
      index = aa:nb:(a-1)*nb+aa;
      Lbo(:,coli) = diag(Tb(:,index)*pinv(Tb(:,index)'*Tb(:,index))*Tb(:,index)');      
   end
   Rbv(:,a) = sqrt(sum(X.^2,1))';
%    Lbv(:,a) = diag(Pb(:,1:a)*Pb(:,1:a)');
    Lbv(:,a)=NaN;
   Lto(:,a) = diag(Tt(:,1:a)*pinv(Tt(:,1:a)'*Tt(:,1:a))*Tt(:,1:a)');
end
echo on;