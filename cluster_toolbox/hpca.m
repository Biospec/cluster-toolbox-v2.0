function [Tb,Pb,Wt,Tt,ssq,Rbo,Rbv,Lbo,Lbv,Lto] = hpca(X,Xin,nPC,tol)
% function [Tb,Pb,Wt,Tt,ssq,Rbo,Rbv,Lbo,Lbv,Lto] = hpca(X,Xin,nPC,tol)
% Multiblock Hierachical PCA
% in  : X (objects x all variables) single, augmented data-block
%       Xin (number of blocks x 2) = begin to end variable index for this block; index for X-block
%       nPC (number of blocks x 1) number of latent variables extracted per block
%       tol (1 x 1) tolerance for convergence (default 1e-12)
%
% out : Tb (objects x number of blocks.max(nPC)) block scores, [t1-block-1 t1-block-2 ...t2-block-1...]
%       Pb (all variables x max(nPC)) block loadings
%       Wt (number of blocks x max(nPC)) super weights
%       Tt (objects x max(nPC)) super scores
%       ssq (max(nPC) x number of blocks) cumulative sum of squares
%       Rbo (objects x max(nPC)) object residuals
%       Rbv (all variables x max(nPC)) block variable residuals
%       Lbo (objects x max(nPC)) object leverages
%       Lbv (all variables x max(nPC)) variable leverages
%		  Lto (objects x max(nPC)) super score object leverages

if nargin == 0
   help hpca
   return
elseif nargin == 3
   tol = 1e-12;
end

maxiter = 2000;
[n,m] = size(X);
nb = size(Xin,1);
[maxPC,maxb] = max(nPC);
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

for a=1:maxPC
   iter = 0;
   Xinit = [];
   for aa=1:nb
      if nPC(aa) >= a
         rowi = Xin(aa,1):Xin(aa,2);
         Xinit = [Xinit X(:,rowi)];
      end
   end
   [v,d] = eigs(Xinit*Xinit',1);
   Tt(:,a) = v;
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
               for aaa=rowi
                  Pb(aaa,a) = Pb(aaa,a)/(Tt(~Xmv(:,aaa),a)'*Tt(~Xmv(:,aaa),a));
               end
               Tb(:,coli) = X(:,rowi)*Pb(rowi,a);
               for aaa=1:n
                  Tb(aaa,coli) = Tb(aaa,coli)/(Pb(~Xmv(aaa,rowi),a)'*Pb(~Xmv(aaa,rowi),a));
               end
               Tb(:,coli) = Tb(:,coli)/norm(Tb(:,coli));
            end
         end
         index = (a-1)*nb+1:a*nb;
         Wt(:,a) = Tb(:,index)'*Tt(:,a)/(Tt(:,a)'*Tt(:,a));
         Tt(:,a) = Tb(:,index)*Wt(:,a)/(Wt(:,a)'*Wt(:,a));
         Tt(:,a) = Tt(:,a)/norm(Tt(:,a));
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
               Tb(:,coli) = X(:,rowi)*Pb(rowi,a)/(Pb(:,a)'*Pb(:,a));
               Tb(:,coli) = Tb(:,coli)/norm(Tb(:,coli));
            end
         end
         index = (a-1)*nb+1:a*nb;
         Wt(:,a) = Tb(:,index)'*Tt(:,a)/(Tt(:,a)'*Tt(:,a));
         Tt(:,a) = Tb(:,index)*Wt(:,a)/(Wt(:,a)'*Wt(:,a));
         Tt(:,a) = Tt(:,a)/norm(Tt(:,a));
      end
   end
   if iter == maxiter
      s = ['WARNING: maximum number of iterations (' num2str(maxiter) ') reached before convergence'];
      disp(s)
   end
   
   X = X - Tt(:,a)*Pb(:,a)';
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
   Lbv(:,a) = diag(Pb(:,1:a)*Pb(:,1:a)');
   Lto(:,a) = diag(Tt(:,1:a)*pinv(Tt(:,1:a)'*Tt(:,1:a))*Tt(:,1:a)');
end