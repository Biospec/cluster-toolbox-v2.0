function output = plsr_boots(data, conc, rep, no_pcs, no_loops)
% output = plsr_boots(data, conc, rep, no_pcs, no_loops)
% data: pre-processed data
% conc: the known concentrations
% rep: an index vector indicating which samples are from the same bio-rep.
%      all the reps of the same biological sample shall have the same
%      number. e.g. the rep_idx of an experiment having 5 biological samples 
%      with 3 analytical reps each should be [1;1;1;2;2;2;3;3;3;4;4;4;5;5;5]. 
%      Set to 1:no. of samples if not specified, i.e. assuming each bio-rep
%      has only one sample.
% no_pcs: the maximum number of PLS component to be extracted, the optimal 
%         number of PCs will be determined by cross-validation, set to 15
%         if not specified.
%no_loops: the number of iterations of bootstrapping resampling, set to
%       1000 if not specified.
% output
%   .Cknown: the know concentrations of each bootstrapping test set
%   .Cknown_av: averaged known concetrations of all the bootstrapping test sets
%   .Cpred: the predicted concentrations of each bootstrapping test sets
%   .Cpred_av: averaged predicted concetrations of all the bootstrapping test sets 
%   .Cpred_std: the std of the predicted concentrations of all the bootstrapping test sets 
%   .PLS_loadings: The PLS loadings
%   .Q2_boots: cross-validated R2 of all the bootstrapping training sets,
%              used for selecting the optimal PLS components.
%   .Q2p: R2 of all the bootstrapping test sets.
%   .no_pcs: the optimal PLS components selected for each bootstrapping
%            iteration.
%   .r2: R2 of the full data set, using the minimal no. of PLS components
%       used by bootstrapping iterations.
%   .rmsec: Root mean square error of the full data set, using the minimal no. 
%          of PLS components used by bootstrapping iterations.
%   .rmsecv_boots: Root mean square error of cross-validation of each of the
%          bootstrapping training set.
%   .rmsep: Root mean square error of the test set of each of the
%          bootstrapping test set.
%   .vip_scores: Variable Importance in Projection for each of the
%          variables.


%Last update: 28/09/2015 - RMSEC and R2 added Yun Xu

    if nargin<2
        help plsr_boots;
        output=NaN;
        return;
    end
    if nargin==2
        rep=[1:size(data,1)]';
        no_pcs=15;
        no_loops=1000;
    end
    if nargin==3
        no_pcs=15;
        no_loops=1000;
    end
    if nargin==4
        no_loops=1000;
    end

    unique_rep=unique(rep);
    no_reps=length(unique_rep);
    h=waitbar(0, 'Please wait...');
    [m,c]=size(conc);
    if c==1
        conc=conc(:);
        unique_conc=unique(conc);
        conc_pred_full=nan(no_loops, length(unique_conc));
    else
        unique_conc=unique(conc,'rows');
        [m,n]=size(unique_conc);
        conc_pred_full=nan(m,n,no_loops);
    end

    for i=1:no_loops
        trn_rep=unique_rep(unique(randi(no_reps,no_reps,1)));
        while length(trn_rep)<4
            trn_rep=unique_rep(unique(randi(no_reps,no_reps,1)));
        end
        trn_idx=find(ismember(rep, trn_rep));
        data_trn=data(trn_idx,:);
        data_tst=data; data_tst(trn_idx,:)=[];
        conc_trn=conc(trn_idx,:);
        conc_tst=conc; conc_tst(trn_idx,:)=[];
        amean=mean(data_trn); 
        cmean=mean(conc_trn);
        no_rep_trn=length(trn_rep);
        trn_rep_idx=rep(trn_idx);
        if no_rep_trn<14
            step=1;
        else
            step=round(no_rep_trn/7);
        end

        for a=1:no_pcs
            predCval=[];Cval=[];
            for ii=1:step:no_rep_trn
                rep_val=trn_rep(ii:min(ii+step-1, no_rep_trn));
                val_idx=find(ismember(trn_rep_idx, rep_val));
                data_val=data_trn(val_idx,:); conc_val=conc_trn(val_idx,:);
                data_trn_inner=data_trn; data_trn_inner(val_idx,:)=[];
                conc_trn_inner=conc_trn; conc_trn_inner(val_idx,:)=[];
                amean_inner=mean(data_trn_inner); 
                cmean_inner=mean(conc_trn_inner);
                data_trn_inner=data_trn_inner-repmat(amean_inner, size(data_trn_inner,1),1);
                conc_trn_inner=conc_trn_inner-repmat(cmean_inner, size(conc_trn_inner,1),1);
                [T,P,Q,W,b]=pls(data_trn_inner,conc_trn_inner,a);
                predC = plspred2(data_val,P,Q,W,b,a,amean_inner, cmean_inner);
                predCval=[predCval; predC]; Cval=[Cval;conc_val];
            end
            rmsecv(a,:)=sqrt(sum((Cval-predCval).^2)/size(Cval,1));
            Q2(a,:)=1-sum((Cval-predCval).^2)./sum((Cval-repmat(mean(Cval),size(Cval,1),1)).^2);
            rmsecv_average(a)=mean(rmsecv(a,:));
            Q2_average(a)=mean(Q2(a,:));
        end
        [opt_Q2, opt_pc]=min(rmsecv_average);
        rmsecv_boots(i,:)=rmsecv(opt_pc,:);
        Q2_boots(i,:)=Q2(opt_pc,:);
        clear rmsecv Q2 rmsecv_average Q2_average
        data_trn=data_trn-repmat(amean, size(data_trn,1),1);
        conc_trn=conc_trn-repmat(cmean, size(conc_trn,1),1);
        [T,P,Q,W,b]=pls(data_trn,conc_trn,opt_pc);
        predC = plspred2(data_tst,P,Q,W,b,opt_pc,amean, cmean);
        rmsep=sqrt(sum((predC-conc_tst).^2)/size(predC,1));
        Q2p = 1-sum((predC-conc_tst).^2)./sum((conc_tst-repmat(mean(conc_tst),size(conc_tst,1),1)).^2);

%        calculating R2 and rmsec for each resampling
%        predCtrn = plspred2(data_trn,P,Q,W,b,opt_pc);
%        predCtrn = predCtrn+repmat(cmean,size(data_trn,1),1);
%        conc_trn=conc_trn+repmat(cmean, size(conc_trn,1),1);
%        rmsec(i,:)=sqrt(sum(predCtrn-conc_trn).^2/size(predCtrn,1));
%        r2(i,:) = 1-sum((predCtrn-conc_trn).^2)./sum((conc_trn-repmat(mean(conc_trn), size(conc_trn,1),1)).^2);

    %     conc_full=[conc_full;conc_tst];
    %     conc_pred_full=[conc_pred_full; predC];
        if c==1
            for ii=1:length(unique_conc)
                conc_pred_full(i, ii)=mean(predC(conc_tst==unique_conc(ii)));
            end
        else
            for ii=1:size(unique_conc,1)
                conc_pred_full(ii,:,i) = mean(predC(ismember(conc_tst, unique_conc(ii,:),'rows'),:));
            end
        end 

        output.Cknown{i}=conc_tst;
        output.Cpred{i}=predC;
        output.rmsep(i,:)=rmsep;
        output.Q2p(i,:)=Q2p;
        output.no_pcs(i)=opt_pc;
        waitbar(i/no_loops,h);
    end

    if c==1
        Cknown_av=unique_conc;
        pred_Cav=nanmean(conc_pred_full); pred_Cav=pred_Cav(:);
        pred_Cstd=nanstd(conc_pred_full); pred_Cstd=pred_Cstd(:);
    else
        Cknown_av=unique_conc;
        pred_Cav=nanmean(conc_pred_full,3);
        pred_Cstd=nanstd(conc_pred_full,0,3);
    end


    close(h);
    output.rmsecv_boots=rmsecv_boots;
    output.Q2_boots=Q2_boots;
    output.Cknown_av=Cknown_av;
    output.Cpred_av=pred_Cav;
    output.Cpred_std=pred_Cstd;

    figure
    if c==1
        plot(Cknown_av, pred_Cav,'^')
        hold on
        h=errorbar(Cknown_av, pred_Cav, pred_Cstd);
        set(h,'linestyle','none');
        b=pinv([ones(size(pred_Cav,1),1) pred_Cav])*Cknown_av;
        set(gca,'XLim',[min(Cknown_av)-mean(Cknown_av)*.05 max(Cknown_av)+mean(Cknown_av)*.05]);
        fn = @(x)[ones(size(x,1),1) x]*b;
        fplot(fn,get(gca,'XLim'));
        h=xlabel('Known concentration'); set(h,'fontsize',14);
        h=ylabel('Predicted concentration'); set(h,'fontsize',14);
        h=title('Averaged known {\itvs} predicted'); set(h,'fontsize',14);
    else
        for i=1:c
            subplot(c,1,i);
            plot(Cknown_av(:,i), pred_Cav(:,i),'^')
            hold on
            h=errorbar(Cknown_av(:,i), pred_Cav(:,i), pred_Cstd(:,i));
            set(h,'linestyle','none');   
            b=pinv([ones(size(pred_Cav(:,i),1),1) pred_Cav(:,i)])*Cknown_av(:,i);
            fn = @(x)[ones(size(x,1),1) x]*b;
            fplot(fn,[min(Cknown_av(:,i)) max(Cknown_av(:,i))]);
            h=xlabel('Known concentration'); set(h,'fontsize',14);
            h=ylabel('Predicted concentration'); set(h,'fontsize',14);
            h=title(['Averaged known {\itvs} predicted, component ' num2str(i)]); 
            set(h,'fontsize',14);
        end
    end
% To minimise the risk of over-fitting, use the smallest no. of
% components across all bootstrapping resamplings for calculating 
% auto-prediction statistics, such as VIP, Loadings, R2 and RMSEC
    no_pcs_final=min(output.no_pcs);
    amean=mean(data); 
    cmean=mean(conc);
    data_mc=data-repmat(amean, size(data,1),1);
    conc_mc=conc-repmat(cmean, size(conc,1),1);
    [T,P,Q,W,b]=pls(data_mc,conc_mc,no_pcs_final);
    [predCtrn,B] = plspred2(data,P,Q,W,b,no_pcs_final,amean, cmean);
    r2 = 1 - sum((predCtrn-conc).^2)./sum((conc-repmat(mean(conc),size(conc,1),1)).^2);
    rmsec = sqrt(sum((predCtrn-conc).^2)/size(conc,1));
    vip_scores = vip(T,P,W,B);
    output.r2=r2;
    output.rmsec = rmsec;
    output.vip_scores=vip_scores;
    output.PLS_loadings=P;    
    output = orderfields(output);
end

function [T,P,Q,W,b,X,Y]=pls(A,C,maxrank)
% function [T,P,Q,W,b,X,Y]=pls(A,C,maxrank);
%A: spectral matrix in training set
%C: concentration matrix in training set
%maxrank: the number of components to be extracted.
%T: scores of spectral matrix
%P: loadings of spectral matrix
%Q: loadings of concentration matrix
%W: weight matrix
%b: coefficient
%X: Residual spectral matrix
%Y: Residual concentration matrix
%This program is for general use, eithe PLS1 or PLS2.


[m,n]=size(A);
maxrank=min([m n maxrank]);

X=A;
% X=[A; diag(ones(n,1)*10)];
Y=C;
% Y=[Y; zeros(n,size(C,2))];
for i=1:maxrank
    XY=X'*Y;
    [u,s,v]=svd(XY,0);
    q=v(:,1);
    w=u(:,1);
    t=X*w;
    p=t'*X/(t'*t);p=p';
    u=Y*q;
    b(i,1)=1/(t'*t)*t'*u;
    X=X-t*p';
    Y=Y-b(i)*t*q';
    
    T(:,i)=t;
    W(:,i)=w;
    Q(:,i)=q;
    P(:,i)=p;
    
end
end

function [C,BV]=plspred2(X,P,Q,W,b,N,amean,cmean,ascal,cscal)
%function [C,BV]=plspred2(X,P,Q,W,b,n,amean, cmean, ascal, cscal);
%Same to plspred. But quick prediction procedure is used.
%please refer to plspred
[m,n]=size(X);
if nargin<5
    help plspred;
    return;
elseif nargin==5
    N=size(P,2);
elseif nargin==8
%     [m,n]=size(X);
    X=X-ones(m,1)*amean;
elseif nargin>8
    X=X-ones(m,1)*amean;
    X=X./(ones(m,1)*ascal);
end
W=W(:,1:min(N, size(W,2)));
P=P(:,1:min(N, size(W,2)));
b=b(1:min(N, size(W,2)));
Q=Q(:,1:min(N, size(W,2)));
R=Q*diag(b);

B=W*inv(P'*W)*R';
C=X*B;
BV=B;


%  [m,n]=size(X);

if nargin==8
    C=C+ones(m,1)*cmean;
end   
if nargin==10
    C=C+ones(m,1)*cmean;
    C=C.*(ones(m,1)*cscal);
end
end


function vip_scores = vip(T,P,w,B)
%   vip_scores = vip(T,P,w,B)
%   T = X-block scores
%    P = X-block loadings
%     w = X-block weights
%     B = regression vectors for each column of y and each number of
%           latent variables (reg) 
%
% OUTPUTS: 
%  vip_scores = a set of column vectors equal in length to the number of
%  variables included in the model. It contains one column of VIP scores
%  for each predicted y-block column.
%
% See Chong & Jun, Chemo. Intell. Lab. Sys. 78 (2005) 103?12.
%


wpw  = (w/(P'*w));
nx   = size(T,1);
ny   = size(B,2);

%pre-calculate some misc. things
TT = sum(T.^2,1);
w_norm = (w*diag(1./sqrt(sum(w.^2,1))));  %normalized weights

for i = 1:ny;
  %calculate regression in terms of scores (T*b = y_hat)
  b  = wpw\B(:,i);

  %calculate weighted T^2
  SS = b.^2.*TT';

  %VIP scores for this y
  vip_scores(:,i) = nx*w_norm.^2*SS./sum(SS);
end
end
