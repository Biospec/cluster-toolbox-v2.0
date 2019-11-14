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


%Last update: 26/07/2016 - permutation test added Yun Xu

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
        perm_idx=reps_perm(rep);
%         perm_idx = randperm(size(data,1));
        data2 = data(perm_idx,:);
        trn_rep=unique_rep(unique(randi(no_reps,fix(no_reps*1),1)));
        while length(trn_rep)<=4 || length(trn_rep)>=no_reps-2
            trn_rep=unique_rep(unique(randi(no_reps,fix(no_reps*1),1)));
        end
        trn_idx=find(ismember(rep, trn_rep));
        data_trn=data(trn_idx,:);
        data_trn2=data2(trn_idx,:);
        data_tst=data; data_tst(trn_idx,:)=[];
        data_tst2=data2; data_tst2(trn_idx,:)=[];
        conc_trn=conc(trn_idx,:); 
        conc_tst=conc; conc_tst(trn_idx,:)=[];
        amean=mean(data_trn); 
        amean2=mean(data_trn2);
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
        [T2,P2,q2,W2,b2]=pls(data_trn2, conc_trn, opt_pc);        
        [predC, B] = plspred2(data_tst,P,Q,W,b,opt_pc,amean, cmean);
        predC2 = plspred2(data_tst2,P2,q2,W2,b2,opt_pc,amean2,cmean);
        rmsep=sqrt(sum((predC-conc_tst).^2)/size(predC,1));
        rmsep2=sqrt(sum((predC2-conc_tst).^2)/size(predC2,1));
        Q2p = 1-sum((predC-conc_tst).^2)./sum((conc_tst-...
            repmat(mean(conc_tst),size(conc_tst,1),1)).^2);
        Q2p_perm = 1-sum((predC2-conc_tst).^2)./sum((conc_tst-...
            repmat(mean(conc_tst),size(conc_tst,1),1)).^2);
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
                conc_pred_full2(i,ii)=mean(predC2(conc_tst==unique_conc(ii)));
            end
        else
            for ii=1:size(unique_conc,1)
                conc_pred_full(ii,:,i) = mean(predC(ismember(conc_tst, unique_conc(ii,:),'rows'),:));
                conc_pred_full2(ii,:,i) = mean(predC2(ismember(conc_tst, unique_conc(ii,:),'rows'),:));
            end
        end 

        output.Cknown{i}=conc_tst;
        output.Cpred{i}=predC;
        output.rmsep(i,:)=rmsep;
        output.Q2p(i,:)=Q2p;
        output.Q2p_perm(i,:)=Q2p_perm;
        output.no_pcs(i)=opt_pc;
        if c==1
            output.reg_coeff(i,:) = B;
        else
            output.reg_coeff{i} = B;
        end
        waitbar(i/no_loops,h);
    end

    if c==1
        Cknown_av=unique_conc;
        pred_Cav=nanmean(conc_pred_full); pred_Cav=pred_Cav(:);
        pred_C2av=nanmean(conc_pred_full2); pred_C2av=pred_C2av(:);
        pred_Cstd=nanstd(conc_pred_full); pred_Cstd=pred_Cstd(:);
        pred_C2std=nanstd(conc_pred_full2); pred_C2std=pred_C2std(:);
    else
        Cknown_av=unique_conc;
        pred_Cav=nanmean(conc_pred_full,3);
        pred_C2av=nanmean(conc_pred_full2,3);
        pred_Cstd=nanstd(conc_pred_full,0,3);
        pred_C2std=nanstd(conc_pred_full2,0,3);
    end


    close(h);
    if c==1
        pval = length(find(output.Q2p_perm > output.Q2p))/no_loops;
    else
        pval = length(find(mean(output.Q2p_perm,2) > mean(output.Q2p,2)))...
            /no_loops;
    end
    output.pval = pval;
    output.rmsecv_boots=rmsecv_boots;
    output.Q2_boots=Q2_boots;
    output.Cknown_av=Cknown_av;
    output.Cpred_av=pred_Cav;
    output.Cpred_std=pred_Cstd;
    output.Cpred_av_perm=pred_C2av;
    output.Cpred_std_perm=pred_C2std;
    figure
    if c==1
        plot(Cknown_av, pred_Cav,'^')
        hold on
        h=errorbar(Cknown_av, pred_Cav, pred_Cstd);
        set(h,'linestyle','none');
        b=pinv([ones(size(pred_Cav,1),1) pred_Cav])*Cknown_av;
        set(gca,'XLim',[min(Cknown_av)-mean(Cknown_av)*.05 max(Cknown_av)+mean(Cknown_av)*.05]);
        h=line([min(Cknown_av) max(Cknown_av)], [min(Cknown_av) max(Cknown_av)]);
        set(h,'linestyle', '--');
%        fn = @(x)[ones(size(x,1),1) x]*b;
        fn = @(x)b(1) + x*b(2);
        fplot(fn,get(gca,'XLim'));
        h=xlabel('Known concentration'); set(h,'fontsize',14);
        h=ylabel('Predicted concentration'); set(h,'fontsize',14);
        h=title('Averaged known {\itvs} predicted (Test set)'); set(h,'fontsize',14);
    else
        for i=1:c
            subplot(c,1,i);
            plot(Cknown_av(:,i), pred_Cav(:,i),'^')
            hold on
            h=errorbar(Cknown_av(:,i), pred_Cav(:,i), pred_Cstd(:,i));
            set(h,'linestyle','none'); 
            h=line([min(Cknown_av) max(CKnown_av)], [min(Cknown_av) max(Cknown_av)]);
            set(h,'linestyle', '--');
            b=pinv([ones(size(pred_Cav(:,i),1),1) pred_Cav(:,i)])*Cknown_av(:,i);
            fn = @(x)[ones(size(x,1),1) x]*b;
            fplot(fn,[min(Cknown_av(:,i)) max(Cknown_av(:,i))]);
            h=xlabel('Known concentration'); set(h,'fontsize',14);
            h=ylabel('Predicted concentration'); set(h,'fontsize',14);
            h=title(['Averaged known {\itvs} predicted, component ' num2str(i) ' (Test set)']); 
            set(h,'fontsize',14);
        end
    end
%     figure
%     if c==1
%         plot(Cknown_av, pred_C2av,'^')
%         hold on
%         h=errorbar(Cknown_av, pred_C2av, pred_C2std);
%         set(h,'linestyle','none');
%         b=pinv([ones(size(pred_C2av,1),1) pred_C2av])*Cknown_av;
%         set(gca,'XLim',[min(Cknown_av)-mean(Cknown_av)*.05 max(Cknown_av)+mean(Cknown_av)*.05]);
%         fn = @(x)[ones(size(x,1),1) x]*b;
%         fplot(fn,get(gca,'XLim'));
%         h=xlabel('Known concentration'); set(h,'fontsize',14);
%         h=ylabel('Predicted concentration'); set(h,'fontsize',14);
%         h=title('Averaged known {\itvs} predicted (NULL data)'); set(h,'fontsize',14);
%     else
%         for i=1:c
%             subplot(c,1,i);
%             plot(Cknown_av(:,i), pred_C2av(:,i),'^')
%             hold on
%             h=errorbar(Cknown_av(:,i), pred_C2av(:,i), pred_C2std(:,i));
%             set(h,'linestyle','none');   
%             b=pinv([ones(size(pred_C2av(:,i),1),1) pred_C2av(:,i)])*Cknown_av(:,i);
%             fn = @(x)[ones(size(x,1),1) x]*b;
%             fplot(fn,[min(Cknown_av(:,i)) max(Cknown_av(:,i))]);
%             h=xlabel('Known concentration'); set(h,'fontsize',14);
%             h=ylabel('Predicted concentration'); set(h,'fontsize',14);
%             h=title(['Averaged known {\itvs} predicted, component ' num2str(i) ' (NULL data)']); 
%             set(h,'fontsize',14);
%         end
%     end
    if c==1
        figure
        hist(output.Q2p,20)
        hold on
        hist(output.Q2p_perm,20)
        h=findobj(gca,'Type','patch');
        set(h(1),'EdgeColor','r');
        set(h(2),'EdgeColor','b');
        hh1=hatchfill(h(1),'single', -45, 3);
        hh2=hatchfill(h(2),'single', 45, 3);
        set(hh1,'Color','r');
        set(hh2,'Color','b');
        legend('Observed distribution', 'Null distribution')
        h=xlabel('Q2p');
        set(hh2,'Color','b')
        set(h,'FontSize',14)
        h=ylabel('No. of hits');
        set(h,'FontSize',14)
        axis('tight')
        set(gca,'YLim',[0 max(get(gca,'YLim'))*1.05])
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


function H = hatchfill(A,STYL,ANGLE,SPACING,FACECOL)
% HATCHFILL Hatching and speckling of patch objects
%   HATCHFILL(A) fills the patch(es) with handle(s) A.
%   A can be a vector of handles or a single handle.
%   If A is a vector, then all objects of A should
%   be part of the same group for predictable results.
%   The hatch consists of black lines angled at
%   45 degrees spaced 5 pixels apart, with no color
%   filling between the lines.
%
%   HATCHFILL(A,STYL) applies STYL pattern with default paramters.
%      - STYL can be 'single' for single lines (the default),
%      'cross' for a double-crossed hatch, 'speckle' for
%     speckling inside the patch boundary, and 'outspeckle' for
%      for speckling outside the boundary. 'fill' will
%      apply only a gray fill and no hatching.
%
%   HATCHFILL(A,STYL,ANGLE,SPACING) applies a hatch/speckle with
%   customized parameters:
%      - ANGLE sets the angle of hatch lines. For speckling, it 
%      controls the width of the speckling region.
%      - SPACING controls the spacing of hatch lines or the
%      density of speckle points.
%      If STYL is 'fill', then ANGLE and SPACING are ignored.
%
%   HATCHFILL(A,STYL,ANGLE,SPACING,FACECOL) allows the user
%   to specify a fill color. (The default is 'none'.)
%
%   H = HATCHFILL(...) returns handles to the line objects
%   comprising the hatch/speckle.
%
%   Examples:
%       Gray region with hatching:
%       hh = hatchfill(a,'cross',45,5,[0.5 0.5 0.5]);
%
%       Speckled region:
%       hatchfill(a,'speckle',7,1);
%
%   NOTE: This function depends on the script hatch_xy.m
%   based on the work of R. Pawlowicz, K. Pankratov, and
%   Iram Weinstein.
%
%   Neil Tandon 11 Jul 2011

    % set defaults:
    if nargin == 1
        STYL = 'single';
        ANGLE = 45;
        SPACING = 5;
        FACECOL = 'none';
    end

    % For backwards compatability:
    if strcmpi(STYL,'none')
        STYL = 'fill';
    end

    if nargin == 2
        if strcmpi(STYL,'single') || strcmpi(STYL,'cross')
            ANGLE = 45;
            SPACING = 5;
            FACECOL = 'none';
        elseif strcmpi(STYL,'speckle') || strcmpi(STYL,'outspeckle')
            ANGLE = 7;
            SPACING = 1;
            FACECOL = 'none';
        elseif strcmpi(STYL,'fill')
            FACECOL = [0.8 0.8 0.8];
        end
    end

    if nargin == 3
        error('Invalid number of input arguments');
    end

    if nargin == 4
        if strcmpi(STYL,'fill')
            FACECOL = [0.8 0.8 0.8];
        else
            FACECOL = 'none';
        end
    end

    if ( ~strcmpi(STYL,'single') && ~strcmpi(STYL,'cross') && ...
         ~strcmpi(STYL,'speckle') && ~strcmpi(STYL,'outspeckle') && ...
         ~strcmpi(STYL,'fill') )
        error(['Invalid style: ',STYL])
    end

    linec = 'k';
    linew = 0.5;
    specksize = 2;

    % axis handle is one or two hierarchical levels up:
    % (Additional check suggested by Dan K)
    hax = get(A(1),'parent');
    is_axes = strcmpi(get(hax,'type'),'axes');
    if ~is_axes
       hax = get(hax,'parent');
    end
    is_axes = strcmpi(get(hax,'type'),'axes');

    x_is_log = 0; y_is_log = 0;
    x_is_reverse = 0; y_is_reverse = 0;

    if is_axes
       axsize_in = get(hax,'position');
       y_is_log = strcmpi(get(hax,'yscale'),'log');
       if y_is_log
           ylims = get(hax,'ylim');
           dy = (ylims(2) - ylims(1))/(log10(ylims(2))-log10(ylims(1)));
           set(hax,'units','pixels');
           axsize = get(hax,'position');
           set(hax,'position',[ axsize(1:3) dy*axsize(4) ]);
           set(hax,'units','normalized')
       end

       x_is_log = strcmpi(get(hax,'xscale'),'log');
       if x_is_log
           xlims = get(hax,'xlim');
           dx = (xlims(2) - xlims(1))/(log10(xlims(2))-log10(xlims(1)));
           set(hax,'units','pixels');
           axsize = get(hax,'position');
           set(hax,'position',[ axsize(1:2) dx*axsize(3) axsize(4) ]);
           set(hax,'units','normalized')
       end

       if strcmp(STYL,'single') || strcmp(STYL,'cross')
          y_is_reverse = strcmpi(get(hax,'ydir'),'reverse');
          if y_is_reverse
              ANGLE = -ANGLE;
          end
          x_is_reverse = strcmpi(get(hax,'xdir'),'reverse');
          if x_is_reverse
              ANGLE = 180-ANGLE;
          end
       end
    end

    % Apply hatch:
    j = 1;
    for k = 1:length(A)
        set(A,'facecolor',FACECOL);
        v = get(A(k),'vertices');
        if any(v(end,:)~=v(1,:))
            v(end+1,:) = v(1,:);
        end
        x = v(:,1);
        if x_is_log
            x = log10(v(:,1));
        end
        y = v(:,2);
        if y_is_log
            y = log10(v(:,2));
        end

        if strcmp(STYL,'fill')
            H = NaN;
            continue
        end

        [xhatch,yhatch] = hatch_xy(x,y,STYL,ANGLE,SPACING);
        if x_is_log
            xhatch = 10.^xhatch;
        end
        if y_is_log
            yhatch = 10.^yhatch;
        end
        if strcmp(STYL,'speckle') || strcmp(STYL,'outspeckle')
            if any(xhatch)
                H(j) = line(xhatch,yhatch,'marker','.','linest','none', ...
                        'markersize',specksize,'color',linec);
                j = j+1;
            end
        elseif strcmp(STYL,'single') || strcmp(STYL,'cross')
            H(j) = line(xhatch,yhatch);
            set(H(j),'color',linec,'linewidth',linew);
            j = j+1;
        end
    end

    if y_is_log || x_is_log
        set(hax,'position',axsize_in);
    end
end

%%%%%%%%%%%%%%%%%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%

function [xi,yi,x,y]=hatch_xy(x,y,varargin);
%
% M_HATCH Draws hatched or speckled interiors to a patch
%       
%    M_HATCH(LON,LAT,STYL,ANGLE,STEP,...line parameters);
%
% INPUTS:
%     X,Y - vectors of points.
%     STYL - style of fill
%     ANGLE,STEP - parameters for style
%
%     E.g.
%                 
%      'single',45,5  - single cross-hatch, 45 degrees,  5 points apart 
%      'cross',40,6   - double cross-hatch at 40 and 90+40, 6 points apart
%      'speckle',7,1  - speckled (inside) boundary of width 7 points, density 1
%                               (density >0, .1 dense 1 OK, 5 sparse)
%      'outspeckle',7,1 - speckled (outside) boundary of width 7 points, density 1
%                               (density >0, .1 dense 1 OK, 5 sparse)
%
%     
%      H=M_HATCH(...) returns handles to hatches/speckles.
%
%      [XI,YI,X,Y]=MHATCH(...) does not draw lines - instead it returns
%      vectors XI,YI of the hatch/speckle info, and X,Y of the original
%      outline modified so the first point==last point (if necessary).
%
%     Note that inside and outside speckling are done quite differently
%     and 'outside' speckling on large coastlines can be very slow.

%
% Hatch Algorithm originally by K. Pankratov, with a bit stolen from 
% Iram Weinsteins 'fancification'. Speckle modifications by R. Pawlowicz.
%
% R Pawlowicz 15/Dec/2005
  
    styl='speckle';
    angle=7;
    step=1/2;

    if length(varargin)>0 & isstr(varargin{1}),
      styl=varargin{1};
      varargin(1)=[];  
    end;
    if length(varargin)>0 & ~isstr(varargin{1}),
      angle=varargin{1};
      varargin(1)=[];  
    end;
    if length(varargin)>0 & ~isstr(varargin{1}),
      step=varargin{1};
      varargin(1)=[];
    end;

    I = zeros(1,length(x));
    %[x,y,I]=m_ll2xy(lon,lat,'clip','patch');


    if x(end)~=x(1) & y(end)~=y(1),
      x=x([1:end 1]);
      y=y([1:end 1]);
      I=I([1:end 1]);
    end;

    if strcmp(styl,'speckle') | strcmp(styl,'outspeckle'),
      angle=angle*(1-I);
    end;

    if size(x,1)~=1,
     x=x(:)';
     angle=angle(:)';
    end;
    if size(y,1)~=1,
     y=y(:)';
    end;


    % Code stolen from Weinstein hatch
    oldu = get(gca,'units');
    set(gca,'units','points');
    sza = get(gca,'pos'); sza = sza(3:4);
    set(gca,'units',oldu)   % Set axes units back

    xlim = get(gca,'xlim');
    ylim = get(gca,'ylim');
    xsc = sza(1)/(xlim(2)-xlim(1)+eps);
    ysc = sza(2)/(ylim(2)-ylim(1)+eps);

    switch lower(styl),
     case 'single',
      [xi,yi]=drawhatch(x,y,angle,step,xsc,ysc,0);
      if nargout<2,
        xi=line(xi,yi,varargin{:});
      end;  
     case 'cross',
      [xi,yi]=drawhatch(x,y,angle,step,xsc,ysc,0);
      [xi2,yi2]=drawhatch(x,y,angle+90,step,xsc,ysc,0);
      xi=[xi,xi2];
      yi=[yi,yi2];
      if nargout<2,
        xi=line(xi,yi,varargin{:});
      end;  
     case 'speckle',
      [xi,yi ]  =drawhatch(x,y,45,   step,xsc,ysc,angle);
      [xi2,yi2 ]=drawhatch(x,y,45+90,step,xsc,ysc,angle);
      xi=[xi,xi2];
      yi=[yi,yi2];
      if nargout<2,
        if any(xi),
          xi=line(xi,yi,'marker','.','linest','none','markersize',2,varargin{:});
        else
          xi=NaN;
        end;    
      end; 
     case 'outspeckle',
      [xi,yi ]  =drawhatch(x,y,45,   step,xsc,ysc,-angle);
      [xi2,yi2 ]=drawhatch(x,y,45+90,step,xsc,ysc,-angle);
      xi=[xi,xi2];
      yi=[yi,yi2];
      inside=logical(inpolygon(xi,yi,x,y)); % logical needed for v6!
      xi(inside)=[];yi(inside)=[];
      if nargout<2,
        if any(xi),
          xi=line(xi,yi,'marker','.','linest','none','markersize',2,varargin{:});
        else
          xi=NaN;
        end;    
      end; 

    end;


    return
end

function [xi,yi]=drawhatch(x,y,angle,step,xsc,ysc,speckle);
%
% This is the guts. 
%

    angle=angle*pi/180;

    % Idea here appears to be to rotate everthing so lines will be
    % horizontal, and scaled so we go in integer steps in 'y' with
    % 'points' being the units in x.
    % Center it for "good behavior".
    ca = cos(angle); sa = sin(angle);
    x0 = mean(x); y0 = mean(y);   
    x = (x-x0)*xsc; y = (y-y0)*ysc;
    yi = x*ca+y*sa;              % Rotation
    y = -x*sa+y*ca;
    x = yi;
    y = y/step;    % Make steps equal to one

    % Compute the coordinates of the hatch line ...............
    yi = ceil(y);
    yd = [diff(yi) 0]; % when diff~=0 we are crossing an integer
    fnd = find(yd);    % indices of crossings
    dm = max(abs(yd)); % max possible #of integers between points


    %
    % This is going to be pretty space-inefficient if the line segments
    % going in have very different lengths. We have one column per line
    % interval and one row per hatch line within that interval.
    %
    A = cumsum( repmat(sign(yd(fnd)),dm,1), 1);

    % Here we interpolate points along all the line segments at the
    % correct intervals.
    fnd1 = find(abs(A)<=abs( repmat(yd(fnd),dm,1) ));
    A  = A+repmat(yi(fnd),dm,1)-(A>0);
    xy = (x(fnd+1)-x(fnd))./(y(fnd+1)-y(fnd));
    xi = repmat(x(fnd),dm,1)+(A-repmat(y(fnd),dm,1) ).*repmat(xy,dm,1);
    yi = A(fnd1);
    xi = xi(fnd1);


     % Sorting points of the hatch line ........................
    %%%yi0 = min(yi); yi1 = max(yi);
    % Sort them in raster order (i.e. by x, then by y)
    % Add '2' to make sure we don't have problems going from a max(xi)
    % to a min(xi) on the next line (yi incremented by one)
    xi0 = min(xi); xi1 = max(xi);
    ci = 2*yi*(xi1-xi0)+xi;
    [ci,num] = sort(ci);
    xi = xi(num); yi = yi(num);


    % if this happens an error has occurred somewhere (we have an odd
    % # of points), and the "fix" is not correct, but for speckling anyway
    % it really doesn't make a difference.
    if rem(length(xi),2)==1, 
      disp('mhatch warning');
      xi = [xi; xi(end)];
      yi = [yi; yi(end)];
    end

     % Organize to pairs and separate by  NaN's ................
    li = length(xi);
    xi = reshape(xi,2,li/2);
    yi = reshape(yi,2,li/2);

    % The speckly part - instead of taking the line we make a point some
    % random distance in.
    if length(speckle)>1 | speckle(1)~=0,

     if length(speckle)>1,
       % Now we get the speckle parameter for each line.

       % First, carry over the speckle parameter for the segment
    %   yd=[0 speckle(1:end-1)];
       yd=[speckle(1:end)];
       A=repmat(yd(fnd),dm,1);
       speckle=A(fnd1);

       % Now give it the same preconditioning as for xi/yi
       speckle=speckle(num);
       if rem(length(speckle),2)==1, 
         speckle = [speckle; speckle(end)];
       end
       speckle=reshape(speckle,2,li/2);

     else
       speckle=[speckle;speckle];
     end;

     % Thin out the points in narrow parts.
     % This keeps everything when abs(dxi)>2*speckle, and then makes
     % it increasingly sparse for smaller intervals.
     oldxi=xi;oldyi=yi;
     dxi=diff(xi);
     nottoosmall=sum(speckle,1)~=0 & rand(1,li/2)<abs(dxi)./(max(sum(speckle,1),eps));
     xi=xi(:,nottoosmall);
     yi=yi(:,nottoosmall);
     dxi=dxi(nottoosmall);
     if size(speckle,2)>1, speckle=speckle(:,nottoosmall); end;
     % Now randomly scatter points (if there any left)
     li=length(dxi);
     if any(li),
       xi(1,:)=xi(1,:)+sign(dxi).*(1-rand(1,li).^0.5).*min(speckle(1,:),abs(dxi) );
       xi(2,:)=xi(2,:)-sign(dxi).*(1-rand(1,li).^0.5).*min(speckle(2,:),abs(dxi) );
       % Remove the 'zero' speckles
       if size(speckle,2)>1,
        xi=xi(speckle~=0);
        yi=yi(speckle~=0);
       end;
      end;

    else
     xi = [xi; ones(1,li/2)*nan];  % Separate the line segments
     yi = [yi; ones(1,li/2)*nan];
    end;
    xi = xi(:)'; yi = yi(:)';

    % Transform back to the original coordinate system
    yi = yi*step;
    xy = xi*ca-yi*sa;
    yi = xi*sa+yi*ca;
    xi = xy/xsc+x0;
    yi = yi/ysc+y0;
end
