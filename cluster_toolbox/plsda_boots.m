function output = plsda_boots(data, label, rep_idx, no_pcs, no_loops, rebalance)

% output = plsda_boots(data, label, rep_idx, no_pcs, no_loops)
%inputs:
%   data: the X block, each row is one sample
%   label: the class labels of the samples
%   rep_idx: a vector of indices indicating which samples are the analytical 
%             replicates of the same biological sample, e.g. the rep_idx of an 
%             experiment having 5 biological samples with 3 analytical reps each 
%            should be [1;1;1;2;2;2;3;3;3;4;4;4;5;5;5]. The bootstrapping 
%             will be performed on biological reps, i.e. all the analytical reps
%             of the same biological sample are to be selected as a whole. 
%             Assuming 1 rep for each biological sample if rep_idx is ignored.
%   no_pcs: the number of maximum PLS components to be extracted, the optimal
%             number of PLS component will be determined automatically using a 
%             inner cross-validation, set to 15 if ignored.
%   no_loops: the number of iterations of bootstrapping resampling, set to
%            1000 if ignored.
%   rebalance: set to 1 to enable SMOTE algorithm to compensate imbalanced 
%            class size, set to 0 or ignoring this to disable SMOTE algorithm.
%
%outputs:
%       .Fcv: F-scores of inner cross-validation
%       .Fscores: F-scores of test set of each bootstrapping iteration
%       .Fscores2: F-scores of permuted test set of each bootstrapping
%          iteration
%       .precision: Precision scores of test sets
%       .sensitivity: sensitivity scores of test sets
%       .specificity: specificity scores of test sets
%       .fnr: false negative rate of test sets
%       .PLS_loadings: PLS loadings
%       .PLS_scores: PLS scores
%       .ccr: correct classification rates (CCR) 
%       .ccr_avg: 
%       .ccr_cv: CCR of inner cross-validation
%       .ccr_perm: correct classification rates of the null models
%       .ccr_perm_avg: averaged CCR of permuted test sets of all 
%           bootstrapping iterations
%       .class_known: known class labels of each bootstrapping
%          iteration
%       .class_pred: predicted class labels of each bootstrapping
%          iteration
%       .conf_mat: confusion matrices stored in a 3-D matrix
%       .conf_mat_perm: confusion matrices of the null models
%       .conf_mat_avg: averaged confusion matrix
%       .conf_mat_perm_avg: averaged confusion matrix of the null
%          models
%       .label_fittness: the correlation between the permuted labels 
%          and known labels
%       .label_pred: raw predicted class labels
%       .opt_pcs: optimal number of PLS components of each
%          bootstrapping iteration
%       .pval: the p-value
%       .sample_idx: the indices of samples used for testing in each
%          bootstrapping iteration. Use "unique(rep_idx(sample_idx))" to 
%          backtrack the replicates ids are if there were replicates.
%       .vip_scores: vip_scores

% Last update: 12/10/2015, Yun Xu


    if nargin<2
        help plsda_boots;
        output=NaN;
        return;
    end
    [no_samples, no_var]=size(data);
    no_class=size(label,2);
    
    if nargin==2
        no_pcs=15;
        rep_idx=1:no_samples;
        rep_idx=rep_idx(:);
        no_loops=1000;
        rebalance=0;
    end        
    
    
    if nargin == 3
        no_pcs=15;
        no_loops=1000;
        rebalance=0;
    end
    
    if nargin==4
        no_loops=1000;
        rebalance=0;
    end
    
    if nargin == 5
        rebalance = 0;
    end
    
    if isempty(rep_idx) 
        rep_idx = 1:no_samples;
    end
    
    if no_class==1
        unique_class=unique(label);
        label_pls2=zeros(no_samples,length(unique_class));
        for i=1:length(unique_class);
            label_pls2(label==unique_class(i),i)=1;
        end
    else
        label_pls2 = label;
    end
    no_class=size(label_pls2,2);
    label=label_pls2;
       
    rep_ids=unique(rep_idx);
    no_reps=length(rep_ids);

    % shuffle replicates ids to avoid unbalanced class distribution in
    % inner CV
    rep_ids2 = rep_ids(randperm(no_reps));
    for i = 1:no_reps
        idx = find(ismember(rep_idx, rep_ids2(i)));
        rep_idx2(idx, 1) = i;
    end
    rep_idx = rep_idx2;

    
    h= waitbar(0,'Please wait...');
    disp('Starting bootstrapping resampling...')
    i=0;
    while i<no_loops
        clear pred_L pred_L2 opt_pc class_tst_pred class_tst_pred2
        trn_reps=rep_ids(unique(randi(no_reps,round(no_reps*1),1)));
        trn_idx=find(ismember(rep_idx, trn_reps));
        tst_idx=find(~ismember(rep_idx, trn_reps));
        rep_idx_trn=rep_idx(trn_idx);
        data_trn=data(trn_idx,:);
        data_tst=data(tst_idx,:);
        no_reps_trn=length(trn_reps);
        no_sample_trn=size(data_trn,1);

        label_trn=label(trn_idx,:); 
        [tmp,class_trn]=max(label_trn,[],2);
        label_tst=label(tst_idx, :);
        [tmp, class_tst]=max(label_tst,[],2);

        if length(unique(class_trn))==no_class && length(unique(class_tst))==no_class
            i=i+1;
        else
            continue;
        end
        
            
     
        %if there were more than 14 reps, perform 7-fold cross-validation
        if no_reps_trn>14*no_class
            fold_step=round(no_reps_trn/7);
        else
            fold_step=1;
        end
        class_val=zeros(no_sample_trn,min(no_pcs,no_sample_trn-fold_step));
        for ii=1:fold_step:no_reps_trn
           val_idx=ii:min(ii+fold_step-1,no_reps_trn);
           val_reps=trn_reps(val_idx);
           val_idx=find(ismember(rep_idx_trn, val_reps));
           data_val=data_trn(val_idx,:);
           label_val=label_trn(val_idx,:);
           class_trn_inner = class_trn;
           class_trn_inner(val_idx) = [];
           data_trn_inner=data_trn; data_trn_inner(val_idx,:)=[];
           label_trn_inner = label_trn; label_trn_inner(val_idx, :) = [];
           [tmp, class_trn_inner]=max(label_trn_inner,[],2);

           if rebalance
               [data_trn_inner, class_trn_inner] = smote_sel(data_trn_inner, ...
                   class_trn_inner, 5);
               label_trn_inner = zeros(size(data_trn_inner, 1), no_class);
               for k = 1:no_class
                   label_trn_inner(class_trn_inner == k, k) = 1;
               end               
           end
           label_centre_inner=mean(label_trn_inner);
           label_trn_inner = label_trn_inner - repmat(label_centre_inner, ...
               size(label_trn_inner, 1), 1);
           for iii = 1:no_class
               idx = find(class_trn_inner == iii);
               if numel(idx) > 1
                   cls_centres_inner(iii,:) = median(data_trn_inner(idx, :));
               elseif numel(idx) == 1
                   cls_centres_inner(iii,:) = data_trn_inner(idx,:);
               else
                   cls_centres_innter(iii,:) = 0;
               end
           end
           data_centre_inner = median(cls_centres_inner);
           data_trn_inner=data_trn_inner-repmat(data_centre_inner, size(data_trn_inner,1),1);
              
           
           [T,P,Q,W,b]=pls(data_trn_inner,label_trn_inner,min(no_pcs,no_sample_trn-fold_step));
           class_val_inner=zeros(size(data_val,1),min(no_pcs,no_sample_trn-fold_step));

           for iii=1:size(W,2)
               pred_L_val=plspred2(data_val,P,Q,W,b,iii, data_centre_inner,...
                   label_centre_inner,ones(1,no_var), ones(1,no_class));
               [tmp,class_val_inner(:,iii)]=max(pred_L_val,[],2);           
           end
           class_val(val_idx, :)=class_val_inner;
        end
        ccr_cv = zeros(min(no_pcs,no_sample_trn-fold_step), 1);
        for ii=1:min(no_pcs,no_sample_trn-fold_step)
            F(ii,:) = fscores(class_trn, class_val(:,ii));
            F_cv(ii) = nanmean(F(ii,:));        
            correct_idx=find(class_val(:,ii)==class_trn);
            ccr_cv(ii)=(length(correct_idx))/length(class_trn);
        end
        if all(isnan(F_cv))
            opt_pcs = find(ccr_cv == max(ccr_cv));
        else
            opt_pcs=find(F_cv==max(F_cv));
        end
        if length(opt_pcs)>1
            opt_pcs=opt_pcs(1);
        elseif isempty(opt_pcs)           
            warning('Cannot find an optimal LVs, might be an unfortunate combination? Discard current iteration.')
            i = i - 1;
            continue
        end
        
        if rebalance
            no_samples_trn = size(data_trn, 1);
            [data_trn, class_trn] = smote_sel(data_trn, class_trn, 5);
            label_trn = zeros(size(data_trn,1), no_class);
            for k = 1:length(unique(class_trn_inner))
                label_trn(class_trn == k, k) = 1;
            end
            no_samples_added = size(data_trn, 1) - no_samples_trn;
            perm_idx = reps_perm(rep_idx_trn); perm_idx = perm_idx(:);
            perm_idx_added = no_samples_trn+1 : no_samples_trn + no_samples_added;
            perm_idx_added = perm_idx_added(:);
            perm_idx = [perm_idx; perm_idx_added(randperm(no_samples_added))];            
        else
            perm_idx = reps_perm(rep_idx_trn);
        end
        
        for iii = 1:no_class
            idx = find(class_trn == iii);
            if numel(idx) > 1
                cls_centres(iii,:) = median(data_trn(idx, :));
            elseif numel(idx) == 1
                cls_centresr(iii,:) = data_trn(idx,:);
            else
                cls_centres(iii,:) = 0;
            end
        end
        data_centre = median(cls_centres_inner);
        data_trn=data_trn-repmat(data_centre, size(data_trn,1),1);
        label_centre=mean(label_trn);
        label_trn = label_trn - repmat(label_centre, size(label_trn, 1), 1);
        
        [T,P,Q,W,b]=pls(data_trn,label_trn,opt_pcs);
        label_trn_perm=label_trn(perm_idx,:); 
        [tmp,class_trn_perm]=max(label_trn_perm,[],2);
        label_fittness=length(find(class_trn_perm==class_trn))/no_sample_trn;
        [T2,P2,Q2,W2,b2]=pls(data_trn,label_trn_perm,opt_pcs);        
        [pred_L, B]=plspred2(data_tst,P,Q,W,b,opt_pcs,data_centre,...
            label_centre,ones(1,no_var), ones(1,no_class));
        vip_scores = vip(T,P,W,B);
        [pred_L2, B2]=plspred2(data_tst,P2,Q2,W2,b2,opt_pcs,data_centre,...
            label_centre,ones(1,no_var), ones(1,no_class));
%         vip_scores=vip(T,P,W,B);
%         vip_scores2=vip(T2,P2,W2,B2);
        [tmp, class_tst_pred]=max(pred_L,[],2);
        [tmp, class_tst_pred2]=max(pred_L2,[],2);
        ccr=length(find(class_tst_pred==class_tst))/length(class_tst);
        [Ftest, precision, recall, FNR, specificity]  = fscores(class_tst, ...
            class_tst_pred);
        ccr2=length(find(class_tst_pred2==class_tst))/length(class_tst); 
        Ftest2 = fscores(class_tst, class_tst_pred2);
        conf_mat=nan(no_class, no_class);
        conf_mat2=nan(size(conf_mat));
        for ii=1:no_class
            for iii=1:no_class
                if ~any(class_tst==ii)
                    continue;
                end
                conf_mat(ii,iii)=length(find(class_tst == ii & class_tst_pred==iii))...
                    /length(find(class_tst==ii));
            end
        end    
        for ii=1:no_class;
            for iii=1:no_class
                if isempty(find(class_tst==ii))
                    continue;
                end
                conf_mat2(ii,iii)=length(find(class_tst == ii & class_tst_pred2==iii))...
                    /length(find(class_tst==ii));
            end
        end
        
        output.label_pred{i} = pred_L;
        output.sample_idx = tst_idx;        
        output.class_known{i}=class_tst;
        output.class_pred{i}=class_tst_pred;        
        output.opt_pcs(i)=opt_pcs;
        output.conf_mat(:,:,i)=conf_mat;        
        output.conf_mat_perm(:,:,i)=conf_mat2;
        output.ccr(i)=ccr;
        output.Fscores(i, :) = Ftest;
        output.Fscores2(i, :) = Ftest2;
        output.precision(i, :) = precision;
        output.sensitivity(i, :) = recall;
        output.specificity(i, :) = specificity;
        output.fnr(i, :) = FNR;
        output.Fcv(i, :) = F(opt_pcs, :);
        output.ccr_perm(i)=ccr2;
        output.label_fittness(i)=label_fittness;
        output.ccr_cv(i)=max(ccr_cv);
        output.vip_scores(:,:,i) = vip_scores;
        output.reg_coeff(:,:,i) = B;
        waitbar(i/no_loops,h,['No. of iterations = ' num2str(i) '/' num2str(no_loops)]);
    end
    delta=output.ccr-output.ccr_perm;
    output.pval=length(find(delta<0))/no_loops;
    output.conf_mat_avg=nanmean(output.conf_mat,3);
    output.conf_mat_perm_avg=nanmean(output.conf_mat_perm,3);
    output.ccr_avg=mean(output.ccr);
    output.ccr_perm_avg=mean(output.ccr_perm);
    close(h)        
    disp('Done!')
    figure
    hist(output.ccr,20)
    hold on
    hist(output.ccr_perm,20)
    h=findobj(gca,'Type','patch');
    set(h(1),'EdgeColor','r');
    set(h(2),'EdgeColor','b');
    hh1=hatchfill(h(1),'single', -45, 3);
    hh2=hatchfill(h(2),'single', 45, 3);
    set(hh1,'Color','r');
    set(hh2,'Color','b');
    legend('Observed distribution', 'Null distribution')
    h=xlabel('Correct classification rate');
    set(hh2,'Color','b')
    set(h,'FontSize',14)
    h=ylabel('No. of hits');
    set(h,'FontSize',14)
    axis('tight')
    set(gca,'YLim',[0 max(get(gca,'YLim'))*1.05])
    
    
    plot_confmat(output.conf_mat_avg);
    
    disp(['Averaged CCR = ' num2str(output.ccr_avg)])
    disp('Averaged confusion matrix is:')
    disp(output.conf_mat_avg);
    no_pcs_final=max(no_class, round(min(output.opt_pcs)));
    
    conc=label_pls2;
    amean=mean(data); 
    cmean=mean(conc);
    data_mc=data-repmat(amean, size(data,1),1);
    conc_mc=conc-repmat(cmean, size(conc,1),1);
    [T,P,Q,W,b]=pls(data_mc,conc_mc,no_pcs_final);
    [~,B] = plspred2(data,P,Q,W,b,no_pcs_final,amean, cmean);
%     vip_scores = vip(T,P,W,B);
%     output.vip_scores=vip_scores;
    output.PLS_loadings=P;
    output.PLS_scores=T;
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

function [C,BV]=plspred2(X,P,Q,W,b,N,amean, cmean,ascal, cscal);
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

function [F, precision, recall, FNR, specificty] = fscores(class_known, class_pred)
    unique_cls = unique(class_known);
    no_cls = numel(unique_cls);
    F = zeros(no_cls, 1);
    for i = 1:no_cls
        TP = numel(find(class_known == i & class_pred == i));
        FP = numel(find(class_known ~= i & class_pred == i));
        FN = numel(find(class_known == i & class_pred ~= i));
        TN = numel(find(class_known ~= i & class_pred ~= i));
        if TP + FP ==0
            precision(i) = 0;
        else
            precision(i) = TP / (TP+FP);
        end
        if TP + FN == 0
            recall(i) = 0;
        else            
            recall(i) = TP / (TP+FN);
        end
        FNR(i) = FN / (FN + TP);
        specificty(i) = TN / (TN + FP);
        
        if precision + recall == 0
            F(i) = 0;
        else
            F(i) = 2 * (precision(i) * recall(i)) / (precision(i) + recall(i));
        end
    end
    
end

function [data_out, label_out] = smote_sel(data, label, k)

    if nargin < 3
        k = 5;
    end

    unique_cls = unique(label);
    no_cls = length(unique_cls);
    data_syn = [];
    label_syn = [];
    for i = 1:no_cls
        sample_no(i, 1) = length(find(label == unique_cls(i)));
    end
    major_clsid = unique_cls(sample_no == max(sample_no));

    for i = 1:no_cls
        if unique_cls(i) == major_clsid
            continue
        end    
        N = round(max(sample_no)/sample_no(i));
        data_work = data(label == unique_cls(i), :);
        NN_idx = nearestneighbour(data_work', 'NumberOfNeighbours', k);
        S = [];
        for ii = 1:sample_no(i)
            syn_sample = populate(data_work(ii, :), N, ...
                data_work(NN_idx(:, ii), :));
            S = [S; syn_sample];        
        end
        new_sample_no = sample_no(i) + size(S,1);
        if new_sample_no > max(sample_no)
            S = S(randperm(size(S, 1)), :);
            S = S(1:(max(sample_no) - sample_no(i)), :);
        end
        data_syn = [data_syn; S];
        label_syn = [label_syn; repmat(unique_cls(i), size(S,1), 1)];
    end

    label_out = [label; label_syn];
    [label_out, sort_idx] = sort(label_out);
    data_out =[data; data_syn];
    data_out = data_out(sort_idx,:);

end

function synthetic = populate(Sample, N, nnarray)
    numattrs = size(nnarray,2);
    synthetic = zeros(N, numattrs);
    k = size(nnarray, 1);
    dif = nnarray - repmat(Sample, k, 1);
    for i = 1:N
        nn = randi(k,1,1);
        gap = rand(1, numattrs);
        synthetic(i,:) = Sample + dif(nn,:).*gap;             
    end
end

function [idx, tri] = nearestneighbour(varargin)
%NEARESTNEIGHBOUR    find nearest neighbours
%   IDX = NEARESTNEIGHBOUR(X) finds the nearest neighbour by Euclidean
%   distance to each point (column) in X from X. X is a matrix with points
%   as columns. IDX is a vector of indices into X, such that X(:, IDX) are
%   the nearest neighbours to X. e.g. the nearest neighbour to X(:, 2) is
%   X(:, IDX(2))
%
%   IDX = NEARESTNEIGHBOUR(P, X) finds the nearest neighbour by Euclidean
%   distance to each point in P from X. P and X are both matrices with the
%   same number of rows, and points are the columns of the matrices. Output
%   is a vector of indices into X such that X(:, IDX) are the nearest
%   neighbours to P
%
%   IDX = NEARESTNEIGHBOUR(I, X) where I is a logical vector or vector of
%   indices, and X has at least two rows, finds the nearest neighbour in X
%   to each of the points X(:, I).
%   I must be a row vector to distinguish it from a single point.
%   If X has only one row, the first input is treated as a set of 1D points
%   rather than a vector of indices
%
%   IDX = NEARESTNEIGHBOUR(..., Property, Value)
%   Calls NEARESTNEIGHBOUR with the indicated parameters set. Property
%   names can be supplied as just the first letters of the property name if
%   this is unambiguous, e.g. NEARESTNEIGHBOUR(..., 'num', 5) is equivalent
%   to NEARESTNEIGHBOUR(..., 'NumberOfNeighbours', 5). Properties are case
%   insensitive, and are as follows:
%      Property:                         Value:
%      ---------                         ------
%         NumberOfNeighbours             natural number, default 1
%            NEARESTNEIGHBOUR(..., 'NumberOfNeighbours', K) finds the closest
%            K points in ascending order to each point, rather than the
%            closest point. If Radius is specified and there are not
%            sufficient numbers, fewer than K neighbours may be returned
%
%         Radius                         positive, default +inf
%            NEARESTNEIGHBOUR(..., 'Radius', R) finds neighbours within
%            radius R. If NumberOfNeighbours is not set, it will find all
%            neighbours within R, otherwise it will find at most
%            NumberOfNeighbours. The IDX matrix is padded with zeros if not
%            all points have the same number of neighbours returned. Note
%            that specifying a radius means that the Delaunay method will
%            not be used.
%
%         DelaunayMode                   {'on', 'off', |'auto'|}
%            DelaunayMode being set to 'on' means NEARESTNEIGHBOUR uses the
%            a Delaunay triangulation with dsearchn to find the points, if
%            possible. Setting it to 'auto' means NEARESTNEIGHBOUR decides
%            whether to use the triangulation, based on efficiency. Note
%            that the Delaunay triangulation will not be used if a radius
%            is specified.
%
%         Triangulation                  Valid triangulation produced by
%                                        delaunay or delaunayn
%            If a triangulation is supplied, NEARESTNEIGHBOUR will attempt
%            to use it (in conjunction with dsearchn) to find the
%            neighbours.
%
%   [IDX, TRI] = NEARESTNEIGHBOUR( ... )
%   If the Delaunay Triangulation is used, TRI is the triangulation of X'.
%   Otherwise, TRI is an empty matrix
%
%   Example:
%
%     % Find the nearest neighbour in X to each column of X
%     x = rand(2, 10);
%     idx = nearestneighbour(x);
%
%     % Find the nearest neighbours to each point in p
%     p = rand(2, 5);
%     x = rand(2, 20);
%     idx = nearestneighbour(p, x)
%
%     % Find the five nearest neighbours to points x(:, [1 6 20]) in x
%     x = rand(4, 1000)
%     idx = nearestneighbour([1 6 20], x, 'NumberOfNeighbours', 5)
%
%     % Find all neighbours within radius of 0.1 of the points in p
%     p = rand(2, 10);
%     x = rand(2, 100);
%     idx = nearestneighbour(p, x, 'r', 0.1)
%
%     % Find at most 10 nearest neighbours to point p from x within a
%     % radius of 0.2
%     p = rand(1, 2);
%     x = rand(2, 30);
%     idx = nearestneighbour(p, x, 'n', 10, 'r', 0.2)
%
%
%   See also DELAUNAYN, DSEARCHN, TSEARCH

%TODO    Allow other metrics than Euclidean distance
%TODO    Implement the Delaunay mode for multiple neighbours

% Copyright 2006 Richard Brown. This code may be freely used and
% distributed, so long as it maintains this copyright line
error(nargchk(1, Inf, nargin, 'struct'));

% Default parameters
userParams.NumberOfNeighbours = []    ; % Finds one
userParams.DelaunayMode       = 'auto'; % {'on', 'off', |'auto'|}
userParams.Triangulation      = []    ;
userParams.Radius             = inf   ;

% Parse inputs
[P, X, fIndexed, userParams] = parseinputs(userParams, varargin{:});

% Special case uses Delaunay triangulation for speed.

% Determine whether to use Delaunay - set fDelaunay true or false
nX  = size(X, 2);
nP  = size(P, 2);
dim = size(X, 1);

switch lower(userParams.DelaunayMode)
    case 'on'
        %TODO Delaunay can't currently be used for finding more than one
        %neighbour
        fDelaunay = userParams.NumberOfNeighbours == 1 && ...
            size(X, 2) > size(X, 1)                    && ...
            ~fIndexed                                  && ...
            userParams.Radius == inf;
    case 'off'
        fDelaunay = false;
    case 'auto'
        fDelaunay = userParams.NumberOfNeighbours == 1 && ...
            ~fIndexed                                  && ...
            size(X, 2) > size(X, 1)                    && ...
            userParams.Radius == inf                   && ...
            ( ~isempty(userParams.Triangulation) || delaunaytest(nX, nP, dim) );
end

% Try doing Delaunay, if fDelaunay.
fDone = false;
if fDelaunay
    tri = userParams.Triangulation;
    if isempty(tri)
        try
            tri   = delaunayn(X');
        catch
            msgId = 'NearestNeighbour:DelaunayFail';
            msg = ['Unable to compute delaunay triangulation, not using it. ',...
                'Set the DelaunayMode parameter to ''off'''];
            warning(msgId, msg);
        end
    end
    if ~isempty(tri)
        try
            idx = dsearchn(X', tri, P')';
            fDone = true;
        catch
            warning('NearestNeighbour:DSearchFail', ...
                'dsearchn failed on triangulation, not using Delaunay');
        end
    end
else % if fDelaunay
    tri = [];
end

% If it didn't use Delaunay triangulation, find the neighbours directly by
% finding minimum distances
if ~fDone
    idx = zeros(userParams.NumberOfNeighbours, size(P, 2));

    % Loop through the set of points P, finding the neighbours
    Y = zeros(size(X));
    for iPoint = 1:size(P, 2)
        x = P(:, iPoint);

        % This is the faster than using repmat based techniques such as
        % Y = X - repmat(x, 1, size(X, 2))
        for i = 1:size(Y, 1)
            Y(i, :) = X(i, :) - x(i);
        end

        % Find the closest points, and remove matches beneath a radius
        dSq = sum(abs(Y).^2, 1);
        iRad = find(dSq < userParams.Radius^2);
        if ~fIndexed
            iSorted = iRad(minn(dSq(iRad), userParams.NumberOfNeighbours));
        else
            iSorted = iRad(minn(dSq(iRad), userParams.NumberOfNeighbours + 1));
            iSorted = iSorted(2:end);
        end

        % Remove any bad ones
        idx(1:length(iSorted), iPoint) = iSorted';
    end
    %while ~isempty(idx) && isequal(idx(end, :), zeros(1, size(idx, 2)))
    %    idx(end, :) = [];
    %end
    idx( all(idx == 0, 2), :) = [];
end % if ~fDone
if isvector(idx)
    idx = idx(:)';
end
end % nearestneighbour

function tf = delaunaytest(nx, np, dim)
switch dim
    case 2
        tf = np > min(1.5 * nx, 400);
    case 3
        tf = np > min(4 * nx  , 1200);
    case 4
        tf = np > min(40 * nx , 5000);

        % if the dimension is higher than 4, it is almost invariably better not
        % to try to use the Delaunay triangulation
    otherwise
        tf = false;
end % switch
end % delaunaytest

function I = minn(x, n)

% Make sure n is no larger than length(x)
n = min(n, length(x));

% Sort the first n
[xsn, I] = sort(x(1:n));

% Go through the rest of the entries, and insert them into the sorted block
% if they are negative enough
for i = (n+1):length(x)
    j = n;
    while j > 0 && x(i) < xsn(j)
        j = j - 1;
    end

    if j < n
        % x(i) should go into the (j+1) position
        xsn = [xsn(1:j), x(i), xsn((j+1):(n-1))];
        I   = [I(1:j), i, I((j+1):(n-1))];
    end
end

end %minn

function [P, X, fIndexed, userParams] = parseinputs(userParams, varargin)
if length(varargin) == 1 || ~isnumeric(varargin{2})
    P           = varargin{1};
    X           = varargin{1};
    fIndexed    = true;
    varargin(1) = [];
else
    P             = varargin{1};
    X             = varargin{2};
    varargin(1:2) = [];

    % Check the dimensions of X and P
    if size(X, 1) ~= 1
        % Check to see whether P is in fact a vector of indices
        if size(P, 1) == 1
            try
                P = X(:, P);
            catch
                error('NearestNeighbour:InvalidIndexVector', ...
                    'Unable to index matrix using index vector');
            end
            fIndexed = true;
        else
            fIndexed = false;
        end % if size(P, 1) == 1
    else % if size(X, 1) ~= 1
        fIndexed = false;
    end

    if ~fIndexed && size(P, 1) ~= size(X, 1)
        error('NearestNeighbour:DimensionMismatch', ...
            'No. of rows of input arrays doesn''t match');
    end
end
% Parse the Property/Value pairs
if rem(length(varargin), 2) ~= 0
    error('NearestNeighbour:propertyValueNotPair', ...
        'Additional arguments must take the form of Property/Value pairs');
end

propertyNames = {'numberofneighbours', 'delaunaymode', 'triangulation', ...
    'radius'};
while length(varargin) ~= 0
    property = varargin{1};
    value    = varargin{2};

    % If the property has been supplied in a shortened form, lengthen it
    iProperty = find(strncmpi(property, propertyNames, length(property)));
    if isempty(iProperty)
        error('NearestNeighbour:InvalidProperty', 'Invalid Property');
    elseif length(iProperty) > 1
        error('NearestNeighbour:AmbiguousProperty', ...
            'Supplied shortened property name is ambiguous');
    end
    property = propertyNames{iProperty};

    switch property
        case 'numberofneighbours'
            if rem(value, 1) ~= 0 || ...
                    value > length(X) - double(fIndexed) || ...
                    value < 1
                error('NearestNeighbour:InvalidNumberOfNeighbours', ...
                    'Number of Neighbours must be an integer, and smaller than the no. of points in X');
            end
            userParams.NumberOfNeighbours = value;

        case 'delaunaymode'
            fOn = strcmpi(value, 'on');
            if strcmpi(value, 'off')
                userParams.DelaunayMode = 'off';
            elseif fOn || strcmpi(value, 'auto')
                if userParams.NumberOfNeighbours ~= 1
                    if fOn
                        warning('NearestNeighbour:TooMuchForDelaunay', ...
                            'Delaunay Triangulation method works only for one neighbour');
                    end
                    userParams.DelaunayMode = 'off';
                elseif size(X, 2) < size(X, 1) + 1
                    if fOn
                        warning('NearestNeighbour:TooFewDelaunayPoints', ...
                            'Insufficient points to compute Delaunay triangulation');
                    end
                    userParams.DelaunayMode = 'off';

                elseif size(X, 1) == 1
                    if fOn
                        warning('NearestNeighbour:DelaunayDimensionOne', ...
                            'Cannot compute Delaunay triangulation for 1D input');
                    end
                    userParams.DelaunayMode = 'off';
                else
                    userParams.DelaunayMode = value;
                end
            else
                warning('NearestNeighbour:InvalidOption', ...
                    'Invalid Option');
            end % if strcmpi(value, 'off')

        case 'radius'
            if isscalar(value) && isnumeric(value) && isreal(value) && value > 0
                userParams.Radius = value;
                if isempty(userParams.NumberOfNeighbours)
                    userParams.NumberOfNeighbours = size(X, 2) - double(fIndexed);
                end
            else
                error('NearestNeighbour:InvalidRadius', ...
                    'Radius must be a positive real number');
            end
    

        case 'triangulation'
            if isnumeric(value) && size(value, 2) == size(X, 1) + 1 && ...
                    all(ismember(1:size(X, 2), value))
                userParams.Triangulation = value;
            else
                error('NearestNeighbour:InvalidTriangulation', ...
                    'Triangulation not a valid Delaunay Triangulation');
            end
    end % switch property

    varargin(1:2) = [];
end % while
if isempty(userParams.NumberOfNeighbours)
    userParams.NumberOfNeighbours = 1;
end
end %parseinputs

