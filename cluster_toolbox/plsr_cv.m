function output = plsr_cv(X, Y, k)
% output = plsr_cv(X, Y, k)
% PLS-R with k-fold double cross-validation
% X:    the data matrix
% Y:    the concentration vector or binary matrix
% k:    is the number of folds, will perform leave-one concentration-out 
%       cross-validation if k is ignored
% output:

no_samples = size(X, 1);
unique_conc = unique(Y,'rows');
no_conc = size(unique_conc,1);
if ~exist('k','var')
    step=1;
else
    step=fix(no_conc/k);
end    
data = X;
no_pcs = 15;
pred_L = zeros(size(Y));
count = 1;
for i=1:step:no_conc
    disp(i);
    tst_reps = unique_conc(i:min(no_conc, i+step-1));
    trn_reps = unique_conc;
    trn_reps(ismember(unique_conc, tst_reps)) = [];
    tst_idx = ismember(Y, tst_reps, 'rows');
    trn_idx = ismember(Y, trn_reps, 'rows');
    data_trn = data(trn_idx,:);
    label_trn = Y(trn_idx,:);
    no_reps_trn = length(trn_reps);
    no_sample_trn = size(data_trn,1);
    data_centre = mean(data_trn); 
    conc_centre = mean(label_trn);
%     data_trn = data_trn-ones(no_sample_trn,1)*data_centre; 
    data_tst=X(tst_idx,:);   
    if no_reps_trn>14*no_conc
        step_trn=round(no_reps_trn/7);
    else
        step_trn=1;
    end
    pred_conc_cv=zeros(size(label_trn));
    for ii=1:step_trn:no_reps_trn
       val_idx=ii:min(ii+step_trn-1,no_reps_trn);
       val_reps=trn_reps(val_idx);
       val_idx=find(ismember(label_trn, val_reps));
       data_val=data_trn(val_idx,:);
       label_val=label_trn(val_idx,:);
       data_trn_inner=data_trn; data_trn_inner(val_idx,:)=[];
       data_centre_inner=mean(data_trn_inner); 
       data_trn_inner=data_trn_inner-repmat(data_centre_inner, ...
           size(data_trn_inner,1),1);
%            data_val=data_val-repmat(data_centre_inner, length(val_idx),1);
       label_trn_inner=label_trn; label_trn_inner(val_idx,:)=[];
       label_centre_inner=mean(label_trn_inner);
       [T,P,Q,W,b] = pls(data_trn_inner,label_trn_inner,...
           min(no_pcs,no_sample_trn-step_trn));
       for iii=1:min(no_pcs,no_sample_trn-step_trn)
           pred_L_val = plspred2(data_val,P,Q,W,b,iii, data_centre_inner,...
               label_centre_inner);
           pred_conc_cv(val_idx,:,iii) = pred_L_val;
       end
    end
    for ii=1:min(no_pcs,no_sample_trn-step_trn)
        rmsecv_inner(ii,:) = sqrt(sum((label_trn - pred_conc_cv(:,:,ii)).^2) ...
            ./size(pred_conc_cv,1));
    end
    if size(rmsecv_inner,2)==1
        opt_pcs=find(rmsecv_inner==min(rmsecv_inner));
    else
        opt_pcs=find(mean(rmsecv_inner,2) == min(mean(rmsecv_inner,2)));
    end
    if length(opt_pcs)>1
        opt_pcs=opt_pcs(1);
    end
    data_trn_cv = data_trn - repmat(data_centre, size(data_trn,1), 1);
    label_trn_cv = label_trn - repmat(conc_centre, size(label_trn,1),1);
    [T,P,Q,W,b]=pls(data_trn_cv,label_trn_cv,opt_pcs);

    [pred_L_cv, B]=plspred2(data_tst,P,Q,W,b,opt_pcs,data_centre,...
        conc_centre);
    pred_L(tst_idx,:) = pred_L_cv;
    rmsecv_in(count,:) = rmsecv_inner;
    opt_pcs_cv(count)=opt_pcs;
    B_cv{count} = B;
    count=count+1;
end

rmsecv_out = sqrt(sum((Y - pred_L).^2)./no_samples);
Q2 = 1 - sum((Y - pred_L).^2)./sum((Y - repmat(mean(Y), no_samples, 1)).^2);
output.rmsecv_inner = rmsecv_in;
output.rmsecv_outer = rmsecv_out;
output.Q2 = Q2;
output.regression_coeff = B_cv;
output.predC = pred_L;
output.opt_pcs = round(mean(opt_pcs_cv));

