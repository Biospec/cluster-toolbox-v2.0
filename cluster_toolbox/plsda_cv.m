function output = plsda_cv(X, Y, k, rep_idx)
% output = plsda_cv(X, Y, k, rep_idx)
% PLS-DA with k-fold cross-validation
% X:    the data matrix
% Y:    the class label vector or binary matrix
% k:    is the number of folds, will perform leave-one-out cross-validation if 
%       k=number of samples
% rep_idx: a vector of indices indicating which samples are the analytical 
%          replicates of the same biological sample, e.g. the rep_idx of an 
%          experiment having 5 biological samples with 3 analytical reps each 
%          should be [1;1;1;2;2;2;3;3;3;4;4;4;5;5;5]. The bootstrapping 
%          will be performed on biological reps, i.e. all the analytical reps
%          of the same biological sample are to be selected as a whole. 
%          Assuming 1 rep for each biological sample if rep_idx is ignored.
% output:
%       .ccr: correct classification rate
%       .ccr_perm: a "NULL" CCR obtained from permuted labels
%       .class_pred: Predicted labels
%       .conf_mat: confusion matrix
%       .class_perm: permuted labels
%       .class_cv_perm: predicted labels using permuted labels as target
%       .conf_mat_perm: the "NULL" confusion matrix from the permuted labels

[no_samples, no_vars] = size(X);
step = round(no_samples/k);
data = X;
no_pcs = 15;
if size(Y,2)==1
    class_known = Y;
    unique_cls = unique(Y);
    no_cls = length(unique_cls);
    Y2 = zeros(no_samples, no_cls);
    for i = 1:no_cls
        Y2(Y==unique_cls(i),i)=1;
    end
else
    Y2 = Y;
    class_known = max(Y2,[],2);
end

if exist('rep_idx','var')
    unique_reps = unique(rep_idx);
    no_reps = length(unique_reps);
else
    rep_idx = 1:no_samples; 
    unique_reps = rep_idx;
    no_reps = no_samples;
end

% Shuffle the data first
perm_idx = reps_perm(rep_idx);
data = data(perm_idx,:);
Y2 = Y2(perm_idx,:);
rep_idx = rep_idx(perm_idx);
% Prepare a NULL label with order of samples permuted
perm_idx = reps_perm(rep_idx);
Yperm=Y2(perm_idx,:); 
[tmp,class_perm]=max(Yperm,[],2);
class_cv = zeros(no_samples,1);   class_cv_perm = zeros(no_samples,1);
for i=1:step:no_reps
    tst_reps = unique_reps(i:min(no_reps, i+step-1));
    trn_reps = unique_reps;
    trn_reps(ismember(unique_reps, tst_reps)) = [];
    tst_idx = ismember(rep_idx, tst_reps);
    trn_idx = ismember(rep_idx, trn_reps);
    data_trn = data(trn_idx,:);
    no_reps_trn = length(trn_reps);
    no_sample_trn = size(data_trn,1);
    data_centre = mean(data_trn); 
    data_trn = data_trn-ones(no_sample_trn,1)*data_centre; 
    label_trn = Y2(trn_idx,:); 
    label_trn_perm = Yperm(trn_idx,:);
    [tmp,class_trn]=max(label_trn,[],2);
    data_tst=X(tst_idx,:);
    label_tst=Y2(tst_idx,:);
    label_tst_perm = Yperm(tst_idx,:);
    [tmp, class_tst]=max(label_tst,[],2);
    [tmp,class_trn_perm]=max(label_trn_perm,[],2);    
    [tmp, class_trn]=max(label_trn,[],2);        
    label_centre=mean(label_trn);
    label_trn=label_trn-ones(no_sample_trn,1)*mean(label_trn);
    if no_reps_trn>14*no_cls
        step_trn=round(no_reps_trn/7);
    else
        step_trn=1;
    end
    class_val=zeros(no_sample_trn,min(no_pcs,no_sample_trn-step_trn));
    for ii=1:step_trn:no_reps_trn
       val_idx=ii:min(ii+step_trn-1,no_reps_trn);
       val_reps=trn_reps(val_idx);
       val_idx=find(ismember(trn_reps, val_reps));
       data_val=data_trn(val_idx,:);
       label_val=label_trn(val_idx,:);
       data_trn_inner=data_trn; data_trn_inner(val_idx,:)=[];
       data_centre_inner=mean(data_trn_inner); 
       data_trn_inner=data_trn_inner-repmat(data_centre_inner, size(data_trn_inner,1),1);
%            data_val=data_val-repmat(data_centre_inner, length(val_idx),1);
       label_trn_inner=label_trn; label_trn_inner(val_idx,:)=[];
       [tmp, class_trn_inner]=max(label_trn_inner,[],2);
       label_centre_inner=mean(label_trn_inner);
       [T,P,Q,W,b]=pls(data_trn_inner,label_trn_inner,min(no_pcs,no_sample_trn-step_trn));
       class_val_inner=zeros(size(data_val,1),min(no_pcs,no_sample_trn-step_trn));
       delta_inner=zeros(size(class_val_inner));
       for iii=1:size(W,2)
           pred_L_val=plspred2(data_val,P,Q,W,b,iii, data_centre_inner,...
               label_centre_inner,ones(1,no_vars), ones(1,no_cls));
           %[~,class_val_inner(:,iii)]=max(pred_L_val,[],2);
           [tmp,class_val_inner(:,iii)]=max(pred_L_val,[],2);           
       end
       class_val(val_idx, :)=class_val_inner;
    end
    for ii=1:min(no_pcs,no_sample_trn-step_trn)
        correct_idx=find(class_val(:,ii)==class_trn);
        ccr_cv(ii)=(length(correct_idx))/length(class_trn);
    end
    opt_pcs=find(ccr_cv==max(ccr_cv));
    if length(opt_pcs)>1
        opt_pcs=opt_pcs(1);
    end
    [T,P,Q,W,b]=pls(data_trn,label_trn,opt_pcs);
    [T2,P2,Q2,W2,b2]=pls(data_trn,label_trn_perm,opt_pcs);
    [pred_L, B]=plspred2(data_tst,P,Q,W,b,opt_pcs,data_centre,...
        label_centre,ones(1,no_vars), ones(1,no_cls));
    [pred_L_perm, B2]=plspred2(data_tst,P2,Q2,W2,b2,opt_pcs,data_centre,...
        label_centre,ones(1,no_vars), ones(1,no_cls));
%         vip_scores=vip(T,P,W,B);
%         vip_scores2=vip(T2,P2,W2,B2);
    [tmp, class_pred] = max(pred_L,[],2);
    [tmp, class_predperm] = max(pred_L_perm,[],2);
    class_cv(tst_idx) = class_pred;
    class_cv_perm(tst_idx) = class_predperm;
end
ccr = length(find(class_cv == class_known))/no_samples;
ccr_perm = length(find(class_cv_perm == class_perm))/no_samples;
conf_mat = zeros(no_cls, no_cls);
conf_mat_perm = zeros(no_cls, no_cls);
for i = 1:no_cls
    for ii = 1:no_cls
        conf_mat(i,ii) = length(find(class_known==i & class_cv==ii)) / ...
            length(find(class_known==i));
        conf_mat_perm(i,ii) = length(find(class_perm==i & class_cv_perm==ii)) / ...
            length(find(class_perm==i));
    end
end
output.ccr = ccr;
output.ccr_perm = ccr_perm;
output.class_pred = class_cv;
output.conf_mat = conf_mat;
output.class_perm = class_perm;
output.class_cv_perm = class_cv_perm;
output.conf_mat_perm = conf_mat_perm;

    