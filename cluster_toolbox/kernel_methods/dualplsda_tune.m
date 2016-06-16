function opt_parameter = dualplsda_tune(X, Y, ker, rep_idx)
% opt_parameter = dualplsda_tune(X, Y, ker, rep_idx)
% Tune kernel parameters and no. of latent variables of K-PLSDA model
% X: the data matrix with each row represents a sample
% Y: the class information vector or binary matrix
% ker: the type of kernel, can be one of the following
%       'linear': linear kernel, will be equivalent to classical PLS-DA
%       'poly': polynomial kernel
%       'rbf': radiant base function kernel
% rep_idx: A vector indicating which samples were the reps of the same
%          sample, i.e. multiple measurements of the same sample.
%          set to 1:no_samples if ignored
%
% By: Yun Xu, 16/06/2016

[m,n]=size(X);

if size(Y,2)==1
    class_known = Y;
    unique_cls = unique(Y);
    no_cls = length(unique_cls);
    Y=zeros(m,no_cls);
    for i=1:no_cls
        Y(class_known == unique_cls(i),i)=1;
    end
else
    class_known = max(Y,[],2);
end
if nargin<4
    rep_idx = 1:m;
    rep_idx = rep_idx(:);
end

   
max_factor=min([20 m n]);
unique_rep = unique(rep_idx);
no_reps_trn = length(unique_rep);

if no_reps_trn>=14
    fold_step = round(no_reps_trn/7);
else
    fold_step = 1;
end

switch ker
    case 'linear'
        for i=1:max_factor
            k=i;
            cv=zeros(m,i);            
            for ii=1:fold_step:m
               if ii > no_reps_trn
                   break;
               end
               val_idx=ii:min(ii+fold_step-1,no_reps_trn);
               val_reps=unique_rep(val_idx);
               val_idx=find(ismember(rep_idx, val_reps));
               Xval=X(val_idx,:);
               Yval=Y(val_idx,:);
               Xtrn=X; 
               Xtrn(val_idx,:)=[];
               Xcentre=mean(Xtrn); 
               Xtrn=Xtrn-repmat(Xcentre, size(Xtrn,1),1);
               Xval=Xval-repmat(Xcentre, length(val_idx),1);
               Ytrn=Y; Ytrn(val_idx,:)=[];
               Ktr = polyker(Xtrn);
               Ktrte = polyker(Xtrn, Xval);
               [alpha,Yhat] = dualpls(Ktr,Ktrte,Ytrn,k);
               [pred_raw, pred_val] = max(Yhat,[],2);
               cv(val_idx,i) = pred_val;
            end
            ccr_cv(i,1) = length(find(cv(:,i)==class_known))/m;
        end
        opt_K = find(ccr_cv==max(ccr_cv));
        if length(opt_K)>1
            opt_K = opt_K(1);
        end
        opt_parameter.K = opt_K;
        opt_parameter.max_ccr_cv = ccr_cv(opt_K);       
        figure
        plot(ccr_cv,'o-');
        h=xlabel('No. of LVs'); set(h,'fontsize',14)
        h=ylabel('Cross Validation CCR'); set(h, 'fontsize',14);
    case 'poly'
        degree = 2:10;
        cv = zeros(m, length(degree), max_factor);
        for i = 1:max_factor
            k = i;
            disp(['Tuning in progress, ' num2str(i) '/' num2str(max_factor) ...
                '...']);
            for ii = 1:length(degree)
                p=degree(ii);
                for iii=1:fold_step:m
                   if iii > no_reps_trn
                       break;
                   end
                   val_idx=iii:min(iii+fold_step-1,no_reps_trn);
                   val_reps=unique_rep(val_idx);
                   val_idx=find(ismember(rep_idx, val_reps));
                   Xval=X(val_idx,:);
                   Yval=Y(val_idx,:);
                   Xtrn=X; 
                   Xtrn(val_idx,:)=[];
                   Xcentre=mean(Xtrn); 
                   Xtrn=Xtrn-repmat(Xcentre, size(Xtrn,1),1);
                   Xval=Xval-repmat(Xcentre, length(val_idx),1);
                   Ytrn=Y; Ytrn(val_idx,:)=[];
                   Ktr = polyker(Xtrn,[], p);
                   Ktrte = polyker(Xtrn, Xval, p);
                   [alpha,Yhat] = dualpls(Ktr,Ktrte,Ytrn,k);
                   [pred_raw, pred_val] = max(Yhat,[],2);
                   cv(val_idx,i,ii) = pred_val;
                end
                ccr_cv(i,ii)=length(find(cv(:,i,ii)==class_known))/m;
            end
        end
        [opt_k, opt_p] = find(ccr_cv == max(max(ccr_cv)));
        opt_k = opt_k(1); opt_p = opt_p(1);
        opt_parameter.opt_k = opt_k;
        opt_parameter.opt_p = degree(opt_p);
        opt_parameter.max_ccr_cv = ccr_cv(opt_k, opt_p);
        opt_parameter.cv_mat = ccr_cv;
        figure
        h = image(ccr_cv);
        set(h,'CDatamapping','scaled');
        colorbar
        set(gca,'XTickLabel',num2str(degree(:)));
        h = xlabel('Degrees'); set(h,'fontsize',14);
        set(gca,'YTick',1:max_factor);
        h = ylabel('No. of LVs'); set(h,'fontsize',14);
    case 'rbf'
        Xnorm = sqrt(norm(X)/m);
        gamma = logspace(log10(Xnorm/100),log10(Xnorm*100),15);
        cv = zeros(m, length(gamma), max_factor);
        for i = 1:max_factor
            k = i; 
            disp(['Tuning in progress, ' num2str(i) '/' num2str(max_factor) ...
                '...']);
            for ii = 1:length(gamma)
                g=gamma(ii);
                for iii=1:fold_step:m
                   if iii > no_reps_trn
                       break;
                   end
                   val_idx=iii:min(iii+fold_step-1,no_reps_trn);
                   val_reps=unique_rep(val_idx);
                   val_idx=find(ismember(rep_idx, val_reps));
                   Xval=X(val_idx,:);
                   Yval=Y(val_idx,:);
                   Xtrn=X; 
                   Xtrn(val_idx,:)=[];
                   Xcentre=mean(Xtrn); 
                   Xtrn=Xtrn-repmat(Xcentre, size(Xtrn,1),1);
                   Xval=Xval-repmat(Xcentre, length(val_idx),1);
                   Ytrn=Y; Ytrn(val_idx,:)=[];
                   Ktr = rbf(Xtrn,g);
                   Ktrte = rbf_two(Xtrn, Xval, g);
                   [alpha,Yhat] = dualpls(Ktr,Ktrte,Ytrn,k);
                   [pred_raw, pred_val] = max(Yhat,[],2);
                   cv(val_idx,i,ii) = pred_val;
                end
                ccr_cv(i,ii)=length(find(cv(:,i,ii)==class_known))/m;
            end
        end
        [opt_k, opt_p] = find(ccr_cv == max(max(ccr_cv)));
        opt_k = opt_k(1); opt_p = opt_p(1);
        opt_parameter.opt_k = opt_k;
        opt_parameter.opt_gamma = gamma(opt_p);
        opt_parameter.max_ccr_cv = ccr_cv(opt_k, opt_p);
        opt_parameter.cv_mat = ccr_cv;
        figure
        h = image(ccr_cv);
        set(h,'CDatamapping','scaled');
        colorbar
        gamma=round(log10(gamma)*10)/10;
        set(gca,'XTick',1:15)
        set(gca,'XTickLabel',num2str(gamma(:)));
        h = xlabel('log10({\gamma})'); set(h,'fontsize',14);
        set(gca,'YTick',1:max_factor);
        h = ylabel('No. of LVs'); set(h,'fontsize',14);
    otherwise
        error('unknown kernel, it has to be linear, poly or rbf')
end
