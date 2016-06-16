function opt_parameter = dualplsr_tune(X, Y, ker, rep_idx)
% opt_parameter = dualplsda_tune(X, Y, ker, rep_idx)
% Tune kernel parameters and no. of latent variables of K-PLS regression model
% X: the data matrix with each row represents a sample
% Y: the class information vector or binary matrix
% ker: the type of kernel, can be one of the following
%       'linear': linear kernel, will be equivalent to classical PLS-DA
%       'poly': polynomial kernel
%       'rbf': radiant base function kernel
% rep_idx: A vector indicating which samples were the reps of the same
%          sample, i.e. multiple measurements of the same sample.
%          set to 1:no_samples if ignored
% By: Yun Xu, 16/06/2016

[m,n]=size(X); n2=size(Y,2);

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
            cv=zeros(m,n2);            
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
               Ytrn=Y; Ytrn(val_idx,:)=[]; Ycentre = mean(Ytrn);
               Ktr = polyker(Xtrn);
               Ktrte = polyker(Xtrn, Xval);
               [alpha,Yhat] = dualpls(Ktr,Ktrte,Ytrn,k);
               Yhat = Yhat + repmat(Ycentre, size(Yhat,1), 1);
               cv(val_idx,:) = Yhat;
            end
            rmsecv(i,:) = sqrt(sum((cv - Y).^2)/m);
            Q2(i,:) = 1 - sum((cv-Y).^2)./sum((Y - repmat(mean(Y), m, 1)).^2);
        end
        opt_K = find(mean(rmsecv,2)==min(mean(rmsecv,2)));        
        if length(opt_K)>1
            opt_K = opt_K(1);
        end
        opt_parameter.K = opt_K;
        opt_parameter.opt_rmsecv = rmsecv(opt_K);
        opt_parameter.opt_Q2 = Q2(opt_K,:);
        figure
        plot(rmsecv,'o-');
        h=xlabel('No. of LVs'); set(h,'fontsize',14)
        h=ylabel('RMSECV'); set(h, 'fontsize',14);
    case 'poly'
        degree = 2:10;
        cv = zeros(m, n2, length(degree), max_factor);
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
                   Ycentre = mean(Ytrn);
                   Ktr = polyker(Xtrn,[], p);
                   Ktrte = polyker(Xtrn, Xval, p);
                   [alpha,Yhat] = dualpls(Ktr,Ktrte,Ytrn,k);
                   Yhat = Yhat + repmat(Ycentre, size(Yhat, 1),1);
                   cv(val_idx,:,i,ii) = Yhat;
                end
                rmsecv(i,ii,:) = sqrt(sum((cv(:,:,i,ii) - Y).^2)./m);                
                Q2(i,ii,:) = 1 - sum((cv(:,:,i,ii) - Y).^2)./sum((Y - ...
                    repmat(mean(Y),m,1)).^2);
            end
        end
        if n2==1
            rmsecv_av = squeeze(rmsecv);
            Q2_av = squeeze(Q2);
        else
            rmsecv_av = mean(rmsecv,3);
            Q2_av = mean(Q2,3);
        end
        [opt_k, opt_p] = find(rmsecv_av == min(min(rmsecv_av)));
        opt_k = opt_k(1); opt_p = opt_p(1);
        opt_parameter.opt_k = opt_k;
        opt_parameter.opt_p = degree(opt_p);
        opt_parameter.opt_rmsecv = squeeze(rmsecv(opt_k, opt_p,:));
        opt_parameter.cv_mat = rmsecv;
        opt_parameter.opt_Q2 = squeeze(Q2(opt_k, opt_p,:));
        figure
        h = image(rmsecv_av);
        set(h,'CDatamapping','scaled');
        colorbar
        set(gca,'XTickLabel',num2str(degree(:)));
        h = xlabel('Degrees'); set(h,'fontsize',14);
        set(gca,'YTick',1:max_factor);
        h = ylabel('No. of LVs'); set(h,'fontsize',14);
    case 'rbf'
        Xnorm = norm(X)/m/2;
        gamma = logspace(log10(Xnorm/100),log10(Xnorm*100),15);
        cv = zeros(m, n2, length(gamma), max_factor);
        for i = 1:max_factor
            k = i; 
            disp(['Tuning in progress, ' num2str(i) '/' num2str(max_factor) ...
                '...']);
            for ii = 1:length(gamma)
                g = gamma(ii);
                for iii = 1:fold_step:m
                   if iii > no_reps_trn
                       break;
                   end
                   val_idx = iii:min(iii+fold_step-1,no_reps_trn);
                   val_reps = unique_rep(val_idx);
                   val_idx = find(ismember(rep_idx, val_reps));
                   Xval = X(val_idx,:);
                   Yval = Y(val_idx,:);
                   Xtrn = X; 
                   Xtrn(val_idx,:) = [];
                   Xcentre = mean(Xtrn); 
                   Xtrn = Xtrn - repmat(Xcentre, size(Xtrn,1),1);
                   Xval = Xval - repmat(Xcentre, length(val_idx),1);
                   Ytrn=Y; Ytrn(val_idx,:)=[];
                   Ycentre = mean(Ytrn);
                   Ktr = rbf_two(Xtrn, Xtrn, g);
                   Ktrte = rbf_two(Xtrn, Xval, g);
                   [alpha,Yhat] = dualpls(Ktr,Ktrte,Ytrn,k);
                   cv(val_idx,:,i,ii) = Yhat;
                end
                rmsecv(i,ii,:) = sqrt(sum((cv(:,:,i,ii) - Y).^2)./m);
                Q2(i,ii,:) = 1 - sum((cv(:,:,i,ii) - Y).^2)./sum((Y - ...
                    repmat(mean(Y),m,1)).^2);
            end
        end
        if n2==1
            rmsecv_av = squeeze(rmsecv);
            Q2_av = squeeze(Q2);
        else
            rmsecv_av = mean(rmsecv,3);
            Q2_av = mean(Q2,3);
        end
        [opt_k, opt_p] = find(rmsecv_av == min(min(rmsecv_av)));
        opt_k = opt_k(1); opt_p = opt_p(1);
        opt_parameter.opt_k = opt_k;
        opt_parameter.opt_gamma = gamma(opt_p);
        opt_parameter.opt_rmsecv = squeeze(rmsecv(opt_k, opt_p,:));
        opt_parameter.cv_mat = rmsecv;
        opt_parameter.opt_Q2 = squeeze(Q2(opt_k, opt_p,:));
        opt_parameter.Q2_mat = Q2;
        figure
        h = image(rmsecv_av);
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
