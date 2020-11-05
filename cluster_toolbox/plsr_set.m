function []=plsr_set(data,conc,set)
%PLSR using a defined training set
%data: pre-processed data
%conc: known concentration for each sample
%set: column telling the code if a spectrum is in the training or test
%     set - use 1 for training and 2 for test


data_tst = data(set==2,:);
data_trn = data(set==1,:);
y_all = conc;
y_tst = conc(set==2,:);
y_trn = conc(set==1,:);

a = 20;
% cross-validation loop on unique concentrations
unique_y_trn = unique(y_trn, 'rows');
for k=1:a % k-fold validation routine
no_pcs=k;
y_pred=nan(size(y_trn));
% excluding the lowest & highest samples in concentrations
for i=2:size(unique_y_trn,1)-1
val_idx=find(ismember(y_trn, unique_y_trn(i,:),'rows'));
data_val = data_trn(val_idx,:);
y_val = y_trn(val_idx,:);
data_trn_inner = data_trn; data_trn_inner(val_idx,:)=[];
y_trn_inner = y_trn; y_trn_inner(val_idx,:)=[];
data_centre=mean(data_trn_inner);
data_mc=data_trn_inner-ones(size(data_trn_inner,1),1)*data_centre;
y_centre=mean(y_trn_inner);
y_mc=y_trn_inner-repmat(y_centre, size(y_trn_inner,1),1);
[T,P,Q,W,b]=pls(data_mc,y_mc,no_pcs);
y_pred_tmp=plspred2(data_val,P,Q,W,b,no_pcs,data_centre, y_centre);
y_pred(find(ismember(y_trn,y_val,'rows')),:)=y_pred_tmp;
end
y_cv_known = y_trn;
y_cv_known(isnan(y_pred(:,1)),:)=[];
y_pred(isnan(y_pred(:,1)),:)=[];
y_pred_cv{k}=y_pred;
% Y_cv(:,a)=Y_pred;
rmsecv(k,:)=sqrt(sum((y_pred-y_cv_known).^2)./size(y_cv_known,1));
Q2(k,:)=1-sum((y_pred-y_cv_known).^2)./...
sum((y_cv_known-repmat(mean(y_cv_known),size(y_cv_known,1),1)).^2);
end

figure
plot(mean(rmsecv,2),'x-')
hold on
xlabel('No. of PLS factors')
ylabel('RMSECV')
box on
hold off

opt_pc=input('optimal number of factors?');

no_pcs=opt_pc;
data_centre=mean(data_trn);
y_centre=mean(y_trn);
data_mc=data_trn-repmat(data_centre, size(data_trn,1),1);
y_mc=y_trn-repmat(y_centre, size(y_trn,1),1);
[T,P,Q,W,b]=pls(data_mc,y_mc,no_pcs);
y_pred_auto=plspred2(data_trn,P,Q,W,b,no_pcs,data_centre,y_centre);
y_pred_test=plspred2(data_tst,P,Q,W,b,no_pcs,data_centre,y_centre);

R2=1-sum((y_pred_auto-y_trn).^2)./...
sum((y_trn-repmat(mean(y_trn), size(y_trn,1),1)).^2);
Q2p=1-sum((y_pred_test-y_tst).^2)./...
sum((y_tst-repmat(mean(y_tst), size(y_tst,1),1)).^2);
Q2cv=Q2(opt_pc,:);
RMSEC=sqrt(sum((y_pred_auto-y_trn).^2)./size(y_trn,1));
RMSEP=sqrt(sum((y_pred_test-y_tst).^2)./size(y_tst,1));
RMSECV=rmsecv(11,:);

results={'factors used:',num2str(no_pcs)
    'R2:',num2str(R2)
    'Q2CV:',num2str(Q2(no_pcs))
    'Q2P:',num2str(Q2p)
    'RMSEC:',num2str(RMSEC)
    'RMSECV:',num2str(rmsecv(no_pcs))
    'RMSEP:',num2str(RMSEP)};


for i=1:7
    str(i,1)=join(results(i,:));
end

dim=[0.15 0.6 0.3 0.3];



figure
hold on
plot(y_trn(:,:), y_pred_auto(:,:),'b^')
plot(y_cv_known(:,i), y_pred_cv{opt_pc}(:,i),'ms')
plot(y_tst(:,i), y_pred_test(:,i),'ro')
line([0 max(y_all(:,i))], [0 max(y_all(:,i))])
xlabel('Known concentration')
ylabel('Predicted concentration')
axis('tight')
legend('Training set','CV set','Test set','Location','SouthEast');
box on
annotation('textbox',dim,'String', str,'FitBoxToText','on')
hold off
end


