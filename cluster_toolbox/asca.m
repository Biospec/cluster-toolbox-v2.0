function output = asca(data, design_mat)

% output = asca(data, design_mat)
% ANOVA-simultaneous component analysis
% data: the m by n data matrix in which each row is a sample and each
%   column is a variable.
% design_mat: a m by c design matrix in which c is the number of factors.
%   different levels of each factor should be coded by different numbers such
%   as 1,2,3 etc.
% by Yun Xu, 2016
% Updated 22/09/2016, changed design_mat from a cell array to a class
%   infomation matrix which is probably more intuitive to use.

design_mat_in = design_mat;
no_factors=size(design_mat_in,2);
no_samples=size(data,1);
clear design_mat;
for i=1:no_factors
    unique_lvls = unique(design_mat_in(:,i));
    no_lvls = length(unique_lvls);
    design_mat_tmp = zeros(no_samples, no_lvls);
    for ii=1:no_lvls
        design_mat_tmp(design_mat_in(:,i) == unique_lvls(ii),ii) = 1;
    end
    design_mat{i} = design_mat_tmp;
end

for i=1:no_factors
    no_lvls(i)=size(design_mat{i},2);
end
D_GM=median(data);
D1=data-repmat(D_GM, no_samples,1);
D_PE=D1;
for i=1:no_factors
    DA_M=zeros(size(data));
    for ii=1:no_lvls(i)
        idx=find(design_mat{i}(:,ii)==1);
        if length(idx)==1
            DA_M(idx,:)=D1(idx,:);
        else
            ka=median(D1(idx,:));
            DA_M(idx,:)=repmat(ka,length(idx),1);
        end
    end
    D1=D1-DA_M;
    D_M{i}=DA_M;
    D_PE=D_PE-DA_M;
end

for i=1:no_factors
    D_test=D_M{i}+D_PE;
%     D_test=D_test-repmat(mean(D_test), size(D_test,1),1);
    [U,S,V]=svd(D_test,'econ');
    scores_test=U*S;
    scores_test=scores_test(:,1:no_lvls(i)-1);
    loadings_test=V(:, 1:no_lvls(i)-1);
%     S=diag(S);
%     ssq_observed(i)=sum(S(:,1:no_lvls(i)-1));
    output.scores_ob{i}=scores_test;
    output.loadings_ob{i}=loadings_test;
    [U,S,V]=svd(D_M{i});
    S=diag(S);
    ssq_observed(i)=sum(S(1:no_lvls(i)-1));
end
output.ssq_ob=ssq_observed;


no_perms=1000; design_mat_perm=design_mat;
disp('Start permutation test...')
h=waitbar(0, 'Please wait...');
for i=1:no_perms
    D1=data-repmat(D_GM, no_samples,1);
    D_PE=D1;
    perm_idx=randperm(no_samples);
    ssq_perm=zeros(1,no_factors);
    design_mat_perm=design_mat;
    for ii=1:no_factors
        design_mat_perm{ii}=design_mat_perm{ii}(perm_idx,:);
    end
    for ii=1:no_factors
        DB_M=zeros(size(data));
        for iii=1:no_lvls(ii)
            idx=find(design_mat_perm{ii}(:,iii)==1);
            if length(idx)==1
                DB_M(idx,:)=D1(idx,:);
            else
                ka=median(D1(idx,:));
                DB_M(idx,:)=repmat(ka,length(idx),1);
            end
        end
        D1=D1-DB_M;
        D_M{ii}=DB_M;
        D_PE=D_PE-DB_M;
    end
    for ii=1:no_factors
%         D_test=D_M{ii}+D_PE;
        [U,S,V]=svd(D_M{ii},'econ');
%         scores_test=U*S;
%         scores_test=scores_test(:,1:no_lvls(i)-1);
%         loadings_test=V(:, 1:no_lvls(i)-1);
        S=diag(S);
        ssq_perm(ii)=sum(S(1:no_lvls(ii)-1));
    end
    output.ssq_perm(i,:)=ssq_perm;
    waitbar(i/no_perms,h)
end

for i=1:no_factors
    output.pval(i,1)=length(find(output.ssq_ob(i)< output.ssq_perm(:,i)))/...
        no_perms;
end
close(h);
disp('Done!')
    