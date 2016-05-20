function output = msca(data, subj_id, maxfac)
% output = msca(data, subj_id)
% multilevel simultaneous component analysis to overcome high between
% subjects varibility
% data: the data matrix
% subj_id: index for each individual subjects
% maxfac: the number of principal components to be extracted
% output 
%       .between_matrix, between-subjects variation matrix
%       .within_matrix, within subjects variation matrix
%       .scores/.loadings/.explv, the scores, loadings and explained
%       variance% of each variation matrix.


unique_subj=unique(subj_id);
no_subj=length(unique_subj);
[m,n]=size(data);
data_mc=data-repmat(mean(data),m,1);
between_mat=zeros(no_subj,n);
within_mat=data_mc;
for i=1:no_subj
    idx=find(subj_id==unique_subj(i));
    local_centre=mean(data_mc(idx,:));
    between_mat(i,:)=local_centre;
    within_mat(idx,:)=within_mat(idx,:)-repmat(local_centre,length(idx),1);
end
output.between_mat=between_mat;
output.within_mat=within_mat;
disp('PCA on between subjects matrix')
[tt_between, pp_between, pr_between] = pca(between_mat, maxfac);
disp('PCA on within subjects matrix')
[tt_within, pp_within, pr_within] = pca(within_mat, maxfac);
output.scores_between = tt_between;
output.loadings_between = pp_between;
output.explv_between = pr_between;
output.scores_within = tt_within;
output.loadings_within = pp_within;
output.explv_within = pr_within;

end

