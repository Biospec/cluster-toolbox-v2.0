function perm_idx = reps_perm(reps)

unique_reps = unique(reps);
no_reps = length(unique_reps);
reps_after_perm = unique_reps(randperm(no_reps));
no_samples = length(reps);
perm_idx = zeros(no_samples, 1);
for i=1:no_reps
    perm_idx(reps==unique_reps(i)) = find(ismember(reps, ...
        reps_after_perm(i)));
end