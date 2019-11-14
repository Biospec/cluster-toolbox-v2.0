function [perm_idx, reps_permuted] = reps_perm(reps)

unique_reps = unique(reps);
no_reps = length(unique_reps);
reps_after_perm = unique_reps(randperm(no_reps));
no_samples = length(reps);
perm_idx = zeros(no_samples, 1);
idx = 1;
for i=1:no_reps
    reps_found = find(reps == reps_after_perm(i));
    no_reps_found = length(reps_found);
    perm_idx(idx:idx + no_reps_found - 1) = reps_found;
    reps_permuted(idx:idx + no_reps_found - 1) = i;
    idx = idx + no_reps_found;
end