function file_path_name = get_file_name(data_set, folds, fold, cov_func, algo, seed)
% creates a file name from the given parameters
    file_path_name  = sprintf('%s_folds_%i_fold_%i_kernel_%s_algo_%s_seed_%d', data_set, folds, fold, cov_func, algo, seed);
    file_path_name(~ismember(file_path_name, ['A':'Z' 'a':'z' '0':'9'])) = '_';
    file_path_name = ['results' filesep file_path_name];
end