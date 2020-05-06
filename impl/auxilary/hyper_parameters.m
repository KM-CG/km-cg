function hyp = hyper_parameters( data_set, folds, fold, cov, X, y, init_hyp, optimizationSteps)
    try
        hyp = load(file_name(data_set, folds, fold, cov));
        hyp = hyp.hyp;
    catch
        warning('Could not find hyper-parameters file. Optimizing...');
        hyp  = minimize(init_hyp, @gp, -optimizationSteps, @infExact, [], cov, [], X, y);
        save(file_name(data_set, folds, fold, cov), 'hyp');
    end
end

function s = file_name(data_set, folds, fold, cov)
    s = sprintf('hyper_parameters/%s_folds_%i_fold_%i_cov_%s', data_set, folds, fold, covfunc2str(cov));
end