seed = 12345
data_sets = {{'linspace_distorted_10_10_uniform_test_set', 0, false}, {'linspace_distorted_1000_1000_uniform_test_set', 0, false}, {'linspace_distorted_100_100_100_uniform_test_set', 0, false}, {'linspace_distorted_30_30_30_30_uniform_test_set', 0, false}};

cov_funcs = {@(varargin) covSEard(varargin{:})};
approx_algos = {KMCG(FOM()), DTC(), FITC(), FOM(), KMCG(TextBookCG()), TextBookCG()};

folds = 2;
fold = 1;

repetitions = 10;

extended_recording = false;

experiment;

