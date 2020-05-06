%% Run experiments with the SOUND dataset
seed = 12345
data_sets = {{'SOUND', 0, false}};
cov_funcs = {@(varargin) covSEard(varargin{:})};
approx_algos = {TextBookCG(), FITC(), KMCG(TextBookCG()), DTC(), FOM(), KMCG(FOM()), SoD()}; %#ok<NASGU>
folds = 2; %#ok<NASGU>
fold = 1; %#ok<NASGU>

repetitions = 10;

extended_recording = false; %#ok<NASGU>

experiment;

%% add one more run for DTC and FITC where inducing inputs are spaced over a regular grid (as was done in the original experiments)
seed = seed + repetitions;
approx_algos = {DTC1Dgrid(), FITC1Dgrid()};
folds = 2;
fold = 1;

repetitions = 1;

extended_recording = false;

experiment;