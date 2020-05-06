seed = 12345
data_sets = {{'MPG', 100, true}, {'AILERONS', 100, true}, {'ABALONE', 100, true}, {'PRECIPITATION', 100, true}, {'ELEVATORS', 100, true}, {'PUMADYN', 100, true}, {'POLETELECOMM', 100, true}, {'TOY', 100, true}};
cov_funcs = {@(varargin) covMaternard(5, varargin{:}), @(varargin) covSEard(varargin{:})};
approx_algos = {TextBookCG(), FITC(), KMCG(TextBookCG()), DTC(), FOM(), KMCG(FOM())};

folds = 2;
fold = 1;

repetitions = 10;

extended_recording = true;

experiment;