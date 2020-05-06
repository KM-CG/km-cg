function [ str ] = covfunc2str( cf )
%COVFUNC2STR Converts the covariance functions used in the experiment into
%a meaningful string used for the results file names.
    str    = func2str(cf);
    prefix = '@(varargin)';
    if ~strcmp(str(1:numel(prefix)), prefix), error('Unclear how to convert the function "%s" into a filename.', str); end
    suffix1  = ',varargin{:})';
    suffix2  = '(varargin{:})';
    if ~strcmp(str(end-numel(suffix1)+1:end), suffix1) ...
            && ~strcmp(str(end-numel(suffix2)+1:end), suffix2)
        error('Unclear how to convert the function "%s" into a filename.', str); 
    end
    % Since the suffixes are of same length we can just use length of nr. 1
    % for removal.
    str = str(numel(prefix)+1:end-numel(suffix1));
    str(~ismember(str, ['A':'Z' 'a':'z' '0':'9'])) = '_';
end

