classdef FITC1Dgrid < FITC
    % FITC implementation for the SOUND dataset
    methods(Static)
        function val = getFileName()
            val = which(mfilename);
        end
    end
    
    methods(Access=protected)        
        function [KXU, KUU] = get_inducing_input_matrices(obj, m)
            XU  = linspace(min(obj.X),max(obj.X),m)';
            KXU = obj.k(obj.X, XU);
            KUU = obj.k(XU);
        end
        
        function KXSU = get_test_inducing_cov(obj, m)
            XU   = linspace(min(obj.X),max(obj.X),m)';
            KXSU = obj.k(XU, obj.Xs);
        end
    end
end

