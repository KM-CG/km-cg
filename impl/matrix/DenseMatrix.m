classdef DenseMatrix < handle
    % Matrix wrapper for dense matrices
    properties
        K
        p
        ip % inverse of the active permutation
    end
    
    methods        
        function self = DenseMatrix(K)
            self.K = K;
        end
        
        function f = get_number_of_inducing_inputs_f(self)
            N = size(self.K, 1);
            f = @(cg_step) ceil(sqrt(cg_step * N));
        end
        
        function setPermutation(self, p)
            self.p = p;
            self.ip = zeros(size(p));
            self.ip(p) = 1:length(p); % store the inverse permutation
        end
        
        function Ka = mtimes(self, a)
            Ka = zeros(size(self.K, 1), size(a, 2));
            Ka(self.ip, :) = self.K * a(self.ip, :);
        end
        
        function d = diag(self)
            d = diag(self.K);
        end
        
        function varargout = subsref(self,s)
            switch s(1).type
                case '()'
                    if length(s) == 1
                        s_ = s;
                        s_.subs{1} = self.p(s_.subs{1});
                        s_.subs{2} = self.p(s_.subs{2});
                        [varargout{1:nargout}] = builtin('subsref',self.K, s_);
                    else
                        [varargout{1:nargout}] = builtin('subsref',self.K, s);
                    end
                otherwise
                    [varargout{1:nargout}] = builtin('subsref',self, s);
            end
        end
    end
end

