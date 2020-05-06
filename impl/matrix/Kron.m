classdef Kron < DenseMatrix
    properties
        missing
        Ks % the kernel matrices per dimension
        dg % the diagonal of the matrix
        id % indices of the elements that are not missing
        N  % size of the training set
        Np % number of training and test points
        k  % kernel
        X  % inputs
    end
    
    methods
        function self = Kron(k, X, Ks, missing)
            self = self@DenseMatrix([]);
            self.k  = k;
            self.X  = X;
            self.Ks = Ks;
            self.dg = 1;
            for d = 1 : numel(Ks)
                dg      = self.dg * diag(Ks{d})';
                self.dg = dg(:);
            end
            clear dg;
            self.N  = size(X, 1);
            self.Np = self.N + numel(missing);
            self.id = setdiff(1:self.Np, missing);
        end
      
        function f = get_number_of_inducing_inputs_f(self)
            f = @(cg_step) cg_step;
        end
        
        function Ka = mtimes(self, a)
            Ka             = zeros(self.N, size(a, 2));
            for m = 1 : size(a, 2)
                Ka(self.ip, m) = MultiKronMVM(self.Ks, self.Np, self.id, a(self.ip, m));
            end
        end
        
        function dg = diag(self)
            dg = self.dg;
        end
        
        function varargout = subsref(self,s)
            % Overriding this method is not a bottleneck
            switch s(1).type
                case '()'
                    if length(s) == 1
                        if(s.subs{1} == ':')
                            [varargout{1}] = self.subref(s.subs{2});
                        else
                            error('This case should not occur!');
                        end
                    else
                        [varargout{1:nargout}] = builtin('subsref',self,s);
                    end
                otherwise
                    [varargout{1:nargout}] = builtin('subsref',self,s);
            end
        end
        
        function sub = subref(self, range)
            sub = self.k(self.X(self.p, :), self.X(self.p(range), :));
        end
        
        function sub = subref2(self, range)
            % this is significantly slower
            M   = length(range);
            sub = zeros(self.N, M);
            e   = zeros(self.N, 1);
            for m = 1 : M
                e(range(m)) = 1;
                sub(:, m)   = self.mtimes(e);
                e(range(m)) = 0;
            end
        end
    end
end

