classdef Toeplitz < DenseMatrix
    properties
        r % The vector that defines the Toeplitz matrix
        f % fft of r
        missing % The indices that make this matrix non-Toeplitz.
        notmissing
        id
    end
    
    methods
        function self = Toeplitz(r, missing)
            self = self@DenseMatrix([]);
            self.r = r;
            self.f=fft([r; 0; flipud(r(2:end))]);
            self.notmissing = setdiff(1:size(r, 1), missing);
        end
      
        function f = get_number_of_inducing_inputs_f(self)
            f = @(cg_step) cg_step;
        end
        
        function setPermutation(self, p)
            setPermutation@DenseMatrix(self, p);
            self.id = zeros(size(self.notmissing)); self.id(self.ip) = self.notmissing;
        end
        
        function Ta = mtimes(self, a)
            % https://de.mathworks.com/matlabcentral/fileexchange/8548-toeplitzmult
            % Brian Borchers
            x = zeros(size(self.r, 1), size(a, 2));
            x(self.notmissing, :) = a(self.ip, :);
            v = ifft(bsxfun(@times, self.f, fft([x; zeros(size(x))])));
            Ta = zeros(size(a));
            Ta(self.ip, :)=real(v(self.notmissing, :));
        end
        
        function d = diag(self)
            d = self.r(1) * ones(numel(self.notmissing, 1), 1);
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
            M    = length(range);
            temp = repmat(self.r, [1 M]);
            idx  = abs(repmat(self.id', [1 M])-self.id(range))+1;
            sub  = temp(idx);
        end
        
        function sub = subref2(self, range)
            % This function is slower than the one above.
            M   = length(range);
            N   = length(self.notmissing);
            idx = zeros(N, M); idx(range, range) = eye(M);
            sub = self.mtimes(idx);
        end
    end
end

