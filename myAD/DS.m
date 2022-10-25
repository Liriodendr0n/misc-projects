classdef DS
    %DS Doubly Recursive Multivariate Automatic Differentiation
    %   10.1080/0025570X.2002.11953128
    
    properties
        f
        df
        n
        m
    end
    
    methods
        %% Constructor
        
        function obj = DS(varargin)
            if nargin == 1          % create variable
                obj.f = varargin{1};
                obj.df = 0;
                obj.m = 0;
                obj.n = 0;
%             elseif nargin == 2
%                 obj.f = varargin{1};
%                 obj.df = varargin{2};
%                 
            elseif nargin == 4      % full access
                obj.f = varargin{1};
                obj.df = varargin{2};
                obj.n = varargin{3};
                obj.m = varargin{4};
            end
        end
        
        %% DS Operators
        
        % traversal operations
        function out = opV(obj)
            %V operator returns "value"
            out = obj.f;
        end
        function out = opD(obj)
            %D operator returns "derivative"
            out = obj.df;
        end
        function out = opL(obj)
            %L operator pops highest order derivatives off the struct
                if obj.m > 1
                    out = DS(opL(obj.f), opL(obj.df), obj.n, obj.m - 1);
                else
                    out = opL(obj.f);
                end
        end
        
        % constant and variable creation
        function out = makeC(obj, Nvar, Mder, val)
            %makeC creates constant structure
                obj.n = Nvar;
                obj.m = Mder;
                if or(Mder == 0, Nvar == 0)
                    out = val;
                elseif or(Mder > 0, Nvar > 0)
                    out = DS(makeC(obj, Nvar-1, Mder, val), makeC(obj, Nvar, Mder-1, zeros(size(val))), obj.n, obj.m);
                end
        end
        function out = makeC1(obj, Nvar, Mder, val)
            %makeC creates constant structure
                obj.n = Nvar;
                obj.m = Mder;
                if or(Mder == 0, Nvar == 0)
                    out = val;
                elseif or(Mder > 0, Nvar > 0)
                    out = DS(makeC(obj, Nvar-1, Mder, val), makeC(obj, Nvar, Mder-1, ones(size(val))), obj.n, obj.m);
                end
        end
        function out = makeDSvar(obj, Nvar, Mder, varnum, val)
            %makeDSvar creates variable structure
            if isa(val, 'DS')
%                 der = DS(1, opD(val), val.n, val.m);
%                 zer = DS(0, 0, val.n, val.m);
%                 der = opD(val);
                der = DSconst(val.n, val.m, 1);
                zer = DSconst(val.n, val.m, 0);
%                 der = 1;
%                 zer = 0;
            else
                der = ones(size(val));
                zer = zeros(size(val));
            end
            if varnum < Nvar
                out = DS(makeDSvar(obj, Nvar-1, Mder, varnum, val), makeC(obj, Nvar, Mder-1, zer), Nvar, Mder);
            elseif varnum == Nvar
                out = DS(makeC(obj, Nvar-1, Mder, val), makeC(obj, Nvar, Mder-1, der), Nvar, Mder);
%             elseif varnum > Nvar
%                 out = DS(makeC(obj, Nvar-1, Mder, val), makeC(obj, Nvar, Mder-1, 1), Nvar, Mder);
            end
        end
        
        % outputs: values, partials, gradients, hessians
        function out = partial(obj, Did)
            %partial extracts a partial derivative value from a DS
            j = length(Did);
            S = sum(Did);
            M = obj.m;
            while j > 0
                while Did(j) > 0
                    obj = opD(obj);
                    Did(j) = Did(j)-1;
%                     disp("D")
                end
                if and(sum(Did) == 0, S == M)
%                     disp("B")
                    break
                else
                obj = opV(obj);
                j = j-1;
%                 disp("V")
                end
            end
            out = obj;
        end
        function out = DSval(obj)
            %DSval extracts the value of a DS
            Did = zeros(obj.n,1);
            out = partial(obj, Did);
        end 
        function out = DSgrad(obj)
            %DSgrad makes the gradient vector of a scalar valued DS
            if length(obj) ~= 1
                disp("Gradients can only be made for scalar valued DSs")
                return
            end
            Did = zeros(obj.n,1);
            out = zeros(obj.n,length(obj));
%             out = DS;
            for i = 1:obj.n
                Did(i) = 1;
%                 partial(obj, Did);
                out(i,1) = partial(obj, Did);
%                 out = [out; partial(obj, Did)];
                Did(i) = 0;
            end
        end
        function out = DSjac(obj)
            %DSjac makes the jacobian matrix of a vector valued DS
            Did = zeros(obj.n,1);
            out = zeros(length(obj), obj.n);
%             out = DS;
            for i = 1:obj.n
                Did(i) = 1;
                partial(obj, Did);
                out(:,i) = partial(obj, Did).';
%                 out = [out, partial(obj, Did).'];
                Did(i) = 0;
            end
        end
        function out = DShess(obj)
            %DShess makes the hessian matrix of a scalar valued DS
%             if length(obj) ~= 1
%                 disp("Hessians can only be made for scalar valued DSs")
%                 return
%             end
            Did = zeros(obj.n,1);
            out = zeros(obj.n);
%             out = DS;
            for i = 1:obj.n
%                 Rout = DS;
                for j = 1:obj.n
                    Did(i) = 1;
                    Did(j) = Did(j)+1;
                    out(i,j) = partial(obj, Did);
%                     Rout = [Rout, partial(obj, Did)];
                    Did(i) = 0;
                    Did(j) = 0;
                end
%                 out = [out; Rout];
            end
        end
        
        %% Arithmetic Operators
        
        % unary operators
        function out = uplus(a)
            out = DS(opV(a), opD(a), a.n, a.m);
        end
        function out = uminus(a)
            out = DS(uminus(opV(a)), uminus(opD(a)), a.n, a.m);
        end
        
        % addition
        function out = plus(a, b)
            if and(isa(a, 'DS'), isa(b, 'DS'))
                out = DS(opV(a) + opV(b), opD(a) + opD(b), a.n, a.m);
            elseif and(isa(a, 'double'), isa(b, 'DS'))
                out = DS(a + opV(b), opD(b), b.n, b.m);
            elseif and(isa(a, 'DS'), isa(b, 'double'))
                out = DS(opV(a) + b, opD(a), a.n, a.m);
            end
        end
        
        % subtraction in terms of addition and unary minus
        function out = minus(a, b)
            out = plus(a, uminus(b));
        end
        
        % multiplication
        function out = times(a, b)
            if and(isa(a, 'DS'), isa(b, 'DS'))
                out = DS(opV(a) .* opV(b), opD(a) .* opL(b) + opD(b) .* opL(a), a.n, a.m);
            elseif and(isa(a, 'double'), isa(b, 'DS'))
                out = DS(a .* opV(b), a .* opD(b), b.n, b.m);
            elseif and(isa(a, 'DS'), isa(b, 'double'))
                out = DS(opV(a) .* b, opD(a) .* b, a.n, a.m);
            end
        end
        
        % division in terms of times and recip
        function out = recip(a)
            out = DS(1./opV(a), -opD(a)./(opL(a).*opL(a)), a.n, a.m);
        end
        function out = rdivide(a, b)
            out = times(a, recip(b));
        end
        function out = ldivide(a, b)
            out = times(recip(a), b);
        end
        
        % powers
        function out = power(a, n)
            out = DS(opV(a).^n, opD(a).*n.*opL(a).^(n-1), a.n, a.m);
        end
        function out = sqrt(a)
            out = power(a, 0.5);
        end
        
        function out = sign(a)
            out = sign(opV(a));
        end
        function out = abs(a)
            out = DS(abs(opV(a)), sign(opL(a)).*opD(a));
            if isreaopL(opV(a)) == false
                disp("Here be dragons! The compelx modulus is not complex differentiable.")
            end
        end
        
        function out = chain(a, b)
            %chain computes a direct product, used in explicit chain rules
            if and(isa(a, 'DS'), isa(b, 'DS'))
                out = DS(opV(a) .* opV(b), opD(a) .* opD(b), a.n, a.m);
            elseif and(isa(a, 'double'), isa(b, 'DS'))
                out = DS(a .* opV(b), a .* opD(b), b.n, b.m);
            elseif and(isa(a, 'DS'), isa(b, 'double'))
                out = DS(opV(a) .* b, opD(a) .* b, a.n, a.m);
            end
        end
        
        %% Logical Operators
        
        % comparisons
        function out = lt(a, b)
            out = lt(opV(a), opV(b));
        end
        function out = gt(a, b)
            out = gt(opV(a), opV(b));
        end
        function out = le(a, b)
            out = le(opV(a), opV(b));
        end
        function out = ge(a, b)
            out = ge(opV(a), opV(b));
        end
        function out = ne(a, b)
            out = ne(opV(a), opV(b));
        end
        function out = eq(a, b)
            out = eq(opV(a), opV(b));
        end
        
        % conditionals
        function out = and(a, b)
            out = and(opV(a), opV(b));
        end
        function out = or(a, b)
            out = or(opV(a), opV(b));
        end
        function out = xor(a, b)
            out = xor(opV(a), opV(b));
        end
        function out = not(a, b)
            out = not(opV(a), opV(b));
        end
        
        % array
        function out = any(a)
            out = any(opV(a));
        end
        function out = all(a)
            out = all(opV(a));
        end
        function out = find(a)
            out = find(opV(a));
        end
        function out = isscalar(a)
            out = isscalar(opV(a));
        end
        function out = isvector(a)
            out = isvector(opV(a));
        end
        function out = ismatrix(a)
            out = ismatrix(opV(a));
        end
        function out = islogical(a)
            out = islogical(opV(a));
        end
        function out = logical(a)
            out = logical(opV(a));
        end
        
        %% Array Operators
        
%         % size BREAKS CONCATENATION?????????
%         function out = size(obj)
%             %size gets size of underlying value node
%             out = size(opV(obj));
%             disp("hello")
%         end
        function out = length(a)
            %size gets size of underlying value node
            out = length(opV(a));
        end
        
        % concatenation
        function [out] = horzcat(a, varargin)
            %horzcat overloads concatenation to propagate it through a DS
            out = a;
            if nargin == 1
                return
            elseif nargin > 2
                out = horzcat(out, horzcat(varargin{1}, varargin{2:end}));
            else
                b = varargin{:};
                if and(isa(a, 'DS'), isa(b, 'DS'))
                    out = DS(horzcat(opV(a), opV(b)), horzcat(opD(a), opD(b)), a.n, a.m);
                elseif and(isa(a, 'double'), isa(b, 'DS'))
                    a = DS(a, zeros(size(a)), b.n, b.m);
                    out = horzcat(a, b);
                elseif and(isa(a, 'DS'), isa(b, 'double'))
                    b = DS(b, zeros(size(b)), a.n, a.m);
                    out = horzcat(a, b);
                end
            end
        end
        function [out] = vertcat(a, varargin)
            %vertcat overloads concatenation to propagate it through a DS
            out = a;
            if nargin == 1
                return
            elseif nargin > 2
                out = vertcat(out, vertcat(varargin{1}, varargin{2:end}));
            else
                b = varargin{:};
                if and(isa(a, 'DS'), isa(b, 'DS'))
                    out = DS(vertcat(opV(a), opV(b)), vertcat(opD(a), opD(b)), a.n, a.m);
                elseif and(isa(a, 'double'), isa(b, 'DS'))
                    a = DS(a, zeros(size(a)), b.n, b.m);
                    out = vertcat(a, b);
                elseif and(isa(a, 'DS'), isa(b, 'double'))
                    b = DS(b, zeros(size(b)), a.n, a.m);
                    out = vertcat(a, b);
                end
            end
        end
        
        % reference and assignment
        % only supported indexing is A( , ) and V(A( , )) or D(A( , ))
        % anything else would mean a DS is treated differently than a tree of reals in a trenchcoat
        function out = subsref(a, s)
            %subsref overloads parenthesis indexing to reach into a DS
            stype = {s.type};
            ssubs = {s.subs};
            sstruct = substruct(stype{1}, ssubs{1});
            parflag = false;    % has the DS been array-indexed yet
            
            if strcmp(stype{1}, '()')
                out = DS(subsref(opV(a), sstruct), subsref(opD(a), sstruct), a.n, a.m);
                stype = stype(2:end);
                ssubs = ssubs(2:end);
                parflag = true; % it now has
            end
            
            if size(stype, 2) ~= 0  % if indices remain after parenthesis
                if parflag == false
                    out = a;    % since the first if didn't cut a slice
                end
                for i = 1:length(ssubs) % apply V/D dot indices
                    if isa(out, 'double')
                        disp("leaf node")
                    end
                    out = out.(ssubs{i});
                end
            end
        end
        
        function out = subsasgn(a, s, b)
            if any(double([s.type]) == 46)
                disp("Subscripted assignment of DS fields is not supported")
                disp("This is a MATLAB limitation")
                return
            end
            out = DS(subsasgn(opV(a), s, opV(b)), subsasgn(opD(a), s, opD(b)), b.n, b.m);
        end

        %% Matrix Operators
        
        % matrix multiplication
        function out = mtimes(a, b)
            if and(isa(a, 'DS'), isa(b, 'DS'))
                out = DS(opV(a) * opV(b), opD(a) * opL(b) + opL(a) * opD(b), a.n, a.m);
            elseif and(isa(a, 'double'), isa(b, 'DS'))
                out = DS(a * opV(b), a * opD(b), b.n, b.m);
            elseif and(isa(a, 'DS'), isa(b, 'double'))
                out = DS(opV(a) * b, opD(a) * b, a.n, a.m);
            end
        end
        
        function out = inv(a)
            out = DS(inv(opV(a)), -inv(opL(a))*opD(a)*inv(opL(a)), a.n, a.m);
        end
        function out = mrdivide(a, b)
            out = mtimes(a, inv(b));
        end
        function out = mldivide(a, b)
            out = mtimes(inv(a), b);
        end
        
        function out = transpose(a)
            out = DS(transpose(opV(a)), transpose(opD(a)), a.n, a.m);
        end
        function out = ctranspose(a)
            out = DS(ctranspose(opV(a)), ctranspose(opD(a)), a.n, a.m);
        end
        
        function out = trace(a)
            out = DS(trace(opV(a)), trace(opD(a)), a.n, a.m);
        end
        
        function out = det(a)
            out = DS(det(opV(a)), det(opV(a)).*trace(inv(opV(a))*opD(a)), a.n, a.m);
        end
        
%         function out = eig(a)
%             [u, ll, v] = eig(opV(a));
%             for i = 1:length(ll)
%                 l(i) = ll(i,i);
%                 dl(i) = trace((u(:,i)*v(:,i)'/(v(:,i)'*u(:,i)))*opD(a));
%             end
%             out = DS(l, dl, a.n, a.m);
%         end

        %% Elementary Functions
        
        % exponential
        function out = exp(a)
            out = DS(exp(opV(a)), exp(opL(a)).*opD(a), a.n, a.m);
        end
        function out = log(a)
            out = DS(log(opV(a)), recip(opL(a)).*opD(a), a.n, a.m);
        end
        
        % circular trig
        function out = sin(a)
            out = DS(sin(opV(a)), cos(opL(a)).*opD(a), a.n, a.m);
        end
        function out = asin(a)
            out = DS(asin(opV(a)), opD(a)./sqrt(1-opL(a).^2), a.n, a.m);
        end
        
        function out = cos(a)
            out = DS(cos(opV(a)), -sin(opL(a)).*opD(a), a.n, a.m);
        end
        function out = acos(a)
            out = DS(acos(opV(a)), -opD(a)./sqrt(1-opL(a).^2), a.n, a.m);
        end
        
        function out = tan(a)
            out = DS(tan(opV(a)), sec(opL(a)).^2.*opD(a), a.n, a.m);
        end
        function out = atan(a)
            out = DS(atan(opV(a)), opD(a)./(1+opL(a).^2), a.n, a.m);
        end
        
        function out = cot(a)
            out = DS(cot(opV(a)), -csc(opL(a)).^2.*opD(a), a.n, a.m);
        end
        function out = acot(a)
            out = DS(acot(opV(a)), -opD(a)./(1+opL(a).^2), a.n, a.m);
        end
        
        function out = sec(a)
            out = DS(sec(opV(a)), sec(opL(a)).*tan(opL(a)).*opD(a), a.n, a.m);
        end
        function out = asec(a)
            out = DS(asec(opV(a)), opD(a)./(opL(a).^2.*sqrt(1-1./(opL(a).^2))), a.n, a.m);
        end
        
        function out = csc(a)
            out = DS(csc(opV(a)), -cot(opL(a)).*csc(opL(a)).*opD(a), a.n, a.m);
        end
        function out = acsc(a)
            out = DS(acsc(opV(a)), -opD(a)./(opL(a).^2.*sqrt(1-1./(opL(a).^2))), a.n, a.m);
        end
        
        % hyperbolic trig
        function out = sinh(a)
            out = DS(sinh(opV(a)), cosh(opL(a)).*opD(a), a.n, a.m);
        end
        function out = asinh(a)
            out = DS(asinh(opV(a)), opD(a)./sqrt(opL(a).^2+1), a.n, a.m);
        end
        
        function out = cosh(a)
            out = DS(cosh(opV(a)), sinh(opL(a)).*opD(a), a.n, a.m);
        end
        function out = acosh(a)
            out = DS(acosh(opV(a)), opD(a)./sqrt(opL(a).^2-1), a.n, a.m);
        end
        
        function out = tanh(a)
            out = DS(tanh(opV(a)), sech(opL(a)).^2.*opD(a), a.n, a.m);
        end
        function out = atanh(a)
            out = DS(atanh(opV(a)), opD(a)./(1-opL(a).^2), a.n, a.m);
        end
        
        function out = coth(a)
            out = DS(coth(opV(a)), -csch(opL(a)).^2.*opD(a), a.n, a.m);
        end
        function out = acoth(a)
            out = Ds(acoth(opV(a)), opD(a)./(1-opL(a).^2), a.n, a.m);
        end
        
        function out = sech(a)
            out = DS(sech(opV(a)), -sech(opL(a)).*tanh(opL(a)).*opD(a), a.n, a.m);
        end
        function out = asech(a)
            out = DS(asech(opV(a)), -opD(a)./(opL(a).^2.*sqrt(1./(opL(a).^2)-1)), a.n, a.m);
        end
        
        function out = csch(a)
            out = DS(csch(opV(a)), -coth(opL(a)).*csch(opL(a)).*opD(a), a.n, a.m);
        end
        function out = acsch(a)
            out = DS(acsch(opV(a)), -opD(a)./(opL(a).^2.*sqrt(1./(opL(a).^2)+1)), a.n, a.m);
        end
        
        %% Special Functions
        
        
        
        %%
        
    end
end

