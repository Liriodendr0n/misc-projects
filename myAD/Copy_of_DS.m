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
            elseif nargin == 4      % full access
                if isscalar(varargin{1})
                    obj.f = varargin{1};
                    obj.df = varargin{2};
                    obj.n = varargin{3};
                    obj.m = varargin{4};
                else
                    F = reshape(varargin{1}, [numel(varargin{1}), 1]);
                    DF = reshape(varargin{2}, [numel(varargin{2}), 1]);
                    for i = 1:numel(F)
                        obj(i,1) = DS(F(i), DF(i) , varargin{3}, varargin{4});
                    end
                    obj = reshape(obj, size(varargin{1}));
%                     obj = arrayfun(@(x)DS(x, varargin{2}, varargin{3}, varargin{4}), varargin{1});
                end
            end
                
        end
        
        %% DS operators
        
        function out = opV(obj)
            %V operator returns "value"
            out = reshape([obj.f], size(obj));
        end
        
        function out = opD(obj)
            %D operator returns "derivative"
            out = reshape([obj.df], size(obj));
        end
        
        function out = opL(obj)
            %L operator pops highest order derivatives off the struct
                N = [obj.n];
                M = [obj.m];
                objf = reshape([obj.f], size(obj));
                objdf = reshape([obj.df], size(obj));
                if M(1) > 1
                    out = DS(opL(objf), opL(objdf), N(1), M(1)-1);
                else
                    out = opL(objf);
                end
        end
        
        function out = makeC(obj, Nvar, Mder, val)
            %makeC creates constant structure
                obj.n = Nvar;
                obj.m = Mder;
                if isscalar(val)
                    if or(Mder == 0, Nvar == 0)
                        out = val;
                    elseif or(Mder > 0, Nvar > 0)
                        out = DS(makeC(obj, Nvar-1, Mder, val), makeC(obj, Nvar, Mder-1, 0), obj.n, obj.m);
                    end
                else
                    Vval = reshape(val, [numel(val), 1]);
                    for i = 1:numel(val)
                        out(i,1) = makeC(obj, Nvar, Mder, Vval(i));
                    end
                    out = reshape(out, size(val));
%                     out = arrayfun(@(x)makeC(DS, Nvar, Mder, x), val);
                end
        end
        
        function out = makeDSvar(obj, Nvar, Mder, varnum, val)
            %makeDSvar creates variable structure
            if varnum < Nvar
                out = DS(makeDSvar(obj, Nvar-1, Mder, varnum, val), makeC(obj, Nvar, Mder-1, 0), Nvar, Mder);
            elseif varnum == Nvar
                out = DS(makeC(obj, Nvar-1, Mder, val), makeC(obj, Nvar, Mder-1, 1), Nvar, Mder);
            elseif varnum > Nvar
                out = DS(makeC(obj, Nvar-1, Mder, val), makeC(obj, Nvar, Mder-1, 1), Nvar, Mder);
            end
        end
        
        function out = partial(obj, Did)
            %partial extracts partial derivatives from a DS
            j = length(Did);
            while j > 0
                while Did(j) > 0
                    obj = opD(obj);
                    Did(j) = Did(j)-1;
                end
                obj = opV(obj);
                j = j-1;
            end
            out = obj;
        end
        
        function out = DSval(obj)
            %DSval extracts the value of a DS
            Did = zeros(obj.n,1);
            out = partial(obj, Did);
        end
        
        function out = DSgrad(obj)
            %DSgrad makes the gradient vector of a DS
            Did = zeros(obj.n,1);
            out = zeros(obj.n,1);
            for i = 1:obj.n
                Did(i) = 1;
                out(i,1) = partial(obj, Did);
                Did(i) = 0;
            end
        end
        
        function out = DShess(obj)
            %DShess makes the hessian matrix of a DS
            Did = zeros(obj.n,1);
            out = zeros(obj.n);
            for i = 1:obj.n
                for j = 1:obj.n
                    Did(i) = 1;
                    Did(j) = Did(j)+1;
                    out(i,j) = partial(obj, Did);
                    Did(i) = 0;
                    Did(j) = 0;
                end
            end
        end
        
        %% Array Operators
        
        % unary operators
        function out = uplus(a)
            an = [a.n];
            am = [a.m];
            out = DS(opV(a), opD(a), an(1), am(1));
        end
        function out = uminus(a)
            an = [a.n];
            am = [a.m];
            out = DS(uminus(opV(a)), uminus(opD(a)), an(1), am(1));
        end
        
        % addition
        function out = plus(a, b)
            if and(isa(a, 'DS'), isa(b, 'DS'))
                an = [a.n];
                am = [a.m];
                out = DS(opV(a) + opV(b), opD(a) + opD(b), an(1), am(1));
            elseif and(isa(a, 'double'), isa(b, 'DS'))
%                 out = DS(a + opV(b), opD(b), b.n, b.m);
                bn = [b.n];
                bm = [b.m];
                out = makeC(DS, bn(1), bm(1), a) + b;
            elseif and(isa(a, 'DS'), isa(b, 'double'))
%                 out = DS(opV(a) + b, opD(a), a.n, a.m);
                an = [a.n];
                am = [a.m];
                out = a + makeC(DS, an(1), am(1), b);
            end
        end
        
        % subtraction in terms of addition and unary minus
        function out = minus(a, b)
            out = plus(a, uminus(b));
        end
        
        % multiplication
        function out = times(a, b)
            if and(isa(a, 'DS'), isa(b, 'DS'))
                an = [a.n];
                am = [a.m];
                out = DS(opV(a) .* opV(b), opD(a) .* opL(b) + opD(b) .* opL(a), an(1), am(1));
            elseif and(isa(a, 'double'), isa(b, 'DS'))
%                 out = DS(a .* opV(b), a .* opD(b), b.n, b.m);
                bn = [b.n];
                bm = [b.m];
                out = makeC(DS, bn(1), bm(1), a) .* b;
            elseif and(isa(a, 'DS'), isa(b, 'double'))
%                 out = DS(opV(a) .* b, opD(a) .* b, a.n, a.m);
                an = [a.n];
                am = [a.m];
                out = a .* makeC(DS, an(1), am(1), b);
            end
        end
        
        % division in terms of times and recip
        function out = recip(a)
            an = [a.n];
            am = [a.m];
            out = DS(1./opV(a), -opD(a)./(opL(a).*opL(a)), an(1), am(1));
        end
        function out = rdivide(a, b)
            out = times(a, recip(b));
        end
        function out = ldivide(a, b)
            out = times(recip(a), b);
        end
        
        % powers
        function out = power(a, n)
            an = [a.n];
            am = [a.m];
            out = DS(opV(a).^n, opD(a).*n.*opL(a).^(n-1), an(1), am(1));
        end
        function out = sqrt(a)
            out = power(a, 0.5);
        end
        
        % equality and inequalities
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
        
        function out = sign(a)
            out = sign(opV(a));
        end
        function out = abs(a)
            out = DS(abs(opV(a)), sign(opL(a)).*opD(a));
            if isreaopL(opV(a)) == false
                disp("Here be dragons! The compelx modulus is not complex differentiable.")
            end
        end
        
        function out = sum(a)
            if isvector(a)
                i = 1;
                out = a(1);
                while i <= length(a)
                    out = out + a(i);
                    i = i+1;
                end
            end
        end
        
        %% Elementary Functions
        
        % exponential
        function out = exp(a)
            an = [a.n];
            am = [a.m];
            out = DS(exp(opV(a)), exp(opL(a)).*opD(a), an(1), am(1));
        end
        function out = log(a)
            an = [a.n];
            am = [a.m];
            out = DS(log(opV(a)), recip(opL(a)).*opD(a), an(1), am(1));
        end
        
        % circular trig
        function out = sin(a)
            an = [a.n];
            am = [a.m];
            out = DS(sin(opV(a)), cos(opL(a)).*opD(a), an(1), am(1));
        end
        function out = asin(a)
            an = [a.n];
            am = [a.m];
            out = DS(asin(opV(a)), opD(a)./sqrt(1-opL(a).^2), an(1), am(1));
        end
        
        function out = cos(a)
            an = [a.n];
            am = [a.m];
            out = DS(cos(opV(a)), -sin(opL(a)).*opD(a), an(1), am(1));
        end
        function out = acos(a)
            an = [a.n];
            am = [a.m];
            out = DS(acos(opV(a)), -opD(a)./sqrt(1-opL(a).^2), an(1), am(1));
        end
        
        function out = tan(a)
            an = [a.n];
            am = [a.m];
            out = DS(tan(opV(a)), sec(opL(a)).^2.*opD(a), an(1), am(1));
        end
        function out = atan(a)
            an = [a.n];
            am = [a.m];
            out = DS(atan(opV(a)), opD(a)./(1+opL(a).^2), an(1), am(1));
        end
        
        function out = cot(a)
            an = [a.n];
            am = [a.m];
            out = DS(cot(opV(a)), -csc(opL(a)).^2.*opD(a), an(1), am(1));
        end
        function out = acot(a)
            an = [a.n];
            am = [a.m];
            out = DS(acot(opV(a)), -opD(a)./(1+opL(a).^2), an(1), am(1));
        end
        
        function out = sec(a)
            an = [a.n];
            am = [a.m];
            out = DS(sec(opV(a)), sec(opL(a)).*tan(opL(a)).*opD(a), an(1), am(1));
        end
        function out = asec(a)
            an = [a.n];
            am = [a.m];
            out = DS(asec(opV(a)), opD(a)./(opL(a).^2.*sqrt(1-1./(opL(a).^2))), an(1), am(1));
        end
        
        function out = csc(a)
            an = [a.n];
            am = [a.m];
            out = DS(csc(opV(a)), -cot(opL(a)).*csc(opL(a)).*opD(a), an(1), am(1));
        end
        function out = acsc(a)
            an = [a.n];
            am = [a.m];
            out = DS(acsc(opV(a)), -opD(a)./(opL(a).^2.*sqrt(1-1./(opL(a).^2))), an(1), am(1));
        end
        
        % hyperbolic trig
        function out = sinh(a)
            an = [a.n];
            am = [a.m];
            out = DS(sinh(opV(a)), cosh(opL(a)).*opD(a), an(1), am(1));
        end
        function out = asinh(a)
            an = [a.n];
            am = [a.m];
            out = DS(asinh(opV(a)), opD(a)./sqrt(opL(a).^2+1), an(1), am(1));
        end
        
        function out = cosh(a)
            an = [a.n];
            am = [a.m];
            out = DS(cosh(opV(a)), sinh(opL(a)).*opD(a), an(1), am(1));
        end
        function out = acosh(a)
            an = [a.n];
            am = [a.m];
            out = DS(acosh(opV(a)), opD(a)./sqrt(opL(a).^2-1), an(1), am(1));
        end
        
        function out = tanh(a)
            an = [a.n];
            am = [a.m];
            out = DS(tanh(opV(a)), sech(opL(a)).^2.*opD(a), an(1), am(1));
        end
        function out = atanh(a)
            an = [a.n];
            am = [a.m];
            out = DS(atanh(opV(a)), opD(a)./(1-opL(a).^2), an(1), am(1));
        end
        
        function out = coth(a)
            an = [a.n];
            am = [a.m];
            out = DS(coth(opV(a)), -csch(opL(a)).^2.*opD(a), an(1), am(1));
        end
        function out = acoth(a)
            an = [a.n];
            am = [a.m];
            out = Ds(acoth(opV(a)), opD(a)./(1-opL(a).^2), an(1), am(1));
        end
        
        function out = sech(a)
            an = [a.n];
            am = [a.m];
            out = DS(sech(opV(a)), -sech(opL(a)).*tanh(opL(a)).*opD(a), an(1), am(1));
        end
        function out = asech(a)
            an = [a.n];
            am = [a.m];
            out = DS(asech(opV(a)), -opD(a)./(opL(a).^2.*sqrt(1./(opL(a).^2)-1)), an(1), am(1));
        end
        
        function out = csch(a)
            an = [a.n];
            am = [a.m];
            out = DS(csch(opV(a)), -coth(opL(a)).*csch(opL(a)).*opD(a), an(1), am(1));
        end
        function out = acsch(a)
            an = [a.n];
            am = [a.m];
            out = DS(acsch(opV(a)), -opD(a)./(opL(a).^2.*sqrt(1./(opL(a).^2)+1)), an(1), am(1));
        end
        
        %% Special Functions
        
        
        
        %%
        
    end
end

