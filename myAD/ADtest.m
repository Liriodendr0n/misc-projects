classdef ADtest
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x
        dx
    end
    
    methods
        %% constructor
        function out = ADtest(varargin)
            if nargin == 1
                out = ADtest(varargin{1}.x, varargin{1}.dx);
            elseif nargin == 2
                out.x = varargin{1};
                out.dx = varargin{2};
            end
        end
        
        %% unary operators
        
        %% binary operators
        function out = plus(a, b)
            if isa(a, 'float')
                a = ADtest(a,0);
            end
            if isa(b, 'float')
                b = ADtest(b,0);
            end
            
%             out.x = a.x + b.x;
%             out.dx = a.dx + b.dx;
            out = ADtest(a.x + b.x, plus(a.dx, b.dx));
        end
        
        function out = minus(a, b)
            if isa(a, 'float')
                a = ADtest(a,0);
            end
            if isa(b, 'float')
                b = ADtest(b,0);
            end
            
            out.x = a.x - b.x;
            out.dx = a.dx - b.dx;
        end

        function out = times(a, b)
            if isa(a, 'float')
                a = ADtest(a,0);
            end
            if isa(b, 'float')
                b = ADtest(b,0);
            end
     
%             out.x = a.x .* b.x;
%             out.dx = a.dx .* b.x + a.x .* b.dx;
            
            out = ADtest(a.x .* b.x, a.dx .* b.x + a.x .* b.dx);
        end
        
        function out = rdivide(a, b)
            if isa(a, 'float')
                a = ADtest(a,0);
            end
            if isa(b, 'float')
                b = ADtest(b,0);
            end
            
            out.x = a.x ./ b.x;
            out.dx = (a.dx .* b.x - a.x .* b.dx) ./ b.x^2;
        end
        
        function out = ldivide(a, b)
            if isa(a, 'float')
                a = ADtest(a,0);
            end
            if isa(b, 'float')
                b = ADtest(b,0);
            end
            out.x = a.x .\ b.x;
            out.dx = (a.dx .* b.x - a.x .* b.dx) .\ b.x^2;
        end
        
        %% unary functions
        function out = sin(a)
%             out.x = sin(a.x);
%             out.dx = cos(a.x) .* a.dx;
            out = ADtest(sin(a.x), cos(a.x) .* a.dx);
        end
        
        function out = cos(a)
%             out.x = cos(a.x);
%             out.dx = -sin(a.x) .* a.dx;
            out = ADtest(cos(a.x), -sin(a.x) .* a.dx);
        end
        
    end
end

