classdef opZeros < opSpot
%OPZEROS   Operator equivalent to zeros function.
%
%   opZeros(M,N) creates an operator corresponding to an M-by-N matrix
%   of zeros. If parameter N is omitted it is set to M.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - Public
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % opZeros  constructor.
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = opZeros(varargin)
          if nargin == 0
             m = 1; n = 1;
          elseif nargin == 1
             if length(varargin{1}) == 2
                m = varargin{1}(1);
                n = varargin{1}(2);
             else
                m = varargin{1};
                n = m;
             end
          elseif nargin == 2
             m = varargin{1};
             n = varargin{2};
          else
             error('Too many input arguments.');
          end
          op = op@opSpot('Zeros',m,n);
          op.sweepflag  = true;
        end % function opZeros
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Double
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function A = double(op)
          A = zeros(size(op));
       end
       
    end % methods - public
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - protected
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = protected )
       
        % Multiplication
        function y = multiply(op,x,mode)
           if mode == 1
              s = op.m;
           else
              s = op.n;
           end
   
           if any(isinf(x) | isnan(x))
              y = ones(s,1, class(x)) * NaN;
           else
              y = zeros(s,1, class(x));
           end
        end % function multiply
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Divide
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = divide(op,b,mode)
            % Sweepable
            x = matldivide(op,b,mode);
        end % divide
      
    end % methods - protected
        
end % classdef
