classdef (HandleCompatible) opOrthogonal_2 < opSpot_2
%opOrthogonal_2   Abstract class for orthogonal operators.
%
%   opOrthogonal_2 methods:
%     opOrthogonal_2 - constructor
%     mldivide     - solves Ax=b  via  x=A'b.

%   NOTE: There's no reason to overload @opSpot/mrdivide because it's simply
%   a wrapper to mldivide.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot
    
   methods

        function op = opOrthogonal_2(type,m,n)
            % opOrthogonal_2   Constructor for the abstract class. 
            op = op@opSpot_2(type,m,n);
        end         

        function x = mldivide(op,b)
            % \ (backslash)  x = op\b
            x = op'*b;
        end

   end % methods
      
end % opOrthogonal_2