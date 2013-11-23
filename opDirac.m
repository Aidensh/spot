classdef opDirac < opOrthogonal   
%OPDIRAC  Dirac basis.
%
%   opDirac(N) creates the square N-by-N identity operator. Without
%   any arguments an operator corresponding to the scalar 1 is
%   created.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - public
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        % Constructor
        function op = opDirac(n)
            if nargin < 1, n = 1; end
            
            if isscalar(n)
                n_out = n;
                mns = {n};
            else
                n_out = prod(n);
                
                mns = cell(1,length(n));
                for ind = 1 : length(n)
                    mns{ind} = n(ind);
                end
            end
            op = op@opOrthogonal('Dirac',n_out,n_out);
            
            op.ms = mns;
            op.ns = mns;
            
            op.isDirac   = true;
            op.sweepflag = true;
        end % constructor

        function A = double(op)
            A = eye(size(op));
        end % double

        function result = xtratests(op)
        %XTRATESTS    User defined tests
        %
        % Just a demo here
            result = true;
            disp('How thoughtful of you to test opDirac!!!');
        end % xtratests

    end % methods - public

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - protected
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiplication
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = multiply(op,x,mode)
            % Nothing to do here
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Divide
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = divide(op,x,mode)
            % Nothing to do here
        end % divide
    end % methods - protected
end % opDirac