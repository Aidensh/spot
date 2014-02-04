classdef opCTranspose_2 < opSpot_2
%opCTranspose_2   Conjugate transpose of an operator.
%
%   opCTranspose_2(OP) returns the conjugate tranpose of OP.
%
%   See also opTranspose, opConj, opReal, opImag.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = opCTranspose_2(A)
            if nargin ~= 1
                error('Exactly one operator must be specified.')
            end

            % Input matrices are immediately cast as opMatrix's.
            if isa(A,'numeric'), A = opMatrix(A); end

            % Check that the input operators are valid.
            if ~isa(A,'opSpot_2')
                error('Input operator is not valid.')
            end

            % Construct operator
            if A.activated
                [m, n] = size(A);
            else
                m = nan;
                n = nan;
            end
            op = op@opSpot_2('CTranspose_2', n, m);
            op.cflag      = A.cflag;
            op.linear     = A.linear;
            op.sweepflag  = true;
            op.children   = {A};
            op.ms         = A.ns;
            op.ns         = A.ms;
            op.activated  = A.activated;
        end % function opCTranspose_2
        
        function dim = takesDim(op,mode)
            % returns the number of dimension the operator operates on the
            % dataContainer.
            if mode == 1
                dim = takesDim(op.children{1},2);
            else
                dim = takesDim(op.children{1},1);
            end
        end
        
        function op = activateOp(op,header,mode)
            if mode == 1
                op.children{1} = activateOp(op.children{1},header,2);
            else % mode == 2
                op.children{1} = activateOp(op.children{1},header,1);
            end
            % things miss out before super-class constructor
            [theM, theN] = size(op.children{1});
            
            % things miss out during super-class constructor
            m = max(0,theN); % notice it's swapped!!
            n = max(0,theM);
            if round(m) ~= m || round(n) ~= n
                warning('SPOT:ambiguousParams',...
                    'Size parameters are not integer.');
                m = floor(m);
                n = floor(n);
            end
            op.m    = m;
            op.n    = n;
            
            % things to correct after super-class constructor
            op.ms   = op.children{1}.ns;
            op.ns   = op.children{1}.ms;
            
            op.activated = true;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            op1 = op.children{1};
            str = char(op1);
            if op1.precedence > op.precedence
                str = ['(', str, ')'];
            end
            str = [str ,''''];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Conj
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function opOut = conj(op)
            opOut = transpose(op.children{1});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CTranspose
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function opOut = ctranspose(op)
            opOut = op.children{1};
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Transpose
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function opOut = transpose(op)
            opOut = conj(op.children{1});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Drandn
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = drandn(op)
            y = rrandn(op.children{1});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Rrandn
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = rrandn(op)
            y = drandn(op.children{1});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Dzeros
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = dzeros(op)
            y = rzeros(op.children{1});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Rzeros
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = rzeros(op)
            y = dzeros(op.children{1});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % headerMod
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function h = headerMod(op,header,mode)
            A = op.children{1};
            if mode == 1
                h = headerMod(A,header,2);
            else
                h = headerMod(A,header,1);
            end
        end % headerMod
       
    end % methods - public

    methods( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            A = op.children{1};
            if mode == 1
                y = applyMultiply(A,x,2);
            else
                y = applyMultiply(A,x,1);
            end
        end % function multiply

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Divide
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = divide(op,x,mode)
            A = op.children{1};
            if mode == 1
                y = divide(A,x,2);
            else
                y = divide(A,x,1);
            end
        end % function divide

    end % methods - protected
end % classdef
