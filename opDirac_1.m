classdef opDirac_1 < opOrthogonal   
%opDirac_1  Dirac basis.
%   opDirac_1(N) creates the square N-by-N identity operator. Without
%   any arguments an operator corresponding to the scalar 1 is
%   created.
%
%   opDirac_1([N1 N2 ...]) creates the square (N1*N2*...)-by-(N1*N2*...) identity operator
%   created. Collasping dimension case.
%
%   opDirac_1({C,dim1,dim2,...}) creates the square (N1*N2*...)-by-(N1*N2*...) identity operator
%   created. Collasping dimension case.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - public
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        % Constructor
        function op = opDirac_1(n)
            if ~(nargin == 1)
                error('opDirac_1 needs exactly 1 input arguement')
            end
            
            if ~iscell(n) && isscalar(n) && ~isa(n,'SeisDataContainer') % opDirac_1(N)
                n_out = n;
                mns = {n};
            elseif ~iscell(n) && ~isscalar(n) && ~isa(n,'SeisDataContainer') % opDirac_1([N1 N2 ...])
                n_out = prod(n);
                
                mns = cell(1,length(n));
                for ind = 1 : length(n)
                    mns{ind} = n(ind);
                end
            elseif SDCpckg.utils.isForContainerInfo(n) % opDirac_1({C,dim1,dim2,...})
                refDim = n; 
                refDim(1) = []; % remove the dataContainer
                refDim = spot.utils.uncell(refDim); % make the cell array into MATLAB array
                
                if ~isempty(refDim) % make sure there are dimensions selected
                    if ~all(refDim <= length(size(n{1}))) % incase over-selecting
                        error('invalid dimension selected for data container')
                    end
                    theNs = size(n{1},refDim);
                else
                    error('invalid dataContainer refrerncing dimension info'); % please make sure exactly 1 dimension of C is selected
                end
                % from this point, same as previous if condition
                n_out = prod(theNs);
                
                mns = cell(1,length(theNs));
                for ind = 1 : length(theNs)
                    mns{ind} = theNs(ind);
                end
            elseif isa(n,'SeisDataContainer')
                error('Remember to wrap the data container and all the dimension selection arguements together into one cell')
            else
                error('Invalid input type of the first arguement for opDirac_1');
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
            disp('How thoughtful of you to test opDirac_1!!!');
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
end % opDirac_1