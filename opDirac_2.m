classdef opDirac_2 < opOrthogonal_2   
%opDirac_2  Dirac basis.
%   opDirac_2(N) creates the square N-by-N identity operator. Without
%   any arguments an operator corresponding to the scalar 1 is
%   created.
%
%   opDirac_2([N1 N2 ...]) creates the square (N1*N2*...)-by-(N1*N2*...) identity operator
%   created. Collasping dimension case.
%
%   opDirac_2({C,dim1,dim2,...}) creates the square (N1*N2*...)-by-(N1*N2*...) identity operator
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
        function op = opDirac_2(varargin)
            activated = true;
            if isempty(varargin)
                n = nan;
                mns = {};
                activated = false;
            elseif ~iscell(varargin{1}) && isscalar(varargin{1}) && ~isa(varargin{1},'SeisDataContainer') % opDirac_2(N)
                n = varargin{1};
                mns = {varargin{1}};
            elseif ~iscell(varargin{1}) && ~isscalar(varargin{1}) && ~isa(varargin{1},'SeisDataContainer') % opDirac_2([N1 N2 ...])
                n = prod(varargin{1});
                
                mns = cell(1,length(varargin{1}));
                for ind = 1 : length(varargin{1})
                    mns{ind} = varargin{1}(ind);
                end
            elseif SDCpckg.utils.isForContainerInfo(varargin{1}) % opDirac_2({C,dim1,dim2,...})
                refDim = varargin{1}; 
                refDim(1) = []; % remove the dataContainer
                refDim = spot.utils.uncell(refDim); % make the cell array into MATLAB array
                
                if ~isempty(refDim) % make sure there are dimensions selected
                    if ~all(refDim <= length(size(varargin{1}{1}))) % incase over-selecting
                        error('invalid dimension selected for data container')
                    end
                    theNs = size(varargin{1}{1},refDim);
                else
                    error('invalid dataContainer refrerncing dimension info'); % please make sure exactly 1 dimension of C is selected
                end
                % from this point, same as previous if condition
                n = prod(theNs);
                
                mns = cell(1,length(theNs));
                for ind = 1 : length(theNs)
                    mns{ind} = theNs(ind);
                end
            elseif iscell(varargin{1}) && isempty(varargin{1})
                n = nan;
                mns = {};
                activated = false;
            elseif isa(varargin{1},'SeisDataContainer')
                error('Remember to wrap the data container and all the dimension selection arguements together into one cell')
            else
                error('Invalid input type of the first arguement for opDirac_2');
            end
            op = op@opOrthogonal_2('Dirac_2',n,n);
            
            op.ms = mns;
            op.ns = mns;
            
            op.isDirac   = true;
            op.sweepflag = true;
            
            op.activated = activated;
        end % constructor
        
        
        function h = headerMod(op,header,mode)
            exsize = header.exsize;
            
            if mode == 1
                h = header; % Copy header
                % Replace old first (collapsed) dimensional sizes with operator sizes.
                h.size(exsize(1,1):exsize(2,1)) = [];
                h.size = [op.ms{:} h.size];
            else
                h = header;
                h.size(exsize(1,1):exsize(2,1)) = [];
                h.size = [op.ns{:} h.size];
            end
            
            exsize_out = 1:length(h.size);
            exsize_out = [exsize_out;exsize_out];
            h.exsize   = exsize_out;
            
            if ~isempty(h.IDHistory) % if the operation history is not empty
                prevHis = peek(h.IDHistory);
                
                prevHisID = prevHis.ID';
                
                if ~strcmp(prevHisID,op.ID) %if IDs are different
                    h.IDHistory = push(h.IDHistory,{...
                        {'ID',op.ID},...
                        {'ClassName',class(op)},...
                        {'mode',mode},...
                        {'opM',op.m},...
                        {'opN',op.n},...
                        {'children',{}}});
                else % if IDs are the same
                    prevHisMode = prevHis.mode;
                    
                    if strcmp(prevHisMode,mode) % if modes are the same
                        h.IDHistory = push(h.IDHistory,{...
                            {'ID',op.ID},...
                            {'ClassName',class(op)},...
                            {'mode',mode},...
                            {'opM',op.m},...
                            {'opN',op.n},...
                            {'children',{}}});
                    else % if modes are opposite
                        h.IDHistory = remove(h.IDHistory,1); % cancel out
                        
                        prevHisOrigin = prevHis.origin;
                        
                        h.origin = prevHisOrigin;
                    end
                end
            else % if it is empty ...
                h.IDHistory = push(h.IDHistory,{...
                    {'ID',op.ID},...
                    {'ClassName',class(op)},...
                    {'mode',mode},...
                    {'opM',op.m},...
                    {'opN',op.n},...
                    {'children',{}}});
            end
            
        end
        
        function op = activateOp(op,header,~)
            % activate the operator right before it operates on container
            % header : the header of the data container
            
            n = header.size;
            if any(strcmp('exsize',fieldnames(header)))
                n = prod(n(header.exsize(1,1):header.exsize(2,1)));
            else
                warning('exsize field not exist in header')
                n = prod(n); % assume container is vectorized
            end
            
            % from opSpot's constructor.
            % make up what the operator missed during it was initialized.
            n = max(0,n);
            if round(n) ~= n
                warning('SPOT:ambiguousParams',...
                    'Size parameters are not integer.');
                n = floor(n);
            end
            op.m    = n;
            op.n    = n;
            op.ms   = {n};
            op.ns   = {n};
            
            op.activated = true; % activation complete
        end
        
        function dim = takesDim(op,~)
            % returns the number of dimension the operator operates on the
            % dataContainer.
            dim = 1;
        end
        
        function A = double(op)
            A = eye(size(op));
        end % double

        function result = xtratests(op)
        %XTRATESTS    User defined tests
        %
        % Just a demo here
            result = true;
            disp('How thoughtful of you to test opDirac_2!!!');
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
end % opDirac_2