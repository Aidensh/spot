classdef opFunction_2 < opSpot_2
%opFunction_2   Wrapper for functions.
%
%   opFunction_2(M,N,FUN,CFLAG,LINFLAG) creates a wrapper for function
%   FUN, which corresponds to an M-by-N operator. The FUN parameter
%   can be one of two types:
%
%   1) A handle to a function of the form FUN(X,MODE), where the
%      operator is applied to X when MODE = 1, and the transpose is
%      applied when MODE = 2;
%   2) A cell array of two function handles: {FUN,FUN_TRANSPOSE},
%      each of which requires only one parameter, X.
%
%   Optional arguments CFLAG and LINFLAG indicate whether the
%   function implements a complex or real operator and whether it
%   is linear or not. The default values are CFLAG=0, LINFLAG=1.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( SetAccess = private )
       funHandle  % Function handles
    end % Properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - Public
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        % Constructor
        %function op = opFunction_2(m,n,funhandle,cflag,linflag)
        function op = opFunction_2(varargin)
            activated = true;
            
            import spot.utils.*
            if nargin < 1
                error('opFunction_2 requires at least a funhandle parameter.');
            end
            
            % first input argument
            if isa(varargin{1},'function_handle') %%opFunction_2(funhandle,...)
                % inactive operator
                m = nan;
                n = nan;
            elseif SDCpckg.utils.isForContainerInfo(varargin{1}) %opFunction_2({C,dim1,dim2},funhandle,...)
                C = varargin{1}{1};
                dim = varargin{1}{2};
                
                m = size(C,dim(1));
                n = size(C,dim(2));
                
                varargin(1) = [];
            elseif iscell(varargin{1}) && isempty(varargin{1}) %opFunction_2({},funhandle,...)
                % inactive operator
                m = nan;
                n = nan;
                
                varargin(1) = [];
            elseif isnumeric(varargin{1}) && length(varargin{1}) == 1 %opFunction_2(m,n,funhandle,...)
                m = varargin{1};
                n = varargin{2};
                
                varargin(2) = [];
                varargin(1) = [];
            elseif isa('SeisDataContainer',varargin{1})
                error('Remember to wrap dataContainer and dimension info into cell')
            else
                error('something weird happened in opFunction')
            end
            
            nargs = length(varargin);
            if nargs < 2 || isempty(varargin{2})
                cflag = 0;
            end
            if nargs < 3 || isempty(varargin{3})
                linflag = 1;
            end
            
            if activated
                if ~spot.utils.isposintscalar(m)||~spot.utils.isposintscalar(n)
                    error('Dimensions of operator must be positive integers.');
                end
            end
            
            funhandle = varargin{1};
            if iscell(funhandle) && length(funhandle) == 2
                if ~isa(funhandle{1},'function_handle') || ...
                   ~isa(funhandle{2},'function_handle')
                    error('Invalid function handle specified.');
                end
                fun = @(x,mode) opFunction_2_intrnl(funhandle,x,mode);

            elseif isa(funhandle,'function_handle')
                fun = @(x,mode) funhandle(x,mode);

            else
                error('Invalid function handle specified.');

            end

            % Construct operator
            op = op@opSpot_2('Function_2',m,n);
            op.cflag     = cflag;
            op.linear    = linflag;
            op.funHandle = fun;
            op.sweepflag = true;
            
            op.activated = activated;
        end % Constructor
        
    end % Methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - protected
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = protected )
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            for u = size(x,2):-1:1 % Loop through multivector
                y(:,u) = op.funHandle(x(:,u),mode);
            end
        end % Multiply
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Divide
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = divide(op,b,mode)
            % Non-sweepable
            x = lsqrdivide(op,b,mode);
        end % divide
    end % Methods
    methods
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
            
            % modify unit, origin, offset
            tempContainer = iCon(header);
            
            [unitOut, ~] = opFunction_2.unitGet(tempContainer,mode);
            h.unit = unitOut;
            h.varUnits = unitOut;
            
            [deltaOut,~] = opFunction_2.deltaGet(tempContainer,mode);
            h.delta = deltaOut;
            
            [offsetOut,~] = opFunction_2.originGet(tempContainer,mode);
            h.origin = offsetOut;
            
            if ~isempty(h.IDHistory) % if the operation history is not empty
                prevHis = peek(h.IDHistory);
                
                prevHisID = prevHis.ID;
                
                if ~strcmp(prevHisID,op.ID) %if IDs are different
                    h.IDHistory = push(h.IDHistory,{...
                        {'ID',op.ID},...
                        {'ClassName',class(op)},...
                        {'mode',mode},...
                        {'opM',op.m},...
                        {'opN',op.n},...
                        {'opMs',op.ms},...
                        {'opNs',op.ns},...
                        {'origin',header.origin}});
                else % if IDs are the same
                    prevHisMode = prevHis.mode;
                    
                    if strcmp(prevHisMode,mode) % if modes are the same
                        h.IDHistory = push(h.IDHistory,{...
                            {'ID',op.ID},...
                            {'ClassName',class(op)},...
                            {'mode',mode},...
                            {'opM',op.m},...
                            {'opN',op.n},...
                            {'opMs',op.ms},...
                            {'opNs',op.ns},...
                            {'origin',header.origin}});
                    else % if modes are opposite
                        h.IDHistory = remove(h.IDHistory,1); % cancel out
                        
                        prevHisOrigin = prevHis.origin;
                        
                        h.origin = prevHisOrigin;
                    end
                end
            else % if history is empty ...
                h.IDHistory = push(h.IDHistory,{...
                    {'ID',op.ID},...
                    {'ClassName',class(op)},...
                    {'mode',mode},...
                    {'opM',op.m},...
                    {'opN',op.n},...
                    {'opMs',op.ms},...
                    {'opNs',op.ns},...
                    {'origin',header.origin}});
            end
        end
        
        function op = activateOp(op,header,~)
            % activate the operator right before it operates on container
            % header : the header of the data container
            
            theSize = header.size;
            if any(strcmp('exsize',fieldnames(header)))
                m = prod(theSize(header.exsize(1,1):header.exsize(2,1)));
                n = prod(theSize(header.exsize(1,2):header.exsize(2,2)));
            else
                warning('exsize field not exist in header')
                m = prod(theSize); % assume container is vectorized
                n = 1;
            end
            
            % things missing before super-class constructor
            if ~spot.utils.isposintscalar(m)||~spot.utils.isposintscalar(m)
                error('Dimensions of operator must be positive integers.');
            end
            
            % things missing inside super-class constructor
            m = max(0,m);
            n = max(0,n);
            if round(m) ~= m || round(n) ~= n
                warning('SPOT:ambiguousParams',...
                    'Size parameters are not integer.');
                m = floor(m);
                n = floor(n);
            end
            op.m    = m;
            op.n    = n;
            op.ms   = {m};
            op.ns   = {n};
            
            op.activated = true;
        end
        
        function dim = takesDim(op,mode)
            % returns the number of dimension the operator operates on the
            % dataContainer.
            dim = 1;
        end
    end
    
    methods( Static)
        function [unitOut, unitIn] = unitGet(C,~)
            unitIn = unit(C);
            unitOut = unitIn;
        end
        function [originOut, originIn] = originGet(C,~)
            originIn = origin(C);
            originOut = originIn;
        end
        function [deltaOut, deltaIn] = deltaGet(C,~)
            deltaIn = delta(C);
            deltaOut = deltaIn;
        end
    end
    
    
        
end % Classdef

%======================================================================

function y = opFunction_2_intrnl(funhandle,x,mode)
    if mode == 1
        fun = funhandle{1};
    else
        fun = funhandle{2};
    end

    % Evaluate the function
    y = fun(x);
end % opFunction_2_internl
