classdef opFunction_1 < opSpot
%opFunction_1   Wrapper for functions.
%
%   opFunction_1(M,N,FUN,CFLAG,LINFLAG) creates a wrapper for function
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
        function op = opFunction_1(m,n,funhandle,cflag,linflag)
            import spot.utils.*
            if nargin < 3
                error('opFunction_1 requires at least three parameters.');
            end
            if nargin < 4 || isempty(cflag)
                cflag = 0;
            end
            if nargin < 5 || isempty(linflag)
                linflag = 1;
            end
            if ~spot.utils.isposintscalar(m)||~spot.utils.isposintscalar(n)
                error('Dimensions of operator must be positive integers.');
            end

            if iscell(funhandle) && length(funhandle) == 2
                if ~isa(funhandle{1},'function_handle') || ...
                   ~isa(funhandle{2},'function_handle')
                    error('Invalid function handle specified.');
                end
                fun = @(x,mode) opFunction_1_intrnl(funhandle,x,mode);

            elseif isa(funhandle,'function_handle')
                fun = @(x,mode) funhandle(x,mode);

            else
                error('Invalid function handle specified.');

            end

            % Construct operator
            op = op@opSpot('Function',m,n);
            op.cflag     = cflag;
            op.linear    = linflag;
            op.funHandle = fun;
            op.sweepflag = true;
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
            
            [unitOut, ~] = opFunction_1.unitGet(tempContainer,mode);
            h.unit = unitOut;
            h.varUnits = unitOut;
            
            [deltaOut,~] = opFunction_1.deltaGet(tempContainer,mode);
            h.delta = deltaOut;
            
            [offsetOut,~] = opFunction_1.originGet(tempContainer,mode);
            h.origin = offsetOut;
            
            if ~isempty(h.IDHistory) % if the operation history is not empty
                prevHis = peek(h.IDHistory);
                
                prevHisID = SDCpckg.utils.getInfoFromCell(prevHis,'ID');
                
                if ~strcmp(prevHisID,op.ID) %if IDs are different
                    push(h.IDHistory,{...
                        {'ID',op.ID},...
                        {'mode',mode},...
                        {'origin',header.origin}});
                else % if IDs are the same
                    prevHisMode = SDCpckg.utils.getInfoFromCell(prevHis,'mode');
                    
                    if strcmp(prevHisMode,mode) % if modes are the same
                        push(h.IDHistory,{...
                            {'ID',op.ID},...
                            {'mode',mode},...
                            {'origin',header.origin}});
                    else % if modes are opposite
                        pop(h.IDHistory); % cancel out
                        
                        prevHisOrigin = SDCpckg.utils.getInfoFromCell(prevHis,'origin');
                        
                        h.origin = prevHisOrigin;
                    end
                end
            else % if history is empty ...
                push(h.IDHistory,{...
                    {'ID',op.ID},...
                    {'mode',mode},...
                    {'origin',header.origin}});
            end
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

function y = opFunction_1_intrnl(funhandle,x,mode)
    if mode == 1
        fun = funhandle{1};
    else
        fun = funhandle{2};
    end

    % Evaluate the function
    y = fun(x);
end % opFunction_1_internl
