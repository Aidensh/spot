classdef opFoG_2 < opSpot_2
%OPFOG   Forms the product of to operators.
%
%   opFoG_2(OP1,OP2) creates an operator that successively applies each
%   of the operators OP1, OP2 on a given input vector. In non-adjoint
%   mode this is done in reverse order.
%
%   The inputs must be either Spot operators or explicit Matlab matrices
%   (including scalars).
%
%   See also opDictionary, opStack, opSum.

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
        function op = opFoG_2(A,B)
            if nargin ~= 2
                error('Exactly two operators must be specified.')
            end
            
            activated = A.activated && B.activated;
            
            if activated
                % Input matrices are immediately cast as opMatrix's.
                if isa(A,'numeric'), A = opMatrix(A); end
                if isa(B,'numeric'), B = opMatrix(B); end
                
                % Check that the input operators are valid.
                if ~( (isa(A,'opSpot_2') || isa(A,'opSpot')) && (isa(B,'opSpot_2') || (isa(B,'opSpot'))))
                    error('One of the operators is not a valid input.')
                end
                
                % Check operator consistency and complexity
                [mA, nA]   = size(A);
                [mB, nB]   = size(B);
                compatible = isscalar(A) || isscalar(B) || nA == mB;
                if ~compatible
                    error('Operators are not compatible in size.');
                end
                
                % Determine size
                
                if isscalar(A) || isscalar(B)
                    m = max(mA,mB);
                    n = max(nA,nB);
                else
                    m = mA;
                    n = nB;
                end
            else % if at least one of A or B is inactivated, the whole thing is inactivated
                m = nan;
                n = nan;
                activated = false;
            end
            % Construct operator
            op = op@opSpot_2('FoG_2', m, n);
            
            if activated
                op.ms         = A.ms;
                op.ns         = B.ns;
                
                % Preprocess children
                if isscalar(A), op.children{1} = opMatrix(double(A)); end
                if isscalar(B), op.children{2} = opMatrix(double(B)); end
            end
            op.cflag      = A.cflag  | B.cflag;
            op.linear     = A.linear | B.linear;
            op.sweepflag  = A.sweepflag & B.sweepflag;
            op.children   = {A, B};
            op.precedence = 3;
            op.activated  = activated;
        end % Constructor

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % double
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function A = double(op)
            C1 = op.children{1};
            C2 = op.children{2};
            A  = double(C1)*double(C2);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % drandn
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function A = drandn(op,varargin)
            A = drandn(op.children{2},varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % rrandn
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function A = rrandn(op,varargin)
            A = rrandn(op.children{1},varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            % Get children
            op1 = op.children{1};
            op2 = op.children{2};

            % Format first operator
            str1 = char(op1);
            if op1.precedence > op.precedence
                str1 = ['(',str1,')'];
            end

            % Format second operator
            str2 = char(op2);
            if op2.precedence > op.precedence
                str2 = ['(',str2,')'];
            end

            % Combine
            str = [str1, ' * ', str2];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % headerMod
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function h = headerMod(op,header,mode)
            
            nChildren = 2;
            % check history
            usedHistory = false;
            if ~isempty(header.IDHistory)
                prevHis = peek(header.IDHistory);
                prevHisID = prevHis.ID;
                prevHisMode = prevHis.mode;
                prevHisChildren = prevHis.children;
            end
            
            if ~isempty(header.IDHistory) && strcmp(prevHisID,op.ID) && ~strcmp(prevHisMode,mode)
                % if has history, same ID, opposite mode
                usedHistory = true;
                
                for ind = 1:nChildren
                    header.IDHistory = push(header.IDHistory,ContainerStack.toCells(prevHisChildren{ind}));
                end
            end
            
            
            if mode == 1
                numHistory_Pre = length(header.IDHistory.data);
                g = headerMod(op.children{2},header,mode);
                h = headerMod(op.children{1},g,mode);
                numHistory_Post = length(header.IDHistory.data);
                
                if ~usedHistory
                    if numHistory_Pre ~= numHistory_Post
                        childrenOpHistory = cell(1,2);
                        childrenOpHistory{2} = peek(h.IDHistory);
                        h.IDHistory = remove(h.IDHistory,1);
                        childrenOpHistory{1} = peek(h.IDHistory);
                        h.IDHistory = remove(h.IDHistory,1);
                        
                        h.IDHistory = push(h.IDHistory,{...
                            {'ID',op.ID},...
                            {'ClassName',class(op)},...
                            {'mode',mode},...
                            {'opM',op.m},...
                            {'opN',op.n},...
                            {'children',childrenOpHistory}});
                    end
                else % if used history; in other words, if cancel out
                    topHistory = peek(h.IDHistory);
                    
                    while topHistory.ID ~= op.ID
                        h.IDHistory = remove(h.IDHistory,1);
                        topHistory = peek(h.IDHistory);
                    end
                    h.IDHistory = remove(h.IDHistory,1); % remove it self
                end
                
            else
                numHistory_Pre = length(header.IDHistory.data);
                g = headerMod(op.children{1},header,mode);
                h = headerMod(op.children{2},g,mode);
                numHistory_Post = length(header.IDHistory.data);
                
                if ~usedHistory
                    if numHistory_Pre ~= numHistory_Post
                        childrenOpHistory = cell(1,2);
                        childrenOpHistory{1} = peek(h.IDHistory);
                        h.IDHistory = remove(h.IDHistory,1);
                        childrenOpHistory{2} = peek(h.IDHistory);
                        h.IDHistory = remove(h.IDHistory,1);
                        
                        
                        childrenOpHistory = fliplr(childrenOpHistory);
                        
                        h.IDHistory = push(h.IDHistory,{...
                            {'ID',op.ID},...
                            {'ClassName',class(op)},...
                            {'mode',mode},...
                            {'opM',op.m},...
                            {'opN',op.n},...
                            {'children',childrenOpHistory}});
                    end
                else % if used history; in other words, if cancel out
                    topHistory = peek(h.IDHistory);
                    
                    while topHistory.ID ~= op.ID
                        h.IDHistory = remove(h.IDHistory,1);
                        topHistory = peek(h.IDHistory);
                    end
                    h.IDHistory = remove(h.IDHistory,1); % remove it self
                end
            end
        end % headerMod
        
        
        function op = activateOp(op,header,mode)
            % activate the operator right before it operates on container
            % header : the header of the data container
            
            hTemp = header;
            indStart = 1; % where to start reading container's implicit dimension
            indEnd = indStart;
            
            nChildren = 2;
            
            % check history
            usedHistory = false;
            if ~isempty(hTemp.IDHistory)
                prevHis = peek(hTemp.IDHistory);
                prevHisID = prevHis.ID;
                prevHisMode = prevHis.mode;
                prevHisChildren = prevHis.children;
            end
            
            if ~isempty(hTemp.IDHistory) && strcmp(prevHisID,op.ID) && ~strcmp(prevHisMode,mode)
                % if has history, same ID, opposite mode
                usedHistory = true;
                
                hTemp.IDHistory = remove(hTemp.IDHistory,1);
                
                for ind = 1:nChildren
                    hTemp.IDHistory = push(hTemp.IDHistory,ContainerStack.toCells(prevHisChildren{ind}));
                end
            end
            % check history / end
            
            % Activate Children operators
            A = op.children{1};
            B = op.children{2};
            
            if mode == 1
                B = activateOp(B,hTemp,mode);
                
                hTemp = headerMod(B,hTemp,mode);
                
                A = activateOp(A,hTemp,mode);
                
            else
                A = activateOp(A,hTemp,mode);
                
                hTemp = headerMod(A,hTemp,mode);
                
                B = activateOp(B,hTemp,mode);
                
            end
            
            % parts missing before super-class constructor
            %-- Input matrices are immediately cast as opMatrix's.
            if isa(A,'numeric'), A = opMatrix(A); end
            if isa(B,'numeric'), B = opMatrix(B); end
            
            %-- Check that the input operators are valid.
            if ~( (isa(A,'opSpot_2') || isa(A,'opSpot')) && (isa(B,'opSpot_2') || (isa(B,'opSpot'))))
                error('One of the operators is not a valid input.')
            end
            
            %-- Check operator consistency and complexity
            [mA, nA]   = size(A);
            [mB, nB]   = size(B);
            compatible = isscalar(A) || isscalar(B) || nA == mB;
            if ~compatible
                error('Operators are not compatible in size.');
            end
            
            %-- Determine size
            
            if isscalar(A) || isscalar(B)
                m = max(mA,mB);
                n = max(nA,nB);
            else
                m = mA;
                n = nB;
            end
            % parts missing during super-class constructor
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
            
                op.children{1} = A;
                op.children{2} = B;
            % parts missing after super-class constructor
            op.ms         = A.ms;
            op.ns         = B.ns;
            
            % Preprocess children
            if isscalar(A), op.children{1} = opMatrix(double(A)); end
            if isscalar(B), op.children{2} = opMatrix(double(B)); end
            
            % activation complete
            op.activated = true; 
        end
        
        function dim = takesDim(op,mode)
            dim = takesDim(op.children{1},mode) + takesDim(op.children{2},mode);
        end

    end % public Methods
       
 
    methods ( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function z = multiply(op,x,mode)
            if mode == 1
                y = applyMultiply(op.children{2},x,mode);
                z = applyMultiply(op.children{1},y,mode);
            else % mode = 2
                y = applyMultiply(op.children{1},x,mode);
                z = applyMultiply(op.children{2},y,mode);
            end
        end % Multiply
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Divide
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = divide(op,b,mode)
            % Needs sweepflag
            if op.sweepflag
                x = matldivide(op,b,mode);
            else
                x = lsqrdivide(op,b,mode);
            end
        end % divide
        
    end % Methods
   
end % Classdef
