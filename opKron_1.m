classdef opKron_1 < opSpot
%opKron_1   Kronecker tensor product.
%
%   opKron_1(OP1,OP2,...OPn) creates an operator that is the Kronecker
%   tensor product of OP1, OP2, ..., OPn.

%   Copyright 2009, Rayan Saab, Ewout van den Berg and 
%   Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        permutation; %Permutation vector of intergers defining the order to
        %use when the operators (children) of the Kronecker product are
        %applied to a data vector.
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Public methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = opKron_1(varargin)
            narg=nargin;
            
            %Test the case where varargin is a list
            if narg == 1
                narg=length(varargin{1});
                varargin=varargin{1};
            end
            
            if narg < 2
                error('At least two operators must be specified')
            end
            
            % Input matrices are immediately cast to opMatrix.
            for i=1:narg
                if isa(varargin{i},'numeric')
                    varargin{i} = opMatrix(varargin{i});
                elseif ~isa(varargin{i},'opSpot')
                    error('One of the operators is not a valid input.')
                end
            end
            
            % Determine operator size and complexity (this code is
            % general for any number of operators)
            opA       = varargin{1};
            [m,n]     = size(opA);
            cflag     = opA.cflag;
            linear    = opA.linear;
            sweepflag = opA.sweepflag;
            
            for i=2:narg
                opA       = varargin{i};
                cflag     = cflag  | opA.cflag;
                linear    = linear & opA.linear;
                sweepflag = sweepflag & opA.sweepflag;
                [mi,ni]   = size(opA);
                m = m * mi; n = n * ni;
            end
            
            % Construct operator
            op = op@opSpot('Kron_1', m, n);
            op.cflag       = cflag;
            op.linear      = linear;
            op.sweepflag   = sweepflag;
            op.children    = varargin;
            op.permutation =(1:narg);
            
            %Evaluate the best permutation to use when a multiplication is
            %applied
            if ~ (m == 0 || n == 0)
                op.permutation=op.best_permutation();
            end
            
            % Setting up implicit dimensions of output vector
            % Flipped
            varargin = fliplr(varargin);
            len      = length(varargin);
            op.ms    = cell(1,len);
            op.ns    = cell(1,len);
            for u = 1:len
                child_op = varargin{u};
                if length(child_op.ms) > 1
                    op.ms{u} = [op.ms{u} child_op.ms(:)'];
                else
                    op.ms{u} = [op.ms{u} child_op.ms{:}];
                end
                if length(child_op.ns) > 1
                    op.ns{u} = [op.ns{u} child_op.ns(:)'];
                else
                    op.ns{u} = [op.ns{u} child_op.ns{:}];
                end
            end
                
        end % Constructor
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function str = char(op)
            str=['Kron(',char(op.children{1})];
            
            % Get operators
            for i=2:length(op.children)
                str=strcat(str,[', ',char(op.children{i})]);
            end
            str=strcat(str,')');
        end % Char
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % headerMod
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function h = headerMod(op,header,mode)
            
            imSize = header.size;
            childrenList = fliplr(op.children); % opKron_1 operates in reverse order
            
            startExWith = 1; % exsize start
            endExWith = 1; % exsize end
            onOp = 1; % working on which child operator
            outSize = []; % this goes to h.size, which will be the new implicit size of the vector.
            onDim = 1; % which output dimension
            
            newHeader = header;
            if mode == 1
                
                while endExWith <= length(imSize)
                    matchSize = prod(imSize(startExWith:endExWith));
                    
                    currentChild = childrenList{onOp};
                    currentChildOpNs = spot.utils.uncell(currentChild.ns);
                    
                    if matchSize == prod(currentChildOpNs) && isscalar(currentChild.ms{1}) && isscalar(currentChild.ns{1}) % if the size matched and has no collasping dimension
                        for ind = 1 : length(currentChildOpNs) % check to see if all implicit dimension match with the operator dimension
                            if ~(currentChild.ns{ind} == imSize(ind + startExWith - 1))
                                error('implicit dimension do not match')
                            end
                            outSize(onDim) = currentChild.ms{ind};
                            onDim = onDim + 1;
                        end
                        
                        corrHeader = spot.data.headerRef(header,startExWith:endExWith);
                        childHeaderMod = headerMod(currentChild,corrHeader,mode);
                        newHeader = spot.data.headerAsgn(newHeader,childHeaderMod,startExWith:endExWith);
                        
                        if onOp == length(childrenList) && length(imSize) > endExWith % if reached last children operator and there are still more implicit dimensions not operated ...
                            startExWith = endExWith + 1;
                            endExWith = length(imSize);
                            
                            diracDim = prod(imSize(startExWith:endExWith)); % handle the remaining implicit dimensions by operating with opDirac
                            outSize(onDim) = diracDim;
                            onDim = onDim + 1;
                        end
                        
                        % prepare for the next round operator
                        startExWith = endExWith + 1;
                        endExWith = startExWith;
                        onOp = onOp + 1;
                    elseif matchSize == prod(currentChildOpNs) &&  length(currentChild.ns(1)) < length(currentChild.ns{1})% matched. Collasping dimension
                        for ind = 1 : length(currentChild.ns{1})
                            if ~(imSize(onDim + ind - 1) == currentChild.ns{1}(ind))
                                error('Size of Collasping Dimension do not match')
                            end
                        end
                        
                        outSize(onDim) = currentChild.ms{:};
                        onDim = onDim + length(currentChild.ns{1}) -1;
                        
                        corrHeader = spot.data.headerRef(header,startExWith:endExWith);
                        childHeaderMod = headerMod(currentChild,corrHeader,mode);
                        newHeader = spot.data.headerAsgn(newHeader,childHeaderMod,startExWith:endExWith);
                        
                        if onOp == length(childrenList) && length(imSize) > endExWith % if reached last children operator and there are still more implicit dimensions not operated ...
                            startExWith = endExWith + 1;
                            endExWith = length(imSize);
                            
                            diracDim = prod(imSize(startExWith:endExWith)); % handle the remaining implicit dimensions by operating with opDirac
                            outSize(onDim) = diracDim;
                            onDim = onDim + 1;
                        end
                        
                        % prepare for the next round operator
                        startExWith = endExWith + 1;
                        endExWith = startExWith;
                        onOp = onOp + 1;
                    elseif length(currentChild.ns(1)) == length(currentChild.ns{1}) && currentChild.ns{1}== imSize(startExWith) && length(currentChild.ms(1)) < length(currentChild.ms{1}) % matched. Expanding dimension
                        for ind = 1 : length(currentChild.ms{:})
                            outSize(onDim + ind - 1) = currentChild.ms{:}(ind);
                        end
                        onDim = onDim + length(currentChild.ms{:}) - 1;
                        
                        corrHeader = spot.data.headerRef(header,startExWith:endExWith);
                        childHeaderMod = headerMod(currentChild,corrHeader,mode);
                        newHeader = spot.data.headerAsgn(newHeader,childHeaderMod,startExWith:endExWith);
                        
                        % prepare for the next round operator
                        startExWith = endExWith + 1;
                        endExWith = startExWith;
                        onOp = onOp + 1;
                    elseif endExWith == length(imSize) % if size do not match and is now at the end of the exsize list
                        error('Unable to match data implicit dimension with opKron_1 children size') % give error: the size can never match up
                    else
                        % keep trying :)
                        endExWith = endExWith + 1;
                        
                    end
                end
            else % mode == 2; transpose case
                
                while endExWith <= length(imSize)
                    matchSize = prod(imSize(startExWith:endExWith));
                    
                    currentChild = childrenList{onOp};
                    currentChildOpMs = spot.utils.uncell(currentChild.ms);
                    
                    if matchSize == prod(currentChildOpMs) && isscalar(currentChild.ms{1}) && isscalar(currentChild.ns{1}) % if the size matched and has no collasping dimension
                        for ind = 1 : length(currentChildOpMs)
                            if ~(currentChild.ms{ind} == imSize(ind + startExWith - 1))
                                error('implicit dimension do not match')
                            end
                            outSize(onDim) = currentChild.ns{ind};
                            onDim = onDim + 1;
                        end
                        
                        corrHeader = spot.data.headerRef(header,startExWith:endExWith);
                        childHeaderMod = headerMod(currentChild,corrHeader,mode);
                        newHeader = spot.data.headerAsgn(newHeader,childHeaderMod,startExWith:endExWith);
                        
                        if onOp == length(childrenList) && length(imSize) > endExWith % if reached last children operator and there are still more implicit dimensions not operated ...
                            startExWith = endExWith + 1;
                            endExWith = length(imSize);
                            
                            diracDim = prod(imSize(startExWith:endExWith)); % handle the remaining implicit dimensions by operating with opDirac
                            outSize(onDim) = diracDim;
                            onDim = onDim + 1;
                        end
                        
                        % prepare for the next round operator
                        startExWith = endExWith + 1;
                        endExWith = startExWith;
                        onOp = onOp + 1;
                    elseif matchSize == prod(currentChildOpMs) &&  length(currentChild.ms(1)) < length(currentChild.ms{1})% matched. Collasping dimension
                        for ind = 1 : length(currentChild.ms{1})
                            if ~(imSize(onDim + ind - 1) == currentChild.ms{1}(ind))
                                error('Size of Collasping Dimension do not match')
                            end
                        end
                        
                        outSize(onDim) = currentChild.ns{:};
                        onDim = onDim + length(currentChild.ms{1}) -1;
                        
                        corrHeader = spot.data.headerRef(header,startExWith:endExWith);
                        childHeaderMod = headerMod(currentChild,corrHeader,mode);
                        newHeader = spot.data.headerAsgn(newHeader,childHeaderMod,startExWith:endExWith);
                        
                        
                        if onOp == length(childrenList) && length(imSize) > endExWith % if reached last children operator and there are still more implicit dimensions not operated ...
                            startExWith = endExWith + 1;
                            endExWith = length(imSize);
                            
                            diracDim = prod(imSize(startExWith:endExWith)); % handle the remaining implicit dimensions by operating with opDirac
                            outSize(onDim) = diracDim;
                            onDim = onDim + 1;
                        end

                        % prepare for the next round operator
                        startExWith = endExWith + 1;
                        endExWith = startExWith;
                        onOp = onOp + 1;
                    elseif length(currentChild.ms(1)) == length(currentChild.ms{1}) && currentChild.ms{1}== imSize(startExWith) && length(currentChild.ns(1)) < length(currentChild.ns{1}) % matched. Expanding dimension
                        for ind = 1 : length(currentChild.ns{:})
                            outSize(onDim + ind - 1) = currentChild.ns{:}(ind);
                        end
                        onDim = onDim + length(currentChild.ns{:}) - 1;
                        
                        
                        corrHeader = spot.data.headerRef(header,startExWith:endExWith);
                        childHeaderMod = headerMod(currentChild,corrHeader,mode);
                        newHeader = spot.data.headerAsgn(newHeader,childHeaderMod,startExWith:endExWith);
                        
                        % prepare for the next round operator
                        startExWith = endExWith + 1;
                        endExWith = startExWith;
                        onOp = onOp + 1;
                    elseif endExWith == length(imSize) % if size do not match and is now at the end of the exsize list
                        error('Unable to match data implicit dimension with opKron_1 children size') % give error: the size can never match up
                    else
                        % keep trying :)
                        endExWith = endExWith + 1;
                        
                    end
                end
                
            end
            
            
            n_out_dims = length(outSize);
            % Preallocate and setup header
            header_out        = header;
            header_out.dims   = n_out_dims;
            header_out.size   = zeros(1,n_out_dims);
            header_out.origin = zeros(1,n_out_dims);
            header_out.delta  = zeros(1,n_out_dims);
            header_out.unit   = zeros(1,n_out_dims);
            header_out.label  = zeros(1,n_out_dims);
            
            
            exsize_out = 1:n_out_dims;
            exsize_out = [exsize_out;exsize_out];
            h = header_out;
            h.size = outSize;
            h.origin = newHeader.origin;
            h.delta = newHeader.delta;
            h.unit = newHeader.unit;
            h.label = newHeader.label;
            h.exsize = exsize_out;
        end % headerMod
    end % Methods
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Protected methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods ( Access = protected )
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            
            % The Kronecker product (KP) is applied to the righthand matrix
            % taking in account the best order to apply the operators.
            % That necessitates to decompose the KP in successive matrix
            % products with terms of type I(a) kron A kron I(b).
            % A is the operator to apply. I(a) and I(b) are identity
            % matrices with respective sizes a and b.
            
            opList       = op.children; % Contains list of opKron_1 children
            ncol         = size(x,2); % Number of columns of 'x'
            nbr_children = length(opList); % Number of children
            
            % Pre-registering of the sizes of opKron_1's children
            sizes = zeros(nbr_children,2);
            for i=1:nbr_children
                sizes(i,:) = size(opList{i});
            end
            
            %%%%%%%%%%%%%%%%%%%%%%Multiplication%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if mode == 1 % Classic mode
                perm = op.permutation; % Permutation to take in account.
                m = op.m; % Height of the resulting matrix
                
                for i = 1:nbr_children
                    % Index of the operator A to consider.
                    index = perm(i);
                    
                    % Calculation of the sizes of the identity matrixs used
                    % in the Kronecker product I(a) kron A kron I(b)
                    
                    % Size of I(a)
                    a = 1;
                    for k = 1:(index-1)
                        if i > find(perm==k)
                            a = a * sizes(k,1);
                        else
                            a = a * sizes(k,2);
                        end
                    end
                    
                    % If 'x' has several columns. The initial matrix I(a)
                    % kron A kron I(b) is replicated 'ncol' (number of
                    % columns of x) times) along the diagonal.
                    if ncol>1
                        a=a*ncol;
                    end
                    
                    % Size of I(b)
                    b = 1;
                    for k = (index+1):nbr_children
                        if i > find(perm==k)
                            b = b * sizes(k,1);
                        else
                            b = b * sizes(k,2);
                        end
                    end
                    
                    % Size of the operator A=opList{index} to apply
                    r = sizes(index,1);
                    c = sizes(index,2);
                    
                    % (I(a) kron A kron I(b)) * x;
                    t = reshape(reshape(x,b,a*c).',c,a*b);
                    x = reshape(applyMultiply(opList{index},t,1)',a,r*b)';
                end
                y = reshape(x,m,ncol);
                
            elseif mode == 2 % Transpose mode
                perm = op.permutation(length(opList):-1:1); % The
                % permutation has to be in the other direction since with
                % transposition, operators' computational costs will be
                % inverted.
                n = op.n; % Height of the resulting matrix
                
                for i = 1:nbr_children
                    %Index of the operator A to consider.
                    index = perm(i);
                    
                    % Calculation of the sizes of the identity matrixs used
                    % in the Kronecker product I(a) kron A kron I(b)
                    
                    %Size of I(a)
                    a = 1;
                    for k = 1:(index-1)
                        if i > find(perm==k)
                            a = a * size(opList{k},2);
                        else
                            a = a * size(opList{k},1);
                        end
                    end
                    
                    % If 'x' has several columns. The initial matrix I(a)
                    % kron A kron I(b) is replicated 'ncol' (number of
                    % columns of x) times) along the diagonal.
                    if ncol>1
                        a = a*ncol;
                    end
                    
                    % Size of I(b)
                    b = 1;
                    for k = (index+1):length(opList)
                        if i > find(perm==k)
                            b = b * size(opList{k},2);
                        else
                            b = b * size(opList{k},1);
                        end
                    end
                    
                    % Size of the operator A=opList{index} to apply
                    r = sizes(index,2);
                    c = sizes(index,1);
                    
                    % (I(a) kron A kron I(b)) * x;
                    t = reshape(reshape(x,b,a*c).',c,a*b);
                    x = reshape(applyMultiply(opList{index},t,2)',a,r*b)';
                end
                y = reshape(x,n,ncol);
            end
        end % Multiply
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Divide
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = divide(op,b,mode)
            % Depends on sweepflag
            if op.sweepflag
                x = matldivide(op,b,mode);
            else
                x = lsqrdivide(op,b,mode);
            end
        end % divide
        
    end %Protected methods
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Private methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = private)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % best_permutation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Returns the best permutation associated to this Kronecker product
        function perm = best_permutation(op)
            list = op.children; % List of 'op''s children
            cost = zeros(length(list),2); % Computational costs of the
            % operators (children of 'op'). This is simply a numeric
            % representation of theirs shapes, which will affect 
            % computation time. Operators with low computational costs 
            % should be applied first.
            for i=1:length(list)
                % Cost = (nbr_rows-nbr_columns) / (size of the operator)
                cost(i,1) = (size(list{i},1)-size(list{i},2))/...
                    (size(list{i},1)*size(list{i},2));
                cost(i,2)=int8(i);
            end
            cost=sortrows(cost)';
            
            perm = cost(2,:);
        end
        
    end % private methods
end % opKron_1
