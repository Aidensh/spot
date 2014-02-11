classdef opCurvelet_2 < opSpot_2
%opCurvelet_2  Two-dimensional curvelet operator.
%
%   opCurvelet_2(M,N,NBSCALES,NBANGLES,FINEST,TTYPE,IS_REAL) creates a
%   two-dimensional curvelet operator for M by N matrices. The curvelet
%   transform is computed using the Curvelab code.
%
%   opCurvelet_2({C,dim1,dim2},NBSCALES,NBANGLES,FINEST,TTYPE,IS_REAL) creates a
%   two-dimensional curvelet operator for M by N matrices. The curvelet
%   transform is computed using the Curvelab code.
%
%   The remaining five parameters are optional; NBSCALES gives the number
%   of scales and is set to max(1,ceil(log2(min(M,N)) - 3)) by default, as
%   suggested by Curvelab. NBANGLES gives the number of angles at the
%   second coarsest level which must be a multiple of four with a minimum
%   of 8. By default NBANGLES is set to 16. FINEST sets whether to include
%   the finest scale of coefficients and is set to 0 by default; set this
%   to 1 to include the finest scale, or to 2 to keep the finest scale but
%   set it to zeros. TTYPE determines the type of transformation; either
%   'WRAP' for a wrapping transform or 'ME' for a Mirror-Extended
%   Transform, it's set to 'WRAP' by default.  IS_REAL sets whether the
%   transform is for real data or complex data.
%
%   See also CURVELAB.

%   An example that actually works is opCurvelet_2(32,32), use this for
%   testing.

%   Nameet Kumar - Oct 2010
%   Copyright 2009, Gilles Hennenfent, Ewout van den Berg and 
%   Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = protected)
        nbscales;
        nbangles;
        finest;   
        header;          %sizes of coefficient vectors
        nbcoeffs;           %total number of coefficients
        dims;           %size of curvelet
        ttype;           %type of transformation
    end % Properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = opCurvelet_2(varargin)
%         function op = opCurvelet_2(m,n,nbscales,nbangles,finest,ttype,...
%                 is_real)
            activated = true;
            if length(varargin) < 1 % if no input for constructor
                m = nan;
                n = nan;
                cn = nan;
                activated = false; % needs activation
            elseif SDCpckg.utils.isForContainerInfo(varargin{1}) %opCurvelet_2({C,dim1,dim2}, ...)
                containerDimInfo = varargin{1};
                
                C = containerDimInfo{1};
                
                if length(containerDimInfo)==1 || length(containerDimInfo{2}) > 2
                   error('invalid dataContainer refrerncing dimension info'); 
                end
                mInd = containerDimInfo{2}(1);
                nInd = containerDimInfo{2}(2);
                
                m = size(C,mInd);
                n = size(C,nInd);
                
                varargin(1) = [];
            elseif isscalar(varargin{1}) && isscalar(varargin{2})
                m = varargin{1};
                n = varargin{2};
                
                varargin(1) = []; % remove m
                varargin(1) = []; % remove n
            elseif iscell(varargin{1}) && isempty(varargin{1}) % if needs activation
                m = nan;
                n = nan;
                cn = nan;
                varargin(1) = [];
                activated = false;
            elseif isscalar(varargin{1}) %opCurvelet_2(NBSCALES,NBANGLES,FINEST,TTYPE,IS_REAL)
                m = nan;
                n = nan;
                cn = nan;
                activated = false; % needs activation
            else
                error('invalid first or second argument')
            end
            
            % the optional arguments
            if isempty(varargin) && ~(isnan(m) || isnan(n))
                theNbscales = max(1,ceil(log2(min(m,n)) - 3));
            elseif isempty(varargin) && (isnan(m) || isnan(n))
                theNbscales = nan;
            else % if ~isempty(varargin)
                theNbscales = varargin{1};
                varargin(1) = []; % remove nbscales
            end
            
            if isempty(varargin)
                theNbangles = 16;
            else
                theNbangles = varargin{1};
                varargin(1) = []; % remove nbangles
            end
            
            if isempty(varargin)
                theFinest = 0;
            else
                theFinest = varargin{1};
                varargin(1) = []; % remove finest
            end
            
            if isempty(varargin)
                theTtype = 'WRAP';
            else
                theTtype = varargin{1};
                varargin(1) = []; % remove nbscales
            end
            
            if isempty(varargin)
                is_real = 1;
            else
                is_real = varargin{1};
            end
            
            assert( strcmp(theTtype,'WRAP') || strcmp(theTtype,'ME'),...
                ['Please ensure ttype is set correctly. Options are'...
                ' "WRAP" for a wrapping transform and "ME" for a'...
                ' mirror-extended transform']);
            assert((any(theFinest == [0 1 2])) && (is_real==0||is_real==1),...
                'Please ensure finest and is_real are appropriate values');
            if activated
                assert( isscalar(theNbscales) && isscalar(theNbangles),...
                'Please ensure nbscales and nbangles are scalar values');
                
                if theFinest==0, assert( theNbscales>1, ['Please ensure that '...
                        'm and n are large enough for nbscales to be '...
                        'greater than 1 while finest is set to 0']);
                end
            end
            % Compute length of curvelet coefficient vector
            if activated
                if strcmp(theTtype,'ME')
                    C = mefcv2(randn(m,n),m,n,theNbscales,theNbangles);
                    
                    hdr{1}{1} = size(C{1}{1});
                    cn = prod(hdr{1}{1});
                    for i = 2:theNbscales
                        nw = length(C{i});
                        hdr{i}{1} = size(C{i}{1});
                        hdr{i}{2} = size(C{i}{nw/2+1});
                        cn = cn + nw/2*prod(hdr{i}{1}) + nw/2*prod(hdr{i}{2});
                    end
                else
                    [tmphdr, cn] = fdct_sizes_mex(m,n,theNbscales,theNbangles,...
                        logical(theFinest));
                    hdr = cell(1,theNbscales);
                    hdr{1} = {[tmphdr{1:2}]};
                    for i = 2:theNbscales - (~theFinest)
                        j = 3 + 5*(i-2);
                        hdr{i}={[tmphdr{j+1:j+2}];[tmphdr{j+3:j+4}];...
                            [tmphdr{j}]};
                    end
                    if ~theFinest,  hdr{end} = {[tmphdr{end-1:end}];1}; end;
                end
            end
            % Construct operator
            op = op@opSpot_2('Curvelet_2', cn, m*n);
            op.cflag     = ~is_real;
            op.nbscales  = theNbscales;
            op.nbangles  = theNbangles;
            op.finest    = theFinest;
            op.ttype     = theTtype;
            op.ns        = {[m n]};
            op.sweepflag = true;
            if activated
                op.header    = hdr;
                op.nbcoeffs  = cn;
                op.dims      = [m,n];
            end
        end % Constructor

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % rrandn             
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overloaded to produce a vector that really falls in the range of
        % op
        function y = rrandn(op)
            y = op.drandn;
            y = multiply(op,y,1);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % headerMod             
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overloaded to modify metadata correctly
        function h = headerMod(op,header,mode)
            exsize = header.exsize;

            if mode == 1
                h = header; % Copy header
                % Replace old first (collapsed) dimensional sizes with 
                % operator sizes.
                h.size(exsize(1,1):exsize(2,1)) = [];
                h.size   = [op.ms{:} h.size];
                h.origin(exsize(1,1):exsize(2,1)) = [];
                h.origin = [0 h.origin];
                h.delta(exsize(1,1):exsize(2,1)) = [];
                h.delta  = [1 h.delta];
                h.label(exsize(1,1):exsize(2,1)) = [];
                h.label  = ['lcurvelet' h.label];
                h.unit(exsize(1,1):exsize(2,1)) = [];
                h.unit   = ['ucurvelet' h.unit];

            else % mode == 2
                h = header;
                h.size(exsize(1,1):exsize(2,1)) = [];
                h.size   = [op.ns{:} h.size];
                h.origin(exsize(1,1):exsize(2,1)) = [];
                h.origin = [0 0 h.origin];
                h.delta(exsize(1,1):exsize(2,1)) = [];
                h.delta  = [1 1 h.delta];
                h.label(exsize(1,1):exsize(2,1)) = [];
                h.label  = ['l1' 'l2' h.label];
                h.unit(exsize(1,1):exsize(2,1)) = [];
                h.unit   = ['u1' 'u2' h.unit];
            end

            % Re-append correct exsize
            exsize_out = 1:length(h.size);
            exsize_out = [exsize_out;exsize_out];
            h.exsize   = exsize_out;
            
            
            % modify unit, origin, offset
            tempContainer = iCon(header);
            
            [unitOut, ~] = opCurvelet_2.unitGet(tempContainer,mode);
            h.unit = unitOut;
            h.varUnits = unitOut{1};
            
            [deltaOut,~] = opCurvelet_2.deltaGet(tempContainer,mode);
            h.delta = deltaOut;
            
            [offsetOut,~] = opCurvelet_2.originGet(tempContainer,mode);
            h.origin = offsetOut;
            % endOf : modify unit, origin, offset
            
            % history
            if ~isempty(h.IDHistory) % if the operation history is not empty
                prevHis = peek(h.IDHistory);
                
                prevHisID = prevHis.ID;
                
                if ~strcmp(prevHisID,op.ID) %if IDs are different
                    h.IDHistory = push(h.IDHistory,{...
                        {'ID',op.ID},...
                        {'ClassName',class(op)},...
                        {'mode',mode},...
                        {'origin',header.origin},...
                        {'opM',op.m},...
                        {'opN',op.n},...
                        {'opMs',op.ms},...
                        {'opNs',op.ns},...
                        {'children',{}}});
                else % if IDs are the same
                    prevHisMode = prevHis.mode;
                    
                    if strcmp(prevHisMode,mode) % if modes are the same
                        h.IDHistory = push(h.IDHistory,{...
                            {'ID',op.ID},...
                            {'ClassName',class(op)},...
                            {'mode',mode},...
                            {'origin',header.origin},...
                            {'opM',op.m},...
                            {'opN',op.n},...
                            {'opMs',op.ms},...
                            {'opNs',op.ns},...
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
                    {'origin',header.origin},...
                    {'opM',op.m},...
                    {'opN',op.n},...
                    {'opMs',op.ms},...
                    {'opNs',op.ns},...
                    {'children',{}}});
            end
            % endOf : history
            
        end % headerMod
        
        function op = activateOp(op,header,mode)
            % activate the operator right before it operates on container
            % header : the header of the data container
            
            % Note: there could be issues with activating opCurvelet.
            % Recommended to have activated when it was constructed.
            
            theSize = header.size;
            
            
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
                
                m = prevHis.opNs{:}(1);
                n = prevHis.opNs{:}(2);
            else % if not match
                if mode == 1
                m = theSize(1);
                n = theSize(2);
                else
                    error('Unable to recover the size for opCurvelet');
                end
            end
            
            
            
            % part missing before super-class constructor
            op.nbscales = max(1,ceil(log2(min(m,n)) - 3));
            
            assert( isscalar(op.nbscales) && isscalar(op.nbangles),...
                'Please ensure nbscales and nbangles are scalar values');
            
            if op.finest==0, assert( op.nbscales>1, ['Please ensure that '...
                    'm and n are large enough for nbscales to be '...
                    'greater than 1 while finest is set to 0']);
            end
            
            
            if strcmp(op.ttype,'ME')
                C = mefcv2(randn(m,n),m,n,op.nbscales,op.nbangles);
                
                hdr{1}{1} = size(C{1}{1});
                cn = prod(hdr{1}{1});
                for i = 2:op.nbscales
                    nw = length(C{i});
                    hdr{i}{1} = size(C{i}{1});
                    hdr{i}{2} = size(C{i}{nw/2+1});
                    cn = cn + nw/2*prod(hdr{i}{1}) + nw/2*prod(hdr{i}{2});
                end
            else
                [tmphdr, cn] = fdct_sizes_mex(m,n,op.nbscales,op.nbangles,...
                    logical(op.finest));
                hdr = cell(1,op.nbscales);
                hdr{1} = {[tmphdr{1:2}]};
                for i = 2:op.nbscales - (~op.finest)
                    j = 3 + 5*(i-2);
                    hdr{i}={[tmphdr{j+1:j+2}];[tmphdr{j+3:j+4}];...
                        [tmphdr{j}]};
                end
                if ~op.finest,  hdr{end} = {[tmphdr{end-1:end}];1}; end;
            end
            
            % part missing during super-class constructor
            theM = max(0,cn);
            theN = max(0,m*n);
            if round(theM) ~= theM || round(theN) ~= theN
                warning('SPOT:ambiguousParams',...
                    'Size parameters are not integer.');
                theM = floor(theM);
                theN = floor(theN);
            end
            op.m    = theM;
            op.n    = theN;
            op.ms   = {theM};
            op.ns   = {theN};
            
            % part missing after super-class constructor
            op.ns        = {[m n]};
            op.header    = hdr;
            op.nbcoeffs  = cn;
            op.dims      = [m,n];
            
            op.activated = true; % activation complete
        end
        
        function dim = takesDim(~,mode)
            % returns the number of dimension the operator operates on the
            % dataContainer.
            if mode == 1
                dim = 2;
            else
                dim = 1;
            end
        end
    end % Public Methods
       
 
    methods ( Access = protected )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = multiply(op,x,mode)
            x_n = size(x,2);
            
            if mode == 1
                for u = x_n:-1:1 % Loop over multivectors
                    x_tmp = x(:,u);
                    % Analysis mode
                    if strcmp(op.ttype,'ME')
                        x_tmp = mefcv2(reshape(x_tmp,op.dims(1),...
                                op.dims(2)),op.dims(1),op.dims(2),...
                                op.nbscales,op.nbangles);
                    else
                        x_tmp = fdct_wrapping_mex(op.dims(1),op.dims(2),...
                                op.nbscales,op.nbangles,...
                                logical(op.finest),...
                                reshape(x_tmp,op.dims(1),op.dims(2)));
                        if op.finest == 2, zero_finest_scale; end
                        if ~op.cflag % real transforms have redundancy 
                            x_tmp = fdct_wrapping_c2r(x_tmp);
                        end
                    end
                    y(:,u) = spot.utils.fdct_c2v(x_tmp,op.nbcoeffs);
                end
                x = y;
            else
                for u = x_n:-1:1 % Loop over multivectors
                    x_tmp = x(:,u);
                    % Synthesis mode  
                    if strcmp(op.ttype,'ME')
                        x_tmp = spot.utils.mefdct_v2c(full(x_tmp),...
                            op.header,op.nbangles);
                        x_tmp = meicv2(x_tmp,op.dims(1),op.dims(2),...
                            op.nbscales,op.nbangles);
                    else
                        x_tmp = spot.utils.fdct_v2c(x_tmp,op.header);
                        if op.finest == 2, zero_finest_scale; end
                        if ~op.cflag
                            x_tmp = fdct_wrapping_r2c(x_tmp);
                        end
                        x_tmp = ifdct_wrapping_mex(op.dims(1),...
                                op.dims(2),op.nbscales,op.nbangles,...
                                logical(op.finest),x_tmp);
                    end
                    if ~op.cflag % real transforms don't have complex
                                 % numbers for the inverse transform
                        x_tmp = real(x_tmp);
                    end
                    y(:,u) = x_tmp(:);
                end
                x = y;
            end

            %%% Nested Function
            function zero_finest_scale
                for i = 1:length(x_tmp{end})
                    x_tmp{end}{i} = zeros( size( x_tmp{end}{i} ) );
                end
            end

        end % Multiply
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Divide
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = divide(op,b,mode)
            % Non-sweepable
            x = lsqrdivide(op,b,mode);
        end % divide

    end % Protected Methods
    
    
    methods(Static)
        function [unitOut, unitIn] = unitGet(C,~)
            unitIn = unit(C);
            unitOut = unitIn; % same for now
        end
        
        function [originOut, originIn] = originGet(C,~)
            originIn = origin(C);
            originOut = zeros(1,length(origin(C))); % zero (for now?)
        end
        
        function [deltaOut, deltaIn] = deltaGet(C,~)
            deltaIn = delta(C);
            deltaOut = zeros(1,length(delta(C))); % zero (for now?)
        end
    end
end % opCurvelet_2