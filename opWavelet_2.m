classdef opWavelet_2 < opWavelet2_2
   %opWavelet_2   Wavelet operator.
   %
   %   opWavelet_2(N) creates a Wavelet transform for 1-dimensional signals
   %   of size N.  The wavelet transformation is computed using the Rice
   %   Wavelet Toolbox.
   %
   %   opWavelet_2(N,FAMILY) additionally specifies the FAMILY for the
   %   wavelet. Supported values for FAMILY are 'Daubechies' and 'Haar'.
   %
   %   opWavelet_2({C,dim1}) creates a Wavelet transform for 1-dimensional signals
   %   of size N.  The wavelet transformation is computed using the Rice
   %   Wavelet Toolbox.
   %
   %   opWavelet_2({C,dim1},FAMILY) additionally specifies the FAMILY for the
   %   wavelet. Supported values for FAMILY are 'Daubechies' and 'Haar'.
   %
   %   opWavelet_2(N,FAMILY,FILTER,LEVELS,REDUNDANT,TYPE) allows for four
   %   additional parameters: FILTER (default 8) specifies the filter
   %   length, which must be even. LEVELS (default 5) gives the number of
   %   levels in the transformation. P does not need to be divisible by
   %   2^LEVELS. However, if LEVELS is bigger than LOG2(P), then LEVELS is
   %   adjusted to be equal to FLOOR(LOG2(P)). The Boolean field REDUNDANT
   %   (default false) indicates whether the wavelet is redundant. TYPE
   %   (default 'min') indictates what type of solution is desired; 'min'
   %   for minimum phase, 'max' for maximum phase, and 'mid' for mid-phase
   %   solutions.
   %
   %   The opWavelet_2 operator is linear but not orthogonal. Therefore, the
   %   transpose of the operator is not the inverse operator. However, the
   %   inverse of the operator can be obtained through a left-inverse
   %   operation. For example:
   %               W = opWavelet_2(...)
   %               y = W*x
   %               if z = W'*y, then z ~= x
   %               but, u = W\y, then u = x
   
   %   Copyright 2007-2009, Rayan Saab, Ewout van den Berg and Michael P. Friedlander
   %
   %   June  6, 2012: Added mirror symmetric extension of signals that are not
   %                  integer multiples of 2^levels.
   %                  Hassan Mansour (hassanm@cs.ubc.ca)
   %   June 25, 2012: Overloaded mldivide function to compute the inverse of
   %                  the operator.
   %                  Hassan Mansour (hassanm@cs.ubc.ca)
   %
   %   See the file COPYING.txt for full copyright information.
   %   Use the command 'spot.gpl' to locate this file.
   
   %   http://www.cs.ubc.ca/labs/scl/spot
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Methods - public
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % opWavelet_2. Constructor.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = opWavelet_2(varargin)
            p = 1;
            activated = true;
            if isempty(varargin)
                n = nan;
                p = nan;
                activated = false;
            elseif ~iscell(varargin{1}) && isscalar(varargin{1}) && ~isa(varargin{1},'SeisDataContainer') % opWavelet_2(N,...)
                n = varargin{1};
                varargin(1) = [];
            elseif SDCpckg.utils.isForContainerInfo(varargin{1}) && ~isempty(varargin{1})% opWavelet_2({C,dim1}, ...)
                C = varargin{1}{1};
                nInd = varargin{1}{2};
                if numel(nInd) > 1
                    error('invalid dataContainer refrerncing dimension info'); % please make sure exactly 1 dimension of C is selected
                end
                n = size(C,nInd);
                varargin(1) = [];
            elseif iscell(varargin{1}) && isempty(varargin{1})% opWavelet_2({}, ...)
                n = nan;
                p = nan;
                varargin(1) = [];
                activated = false;
            elseif ischar(varargin{1}) %opWavelet_2(FAMILY,FILTER,LEVELS,REDUNDANT,TYPE)
                % dynamic operator
                n = nan;
                p = nan;
                activated = false;
            else
                error('Invalid input type of the first arguement for opWavelet_2');
            end
            op = op@opWavelet2_2(n,p,varargin{:});
            op.type      = 'Wavelet_2';
            op.sweepflag = true;
            op.activated = activated;
        end % function opWavelet_2
        
        function op = activateOp(op,header,mode)
            % activate the operator right before it operates on container
            % header : the header of the data container
            
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
                p = prevHis.opN;
            else
                
                theSize = header.size;
                if any(strcmp('exsize',fieldnames(header)))
                    p = prod(theSize(header.exsize(1,1):header.exsize(2,1)));
                else
                    warning('exsize field not exist in header')
                    p = prod(theSize); % assume container is vectorized
                end
            end
            q = 1;
            
            isRedundant = op.redundant;
            theLevels = op.levels;
            % parts missing before calling super-class constructor
            if isRedundant
                if p == 1 || q == 1
                    theNseg =   theLevels + 1;
                else
                    theNseg = 3*theLevels + 1;
                end
                n = p*q;
                
                % find coefficient dimensions
                [pext, qext, theLevels] = CoeffDims(p, q, theLevels);
                
                p = pext*qext*theNseg;
            else
                theNseg = [];
                n    = p*q;
                
                % find coefficient dimensions
                [pext, qext, theLevels] = CoeffDims(p, q, theLevels);
                
                m = pext*qext;
            end
            op.levels = theLevels;
            % parts missing in during super-class constructor
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
            
            % parts missing after super-class constructor
            op.signal_dims = [p, q];
            op.coeff_dims = [pext, qext];
            op.redundant = isRedundant;
            op.nseg = theNseg;
            
            
            op.activated = true;
        end
        function dim = takesDim(~,~)
            % returns the number of dimension the operator operates on the
            % dataContainer.
            dim = 1;
        end
        
        
   end % methods - public

end % opWavelet_2

function [pext, qext, levels] = CoeffDims(p, q, levels)
            
            if p >= 2^levels
                plevels = levels;
                if q >= 2^levels
                    qext = ceil(q/(2^levels))*2^levels;
                elseif q > 1
                    qlevels = floor(log2(q));
                    levels = min(plevels,qlevels);
                    qext = ceil(q/(2^levels))*2^levels;
                else
                    qext = q;
                end
                pext = ceil(p/(2^levels))*2^levels;
            elseif p > 1
                plevels = floor(log2(p));
                if q >= 2^levels
                    levels = min(levels,plevels);
                    qext = ceil(q/(2^levels))*2^levels;
                elseif q > 1
                    qlevels = floor(log2(q));
                    levels = min(plevels,qlevels);
                    qext = ceil(q/(2^levels))*2^levels;
                else
                    levels = min(levels,plevels);
                    qext = q;
                end
                pext = ceil(p/(2^levels))*2^levels;
            else
                pext = p;
                if q >= 2^levels
                    qext = ceil(q/(2^levels))*2^levels;
                elseif q > 1
                    levels = floor(log2(q));
                    qext = ceil(q/(2^levels))*2^levels;
                else
                    qext = q;
                end
            end
end