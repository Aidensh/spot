classdef opWavelet2_1 < opSpot
    %OPWAVELET   Wavelet operator.
    %
    %   opWavelet(P,Q,FAMILY) creates a Wavelet operator of given FAMILY for
    %   signals of size P-by-1. The wavelet transformation is computed using
    %   the Rice Wavelet Toolbox.  The values supported for FAMILY are
    %   'Daubechies' and 'Haar'.
    %
    %   opWavelet(P,Q,FAMILY,FILTER,LEVELS,REDUNDANT,TYPE) allows for four
    %   additional parameters: FILTER (default 8) specifies the filter length,
    %   which must be even. LEVELS (default 5) gives the number of levels in
    %   the transformation. P and Q do not need to be divisible by 2^LEVELS.
    %   However, if LEVELS is bigger than LOG2(MIN(P,Q)), then LEVELS is
    %   adjusted to be equal to FLOOR(LOG2(MIN(P,Q))).
    %   The Boolean field REDUNDANT (default false) indicates whether the wavelet
    %   is redundant. TYPE (default 'min') indictates what type of solution is
    %   desired; 'min' for minimum phase, 'max' for maximum phase, and 'mid'
    %   for mid-phase solutions.
    %
    %   opWavelet({C,dim1,dim2},[optional inputs]). dim1 and dim2 are the
    %   index to the dimension working on. P would be the size of (dim1)-th
    %   dimension and Q would be the size of (dim2)-th dimension
    %
    %   The opWavelet operator is linear but not orthogonal. Therefore, the
    %   transpose of the operator is not the inverse operator. However, the
    %   inverse of the operator can be obtained through a left-inverse
    %   operation. For example:
    %               W = opWavelet(...)
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
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties( SetAccess = private, GetAccess = public )
        family     = 'Daubechies';   % Wavelet family
        lenFilter  = 8;              % Filter length
        filter                       % Filter computed by daubcqf
        levels     = 5;              % Number of levels
        typeFilter = 'min'
        redundant  = false;          % Redundant flag
        nseg
        signal_dims                  % Dimensions of the signal domain
        coeff_dims                   % Dimensions of extended coefficients
        funHandle                    % Multiplication function
        funHandle2                   % Divide function
    end % Properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - public
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % opWavelet. Constructor.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %function op = opWavelet2_1(p,q,family,lenFilter,levels,redundant,typeFilter)
        %function op = opWavelet2_1({C,dim1,dim2},family,lenFilter,levels,redundant,typeFilter)
        function op = opWavelet2_1(varargin)
            
            % first (and second) input
            if length(varargin) < 1 % if no input for constructor
                error('need some input for opDFT') % blow it up (for now)
            elseif ~SDCpckg.utils.isForContainerInfo(varargin{1}) && isscalar(varargin{1}) % opWavelet2_1(P,Q,[...])
                p = varargin{1};
                q = varargin{2};
                
                varargin(2) = [];
                varargin(1) = [];
            elseif SDCpckg.utils.isForContainerInfo(varargin{1}) % opWavelet2_1({C,dim1,dim2},[...])
                p = isize(varargin{1}{1},varargin{1}{2}(1));
                q = isize(varargin{1}{1},varargin{1}{2}(2));
                
                varargin(1) = [];
            elseif isa(varargin{1},'SeisDataContainer')
                error('Remember to wrap the data container and all the dimension selection arguements together into one cell')
            else
                error('Invalid input type of the first arguement for opWavelet2_1');
            end
            
            lenInput = length(varargin);
            
            if lenInput < 1 || isempty(varargin{1}) % originally : if nargin <= 2 || isempty(family)
                theFamily = 'Daubechies';
            else
                theFamily = varargin{1};
            end
            if lenInput < 3 || isempty(varargin{3}) % originally : if nargin < 5 || isempty(levels)
                theLevels = 5;
            else
                theLevels = varargin{3};
            end
            if lenInput >= 4 && varargin{4} % originally : if nargin >= 6 && redundant
                if p == 1 || q == 1
                    theNseg =   theLevels + 1;
                else
                    theNseg = 3*theLevels + 1;
                end
                n = p*q;

                % find coefficient dimensions
                [pext, qext, theLevels] = CoeffDims(p, q, theLevels);

                m = pext*qext*theNseg;

                theRedundant = true;
            else
                theNseg = [];
                n    = p*q;

                % find coefficient dimensions
                [pext, qext, theLevels] = CoeffDims(p, q, theLevels);

                m = pext*qext;

                theRedundant = false;
            end

            op = op@opSpot('Wavelet2_1', m, n);
            op.signal_dims = [p, q];
            op.coeff_dims = [pext, qext];
            op.levels = theLevels;
            op.redundant = theRedundant;
            op.nseg = theNseg;

            if lenInput >= 2 && ~isempty(varargin{2})
                op.lenFilter = varargin{2};
            end
            if lenInput >= 5 && ischar(varargin{5})
                op.typeFilter  = varargin{5};
            end
            switch lower(theFamily)
                case {'daubechies'}
                    op.family = 'Daubechies';
                    op.filter = spot.rwt.daubcqf(op.lenFilter,op.typeFilter);

                case {'haar'}
                    op.family = 'Haar';
                    op.filter = spot.rwt.daubcqf(0);

                otherwise
                    error('Wavelet family %s is unknown.', theFamily);
            end

            % Initialize function handle
            if theRedundant
                op.funHandle = @multiply_redundant_intrnl;
                op.funHandle2 = @divide_redundant_intrnl;
            else
                op.funHandle = @multiply_intrnl;
                op.funHandle2 = @divide_intrnl;
            end
            op.sweepflag = true;

        end % constructor
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % headerMod
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function h = headerMod(op,header,mode)
            
            %default headerMod
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
            % endOf : default deaderMod
            
            % modify unit, origin, offset
            tempContainer = iCon(header);
            
            [unitOut, ~] = opWavelet2_1.unitGet(tempContainer,mode);
            h.unit = unitOut;
            %h.varUnits = unitOut;
            
            [deltaOut,~] = opWavelet2_1.deltaGet(tempContainer,mode);
            h.delta = deltaOut;
            
            [originOut,~] = opWavelet2_1.originGet(tempContainer,mode);
            h.origin = originOut;
            % endOf : modify unit, origin, offset
            
            % operator history
            if ~isempty(h.IDHistory) % if the operation history is not empty
                prevHis = peek(h.IDHistory);
                
                prevHisID = SDCpckg.utils.getInfoFromCell(prevHis,'ID');
                
                if ~strcmp(prevHisID,op.ID) %if IDs are different
                    push(h.IDHistory,{...
                        {'ID',op.ID},...
                        {'mode',mode},...
                        {'origin',header.origin},...
                        {'delta',header.delta}});
                else % if IDs are the same
                    prevHisMode = SDCpckg.utils.getInfoFromCell(prevHis,'mode');
                    
                    if strcmp(prevHisMode,mode) % if modes are the same
                        push(h.IDHistory,{...
                            {'ID',op.ID},...
                            {'mode',mode},...
                            {'origin',header.origin},...
                            {'delta',header.delta}});
                    else % if modes are opposite
                        pop(h.IDHistory); % cancel out
                        
                        % info recovery
                        h.origin = SDCpckg.utils.getInfoFromCell(prevHis,'origin');
                        h.delta = SDCpckg.utils.getInfoFromCell(prevHis,'delta');
                    end
                end
            else % if it is empty ...
                push(h.IDHistory,{...
                    {'ID',op.ID},...
                    {'mode',mode},...
                    {'origin',header.origin},...
                    {'delta',header.delta}});
            end
            % endOf : operator history
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Divide
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = mldivide(op,x)
            if isa(x,'SeisDataContainer')
                y = dataDivide(x,op);
            else
                y = op.funHandle2(op,x);
            end
        end % function mldivide

    end % methods - public

    methods( Access = private )


      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % matvec.  Application of Wavlet operator.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function y = multiply_intrnl(op,x,mode)
         p = op.signal_dims(1);
         q = op.signal_dims(2);
         pext = op.coeff_dims(1);
         qext = op.coeff_dims(2);

         levels = op.levels; filter = op.filter;
         if issparse(x), x = full(x); end

         % apply matvec operation
         R = opExtend(p,q,pext,qext);

         if mode == 1

            % extend the signal
            xext = R*x;

            % reshape the extended signal
            Xmat = reshape(xext,pext,qext);

            if isreal(x)
               y  = spot.rwt.mdwt(Xmat, filter, levels);
            else
               y1 = spot.rwt.mdwt(real(Xmat), filter, levels);
               y2 = spot.rwt.mdwt(imag(Xmat), filter, levels);
               y  = y1 + 1i * y2;
            end
            y = y(:);
         else % mode == 2
            Xmat = reshape(x,pext,qext);
            if isreal(x)
               y = spot.rwt.midwt(Xmat, filter, levels);
            else
               y1 = spot.rwt.midwt(real(Xmat), filter, levels);
               y2 = spot.rwt.midwt(imag(Xmat), filter, levels);
               y  = y1 + 1i * y2;
            end

            % apply adjoint of extension operator
            y = R'*y(:);

         end
      end % function matvec

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % matvec_redundant.  Application of redundant Wavlet operator.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function y = multiply_redundant_intrnl(op,x,mode)
         p = op.signal_dims(1);
         q = op.signal_dims(2);
         pext = op.coeff_dims(1);
         qext = op.coeff_dims(2);

         nseg = op.nseg;
         levels = op.levels; filter = op.filter;
         if issparse(x), x = full(x); end

         R = opExtend(p,q,pext,qext);

         if mode == 1

            % extend the signal
            xext = R*x;

            % reshape the extended signal
            Xmat = reshape(xext,pext,qext);

            if isreal(x)
               [yl,yh] = spot.rwt.mrdwt(Xmat, filter, levels);
               y = [yl,yh];
            else
               [yl1,yh1] = spot.rwt.mrdwt(real(Xmat), filter, levels);
               [yl2,yh2] = spot.rwt.mrdwt(imag(Xmat), filter, levels);
               y = [yl1,yh1] + 1i * [yl2,yh2];
            end
            y = y(:);
         else % mode == 2
            xl = reshape(x(1:pext*qext),pext,qext);
            xh = reshape(x(pext*qext+1:end),pext,(nseg-1)*qext);
            if isreal(x)
               y = spot.rwt.mirdwt(xl, xh, filter, levels);
            else
               y1 = spot.rwt.mirdwt(real(xl), real(xh), filter, levels);
               y2 = spot.rwt.mirdwt(imag(xl), imag(xh), filter, levels);
               y = y1 + 1i * y2;
            end

            % apply adjoint of extension operator
            y = R'*y(:);

         end
      end % function matvec_redundant


      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % divide_intrnl.  Application of redundant Wavlet operator.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function y = divide_intrnl(op,x)
         p = op.signal_dims(1);
         q = op.signal_dims(2);
         pext = op.coeff_dims(1);
         qext = op.coeff_dims(2);

         levels = op.levels; filter = op.filter;
         if issparse(x), x = full(x); end


         Xmat = reshape(x,pext,qext);
         if isreal(x)
            y = spot.rwt.midwt(Xmat, filter, levels);
         else
            y1 = spot.rwt.midwt(real(Xmat), filter, levels);
            y2 = spot.rwt.midwt(imag(Xmat), filter, levels);
            y  = y1 + 1i * y2;
         end


         % clip signal back to original dimensions
         y = y(1:p, 1:q);

         y = y(:);
      end % function divide

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % divide_intrnl.  Application of redundant Wavlet operator.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function y = divide_redundant_intrnl(op,x)
         p = op.signal_dims(1);
         q = op.signal_dims(2);
         pext = op.coeff_dims(1);
         qext = op.coeff_dims(2);

         nseg = op.nseg;
         levels = op.levels; filter = op.filter;
         if issparse(x), x = full(x); end


         xl = reshape(x(1:pext*qext),pext,qext);
         xh = reshape(x(pext*qext+1:end),pext,(nseg-1)*qext);
         if isreal(x)
            y = spot.rwt.mirdwt(xl, xh, filter, levels);
         else
            y1 = spot.rwt.mirdwt(real(xl), real(xh), filter, levels);
            y2 = spot.rwt.mirdwt(imag(xl), imag(xh), filter, levels);
            y = y1 + 1i * y2;
         end


         % clip signal back to original dimensions
         y = y(1:p, 1:q);

         y = y(:);


      end % function divide



    end % methods - private

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - protected
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = protected )

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Multiply
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function y = multiply(op,x,mode)
            for u = size(x,2):-1:1
                y(:,u) = op.funHandle(op,x(:,u),mode);
            end
        end % function multiply

    end % methods - protected
    
    methods(Static)
        function [unitOut, unitIn] = unitGet(C,~)
            unitIn = unit(C);
            unitOut = unitIn; % same for now
        end
        
        function [originOut, originIn] = originGet(C,~)
            originIn = origin(C);
            originOut = 0; % zero (for now?)
        end
        
        function [deltaOut, deltaIn] = deltaGet(C,~)
            deltaIn = delta(C);
            deltaOut = 0; % zero (for now?)
        end
    end

end % classdef

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