classdef opDFT_2 < opOrthogonal_2
%opDFT_2  Fast Fourier transform (DFT).
%
%   opDFT_2(M) create a unitary one-dimensional discrete Fourier
%   transform (DFT) for vectors of length M.
%
%   opDFT_2(M,CENTERED), with the CENTERED flag set to true, creates a
%   unitary DFT that shifts the zero-frequency component to the center
%   of the spectrum.
%
%   opDFT_2({C,dim1}) create an unitary one-dimensional discrete Fourier
%   transform (DFT) for vectors of length size(C).
%
%   opDFT_2({C,dim1},CENTERED), with the CENTERED flag set to true, creates a
%   unitary DFT that shifts the zero-frequency component to the center
%   of the spectrum.


%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = private )
        funHandle % Multiplication function
    end % private properties

    properties ( SetAccess = private, GetAccess = public )
        centered
    end % set-private get-public properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - public
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = opDFT_2(varargin)
            
            % work on and check the first input arguement
            if length(varargin) < 1 % if no input for constructor
                m = nan;
            elseif ~SDCpckg.utils.isForContainerInfo(varargin{1}) && isscalar(varargin{1}) % opDFT_2(M) or opDFT_2(M,CENTERED)
                if nargin > 2
                    error('Invalid number of arguments to opDFT_2 with first input double.');
                end
                m = varargin{1};
                if ~(m == round(m) && m > 0)
                    error('Input M has to be a positive integer')
                end
                varargin(1)=[];
            elseif SDCpckg.utils.isForContainerInfo(varargin{1}) % opDFT_2({C,dim1}) or opDFT_2({C,dim1},CENTERED)
                if nargin > 2
                   error('Invalid number of arguments to opDFT_2 with input SeisDataContainer.');
                end
                containerDimInfo = varargin{1}; % {C, dim1}
                
                C = varargin{1}{1}; % extract the DataContainer
                
                % now work with picking dimension. Note opDFT_2 allows
                % picking only one of the dimensions.
                if length(containerDimInfo) == 2 && numel(containerDimInfo{2}) == 1 % {C,dim1}
                    refDim = containerDimInfo{2};
                    if ~all(refDim <= length(size(containerDimInfo{1}))) % make sure not over select
                        error('invalid dimension selected for data container')
                    end
                    m = size(C,containerDimInfo{2});
                else
                    error('invalid dataContainer refrerncing dimension info'); % please make sure exactly 1 dimension of C is selected
                end
                varargin(1)=[];
            elseif isa(varargin{1},'SeisDataContainer')
                error('Remember to wrap the data container and all the dimension selection arguements together into one cell')
            elseif iscell(varargin{1}) && isempty(varargin{1})
                % centered will be handled in the next part
                m = nan;
                varargin(1) = [];
            else
                error('Invalid input type of the first arguement for opDFT_2');
            end
            
            
            % check centered arguement
            centred = false; % set centred to false by default
            if length(varargin) == 1 % if centered arguement is specified (remember varargin(1) is removed earlier)
                arg2 = varargin{1};
                if ~(islogical(arg2) || arg2 == 1 || arg2 == 0) % make sure the 2nd (centred) arguement makes sense. True or False or 1 or 0.
                    error(strcat('Invalid centered arguement of type "',class(arg2),'"'))
                end
                centred = true == arg2;
                
            elseif length(varargin) > 2
                error('opDFT_2 should not have more than 2 input arguments')
            end
            
            % create the operator
            op = op@opOrthogonal_2('DFT_2',m,m);
            op.centered  = centred;
            op.cflag     = true;
            op.sweepflag = true;
            % Create function handle
            if centred
                op.funHandle = @opDFT_2_centered_intrnl;
            else
                op.funHandle = @opDFT_2_intrnl;
            end
            
        end % constructor
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % headerMod
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            
            [unitOut, ~] = opDFT_2.unitGet(tempContainer,mode);
            h.unit = {unitOut};
            h.varUnits = unitOut;
            
            [deltaOut,~] = opDFT_2.deltaGet(tempContainer,mode);
            h.delta = deltaOut;
            
            [offsetOut,~] = opDFT_2.originGet(tempContainer,mode);
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
                        {'origin',header.origin},...
                        {'children',2}});
                else % if IDs are the same
                    prevHisMode = prevHis.mode;
                    
                    if strcmp(prevHisMode,mode) % if modes are the same
                        h.IDHistory = push(h.IDHistory,{...
                            {'ID',op.ID},...
                            {'ClassName',class(op)},...
                            {'mode',mode},...
                            {'origin',header.origin}
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
                    {'origin',header.origin},...
                    {'opM',op.m},...
                    {'opN',op.n},...
                    {'children',{}}});
            end
        end
        function op = activateOp(op,header,~)
            % activate the operator right before it operates on container
            % header : the header of the data container
            m = header.size;
            if any(strcmp('exsize',fieldnames(header)))
                m = prod(m(header.exsize(1,1):header.exsize(2,1)));
            else
                warning('exsize field not exist in header')
                m = prod(m); % assume container is vectorized
            end
            
            % from opSpot's constructor.
            % make up what the operator missed during it was initialized.
            m = max(0,m);
            if round(m) ~= m
                warning('SPOT:ambiguousParams',...
                    'Size parameters are not integer.');
                m = floor(m);
            end
            op.m    = m;
            op.n    = m;
            op.ms   = {m};
            op.ns   = {m};
            
            op.activated = true; % activation complete
        end
        function dim = takesDim(~,~)
            % returns the number of dimension the operator operates on the
            % dataContainer.
            dim = 1;
        end
    end % methods - public

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - protected
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = protected )
        % Multiplication
        function y = multiply(op,x,mode)
            y = op.funHandle(op,x,mode);
        end % Multiply

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Divide
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function x = divide(op,b,mode)
            % Sweepable
            x = matldivide(op,b,mode);
        end % divide
    end % protected methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - private
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = private )
        function y = opDFT_2_intrnl(op,x,mode)
            % One-dimensional DFT
            n = op.n;
            if mode == 1
                % Analysis
                y = fft(full(x));
                y = y / sqrt(n);
            else
                % Synthesis
                y = ifft(full(x));
                y = y * sqrt(n);
            end
        end % opDFT_2_intrnl

        function y = opDFT_2_centered_intrnl(op,x,mode)
            % One-dimensional DFT - Centered
            n = op.n;
            if mode == 1
                y = fftshift(fft(full(x)));
                y = y / sqrt(n);
            else
                y = ifft(ifftshift(full(x)));
                y = y * sqrt(n);
            end
        end % opDFT_2_centered_intrnl
        
        
    end % private methods
    
    methods ( Static )
        % Unit handling
        function [unitOut, unitIn] = unitGet(C,mode)
            % C is the SeisDataContainer
            % mode = 1 for regular fft; mode = 2 for inverse fft
            
            % unit : get
            theUnit = unit(C);
            
            % varUnit : get
            theVarUnit = varUnits(C);
            if ~(ischar(theVarUnit))
                error('something wrong with opDFT_2 varUnit')
            end
            
            % unit : check if it is recognized by opDFT_2
            if mode == 1
                unitRecognized = spLength.existUnit(theUnit) || spTime.existUnit(theUnit);
            else
                unitRecognized = spLengthInv.existUnit(theUnit) || spFrequency.existUnit(theUnit);
            end
            
            % varUnit : check if it is recognized by opDFT_2
            if mode == 1
                varUnitRecognized = spLengthInv.existUnit(theVarUnit) || spFrequency.existUnit(theVarUnit);
            else
                varUnitRecognized = spLength.existUnit(theVarUnit) || spTime.existUnit(theVarUnit);
            end
            
            if unitRecognized && varUnitRecognized
                if strcmp(theUnit, theVarUnit)
                    % good, do nothing
                else % if different
                    % will just use unit field
                    warning('SeisDataContainer:unit','unit and varUnit are both recognized but different, so assumed unit field as the unit')
                end
            elseif unitRecognized && ~varUnitRecognized
                % will just use unit field
            elseif ~unitRecognized && varUnitRecognized
                % recover unit from varUnit
                theUnit = theVarUnit;
            else % ~inUnitRecognized && ~outUnitRecognized
                % opDFT_2 isn't sure what it's doing. Maybe give warning for
                % this or blow the program up?
                warning('SeisDataContainer:unit','opDFT_2 detects no recognizable unit from Leave it as "(someUnit)".')
                theUnit = '(someUnit)';
            end
            
            unitIn = theUnit;
            
            if mode == 1
                if spLength.existUnit(unitIn)
                    unitOut = '1/m';
                elseif spTime.existUnit(unitIn)
                    unitOut = 'Hz';
                else
                    unitOut = '(someUnit)';
                end
            else
                if spLengthInv.existUnit(unitIn)
                    unitOut = 'm';
                elseif spFrequency.existUnit(unitIn)
                    unitOut = 's';
                else
                    unitOut = '(someUnit)';
                end
            end
            
            
        end
        
        function [originOut, originIn] = originGet(C,~)
            originIn = origin(C,1);
            originOut = zeros(1,length(origin(C,1)));
        end
        
        function [deltaOut, deltaIn] = deltaGet(C,~)
            deltaIn = delta(C,1);
            deltaOut = 1/((length(C) - 1) * delta(C,1));
        end
    end
end % opDFT_2