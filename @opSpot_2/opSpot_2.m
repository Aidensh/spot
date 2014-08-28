classdef (HandleCompatible) opSpot_2
    %opSpot_2  Spot operator super class.
    %
    %   A = opSpot_2  creates an empty Spot operator.
    %
    %   A = opSpot_2(type,m,n)  creates a Spot operator named TYPE, of size
    %   M-by-N. CFLAG is set when the operator is
    %   complex. The TYPE and DATA fields provide the type of the operator
    %   (string) and additional data for printing.
    
    %   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
    %   See the file COPYING.txt for full copyright information.
    %   Use the command 'spot.gpl' to locate this file.
    
    %   http://www.cs.ubc.ca/labs/scl/spot
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties( SetAccess = protected )
        linear     = 1;     % Flag the op. as linear (1) or nonlinear (0)
        counter
        m          = 0;     % No. of rows
        n          = 0;     % No. of columns
        ms         = {};    % Vector of implicit rows
        ns         = {};    % Vector of implicit cols
        type       = '';
        cflag      = false; % Complexity of underlying operator
        children   = {};    % Constituent operators (for a meta operator)
        precedence = 1;
        sweepflag  = false; % whether we can do a sweep multiply, A*B
        isDirac    = false; % Whether we can skip this operator
        weights;            % weights for meta operators
        
        activated  = true;  % whether the operator is initialized. Assume true
        
        %metadata handling rule class. Default handler does nothing
        metadataRule;
    end
    
    properties( Dependent = true, SetAccess = private )
        nprods
    end
    
    properties(SetAccess = private, GetAccess = public)
        ID
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Public methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        
        function op = opSpot_2(type,m,n)
            %opSpot_2  Constructor.
            if nargin == 0
                % Relax -- empty constructor.
            elseif nargin == 3
                if ~(isnan(m) && isnan(n))
                    m = max(0,m);
                    n = max(0,n);
                    if round(m) ~= m || round(n) ~= n
                        warning('SPOT:ambiguousParams',...
                            'Size parameters are not integer.');
                        m = floor(m);
                        n = floor(n);
                    end
                    op.type = type;
                    op.m    = m;
                    op.n    = n;
                    op.ms   = {m};
                    op.ns   = {n};
                    op.counter = spot.counter();
                elseif isnan(m) && isnan(n) % non-activated operator
                    op.activated = false;
                    op.type = type;
                    op.counter = spot.counter();
                else
                    error('Something wrong with m,n in opSpot constructor.')
                end
            else
                error('Unsupported use of Spot constructor.');
            end
            op.ID = num2str(rand);
            op.metadataRule = sRule(strcat('op',type));
        end % function opSpot_2
        
        function nprods = get.nprods(op)
            %get.nprods  Get a count of the products with the operator.
            nprods = [op.counter.mode1, op.counter.mode2];
        end % function get.Nprods
        
    end % methods - public
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Protected methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = protected )
        
        function y = applyMultiply(op,x,mode)
            if size(x,2) > 1
%                 warning('opSpot_2:applyMultiply',['The sweeping function in ',...
%                     'applyMultiply is being ',...
%                 'discontinued due to performance issues. Please put all the ',...
%                 'multivector support in your operators multiply function and recode',...
%                 'accordingly']);
            end
        
            op.counter.plus1(mode);
            
            % For border case: empty x
            if isempty(x)
                if mode == 1
                    y = zeros(op.m,0,class(x));
                else
                    y = zeros(op.n,0,class(x));
                end
                return
            end
            
            if op.sweepflag
                y = op.multiply(x,mode);
            else
                x_n = size(x,2);
                
                if x_n == 1
                    y = op.multiply(x,mode);
                else
                    for i=x_n:-1:1
                        y(:,i) = op.multiply(x(:,i),mode);
                    end
                end
            end
        end % applyMultiply
        
        function y = applyDivide(op,x,mode)
            y = op.divide(x,mode);
        end
        
        % Signature of external protected functions (In class folder)
        y = divide(op,x,mode);
    end % methods - protected
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Abstract methods -- must be implemented by subclass.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Abstract, Access = protected )
        y = multiply(op,x,mode)
    end % methods - abstract
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - Static
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Static )
%         function fh = default_varName_handler()
%             fh = @(header,mode)header.varName;
%         end
%         function fh = default_varUnits_handler()
%             fh = @(header,mode)header.varUnits;
%         end
%         function fh = default_origin_handler()
%             fh = @(header,mode)header.origin;
%         end
%         function fh = default_delta_handler()
%             fh = @(header,mode)header.delta;
%         end
%         function fh = default_label_handler()
%             fh = @(header,mode)header.label;
%         end
%         function fh = default_unit_handler()
%             fh = @(header,mode)header.unit;
%         end
    end
end % classdef
