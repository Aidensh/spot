classdef opDFT < opOrthogonal
%OPDFT  Fast Fourier transform (DFT).
%
%   opDFT(M) create a unitary one-dimensional discrete Fourier
%   transform (DFT) for vectors of length M.
%
%   opDFT(M,CENTERED), with the CENTERED flag set to true, creates a
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
        function op = opDFT(m,centered)
            if nargin < 1 || nargin > 2
                error('Invalid number of arguments to opDFT.');
            end
            
            centered = nargin == 2 && centered;
            
            if  ~isscalar(m) || m~=round(m) || m <= 0
                error(['First argument to opDFT has to be a '...
                       'positive integer.']);
            end

            op = op@opOrthogonal('DFT',m,m);
            op.centered  = centered;
            op.cflag     = true;
            op.sweepflag = true;

            % Create function handle
            if centered
                op.funHandle = @opDFT_centered_intrnl;
            else
                op.funHandle = @opDFT_intrnl;
            end
        end % constructor
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
        function x = divide(op,x,mode)
            % Sweepable
            x = matldivide(op,x,mode);
        end % divide
    end % protected methods

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods - private
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods( Access = private )
        function x = opDFT_intrnl(op,x,mode)
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
            x = y;
        end % opDFT_intrnl

        function x = opDFT_centered_intrnl(op,x,mode)
            % One-dimensional DFT - Centered
            n = op.n;
            if mode == 1
                y = fftshift(fft(full(x)));
                y = y / sqrt(n);
            else
                y = ifft(ifftshift(full(x)));
                y = y * sqrt(n);
            end
            x = y;
        end % opDFT_centered_intrnl
    end % private methods
end % opDFT