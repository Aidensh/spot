classdef opDCT2 < opOrthogonal
%OPDCT2  Two-dimensional discrete cosine transform (DCT).
%
%   opDCT2(M,N) creates a two-dimensional discrete cosine transform
%   operator for matrices of size M-by-N. Input and output of the
%   matrices is done in vectorized form. When N is omitted it is set
%   to M by default.

%   Copyright 2009, Ewout van den Berg and Michael P. Friedlander
%   See the file COPYING.txt for full copyright information.
%   Use the command 'spot.gpl' to locate this file.

%   http://www.cs.ubc.ca/labs/scl/spot
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Properties
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   properties( SetAccess = private )
      inputdims;    % Dimensions of the input
   end % properties - private
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Methods - public
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % opDCT2. Constructor.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function op = opDCT2(m,n)
         if nargin < 1
            error('Too few arguments');
         elseif nargin == 1
            n = m;
         end
         if ~isscalar(m) || m~=round(m) || m <= 0
            error('First argument to opDCT2 must be a positive integer.');
         end
         if ~isscalar(n) || n~=round(n) || n <= 0
            error('Second argument to opDCT2 must be a positive integer.');
         end
         op = op@opOrthogonal('DCT2',m*n,m*n);
         op.inputdims = [m,n];
          op.sweepflag  = true;
      end % function opDCT2
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % spy. Really only a pedagogical tool, and only practical to
      % execute for DCTs that have less than, say, a dozen columns.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function spy(op)
         
         colormap('gray');
         mSig = op.inputdims(1);
         nSig = op.inputdims(2);
         x    = zeros(mSig*nSig,1);
         
         % Create image with single pixel.
         k = 0;
         for i=1:mSig
            for j=1:nSig
               k    = k + 1;
               x(k) = 1;
               y    = op'*x;
               x(k) = 0;
               subplot(nSig,mSig,k);
               imagesc(reshape(y,mSig,nSig));
               axis square; axis off;
            end
         end
      end % function plot
      
   end % methods - public
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Methods - protected
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods( Access = protected )
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % multiply.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function y = multiply(op,x,mode)
         m = op.inputdims(1);
         n = op.inputdims(2);
         q = size(x,2);
         if mode == 1
             if isscalar(op)
                % special case: allocate result size of x
                y = zeros(size(x),class(x));
             else
                y = zeros(op.m,q, class(x));
             end
             
             for u = 1:q
                y_tmp  = spot.utils.dct(full(reshape(x(:,u),m,n)));
                y_tmp  = spot.utils.dct(y_tmp')';
                y(:,u) = y_tmp(:);
             end
         else
             if isscalar(op)
                % special case: allocate result size of x
                y = zeros(size(x),class(x));
             else
                y = zeros(op.n,q, class(x));
             end
             
             for u = 1:q
                y_tmp = spot.utils.idct(full(reshape(x(:,u),m,n)));
                y_tmp = spot.utils.idct(y_tmp')';
                y(:,u) = y_tmp(:);
             end
         end
      end % function multiply
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Divide
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function x = divide(op,b,mode)
          % Non-sweepable
          x = lsqrdivide(op,b,mode);
      end % divide
      
   end % methods - protected
   
end % classdef
