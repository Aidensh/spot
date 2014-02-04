classdef opWavelet_1 < opWavelet2_1
   %opWavelet_1   Wavelet operator.
   %
   %   opWavelet_1(N) creates a Wavelet transform for 1-dimensional signals
   %   of size N.  The wavelet transformation is computed using the Rice
   %   Wavelet Toolbox.
   %
   %   opWavelet_1(N,FAMILY) additionally specifies the FAMILY for the
   %   wavelet. Supported values for FAMILY are 'Daubechies' and 'Haar'.
   %
   %   opWavelet_1({C,dim1}) creates a Wavelet transform for 1-dimensional signals
   %   of size N.  The wavelet transformation is computed using the Rice
   %   Wavelet Toolbox.
   %
   %   opWavelet_1({C,dim1},FAMILY) additionally specifies the FAMILY for the
   %   wavelet. Supported values for FAMILY are 'Daubechies' and 'Haar'.
   %
   %   opWavelet_1(N,FAMILY,FILTER,LEVELS,REDUNDANT,TYPE) allows for four
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
   %   The opWavelet_1 operator is linear but not orthogonal. Therefore, the
   %   transpose of the operator is not the inverse operator. However, the
   %   inverse of the operator can be obtained through a left-inverse
   %   operation. For example:
   %               W = opWavelet_1(...)
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
        % opWavelet_1. Constructor.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function op = opWavelet_1(n, varargin)
            if ~iscell(n) && isscalar(n) && ~isa(n,'SeisDataContainer') % opWavelet_1(N,...)
                theN = n;
            elseif SDCpckg.utils.isForContainerInfo(n) % opWavelet_1({C,dim1}, ...)
                C = n{1};
                nInd = n{2};
                if numel(nInd) > 1
                    error('invalid dataContainer refrerncing dimension info'); % please make sure exactly 1 dimension of C is selected
                end
                theN = size(C,nInd);
                
            else
                error('Invalid input type of the first arguement for opWavelet_1');
            end
            op = op@opWavelet2_1(theN,1,varargin{:});
            op.type      = 'Wavelet_1';
            op.sweepflag = true;

        end % function opWavelet_1

   end % methods - public

end % opWavelet_1