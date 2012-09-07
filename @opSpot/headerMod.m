function h = headerMod(op,xmeta,header,mode)
%HEADERMOD  Header Modifying Function for SeisDataContainer compatibility
%
%   h = headerMod(op,xmeta,header,mode) returns a modified header 
%   corresponding to the actions the spot operator would have on the 
%   metadata of the SeisDataContainer.
%
%   op = Spot operator
%   xmeta = explicit metadata stored on the datacon
%   header = header of the

% By default this does nothing to the header whatsoever.
if mode == 1
    h = header;
else
    h = header;
end