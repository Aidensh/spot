function h = headerMod(op,header,mode)
%HEADERMOD  Header Modifying Function for SeisDataContainer compatibility
%
%   h = headerMod(op,header,mode) returns a modified header corresponding
%   to the actions the spot operator would have on the metadata of the
%   SeisDataContainer.

% By default this does nothing whatsoever.
h = header;