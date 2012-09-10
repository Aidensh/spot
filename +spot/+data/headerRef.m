function h = headerRef(header,index)
%HEADERREF  "Subsref" function for SeisDataContainer headers
%
%   h = headerRef(header,index) basically returns a header with a subsrefed
%   range of its header entities

assert(isnumeric(index), 'index must be a numeric vector!');

dims   = header.dims;
size   = header.size;
origin = header.origin;
delta  = header.delta;
unit   = header.unit;
label  = header.label;