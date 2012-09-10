function h = headerRef(header,index)
%HEADERREF  "Subsref" function for SeisDataContainer headers
%
%   h = headerRef(header,index) basically returns a header with a subsrefed
%   range (in the dimensional sense) of its header entities

% Index checking
assert(isnumeric(index), 'index must be numeric!');
assert(isvector(index), 'index must be a vector!');

size    = header.size;
origin  = header.origin;
delta   = header.delta;
unit    = header.unit;
label   = header.label;

% Check and remove singleton dimension
assert(length(origin) > 1, 'dims have to be at least 2');

% Check index range
assert(index(1) >= 1 && index(end) < length(origin), 'index out of bounds');

% Update dims
dims = length(size(index));

% Fill in the new ranges
h        = header;
h.dims   = dims;
h.size   = size(index);
h.origin = origin(index);
h.delta  = delta(index);
h.unit   = unit(index);
h.label  = label(index);