function h = headerRef(header,index)
%HEADERREF  "Subsref" function for SeisDataContainer headers
%
%   h = headerRef(header,index) basically returns a header with a subsrefed
%   range (in the dimensional sense) of its header entities

% Index checking
index = index(1):index(end);
assert(isnumeric(index), 'index must be numeric!');
assert(isvector(index), 'index must be a vector!');
% singleton = false;

sizes   = header.size;
origin  = header.origin;
delta   = header.delta;
unit    = header.unit;
label   = header.label;
% if length(sizes) == 1, sizes = [sizes 1]; singleton = true; end

% Check index range
if ~(index(1) >= 1 && index(end) <= length(sizes))
    warning( 'index out of bounds');
end

% Update dims
% if singleton
%     dims = 1;
%     index = 1;
%     sizes = sizes(1);
% else
    dims = length(sizes(index));
% end

% Fill in the new ranges
h        = header;
h.dims   = dims;
h.size   = sizes(index);
h.origin = origin(index);
h.delta  = delta(index);
h.unit   = unit(index);
h.label  = label(index);