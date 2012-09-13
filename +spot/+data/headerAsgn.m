function h = headerAsgn(h,header,index)
%HEADERASGN  "SubsAsgn" function for SeisDataContainer headers
%
%   h = headerRef(header,index) basically returns a header with a subsrefed
%   range (in the dimensional sense) of its header entities

% Index checking
assert(isnumeric(index), 'index must be numeric!');
assert(isvector(index), 'index must be a vector!');

% Check and remove singleton dimension
assert(length(h.origin) > 1, 'dims have to be at least 2');

% Check index range
assert(index(1) >= 1 && index(end) <= length(h.size), 'index out of bounds');
% assert(length(index) == length(header.origin),...
%    'header chunk must be same length as index');

% Remove singleton dimension
size = header.size;
if length(size) == length(header.origin) + 1 && size(end) == 1 
    size = size(1:end-1);
end

% Fill in the new ranges
h.size   = [h.size(1:index-1) size h.size(index+1:end)];
h.origin = [h.origin(1:index-1) header.origin h.origin(index+1:end)];
h.delta  = [h.delta(1:index-1) header.delta h.delta(index+1:end)];
h.unit   = [h.unit(1:index-1) header.unit h.unit(index+1:end)];
h.label  = [h.label(1:index-1) header.label h.label(index+1:end)];