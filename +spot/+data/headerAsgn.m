function h = headerAsgn(h,header,index)
%HEADERASGN  "SubsAsgn" function for SeisDataContainer headers
%
%   h = headerRef(header,index) basically returns a header with a subsrefed
%   range (in the dimensional sense) of its header entities

% Index checking
assert(isnumeric(index), 'index must be numeric!');
assert(isvector(index), 'index must be a vector!');

% Check index range
assert(index(1) >= 1, 'index out of bounds');
% assert(length(index) == length(header.origin),...
%    'header chunk must be same length as index');

% Remove singleton dimension
size = header.size;
if length(size) == length(header.origin) + 1 && size(end) == 1 
    size = size(1:end-1);
end

% Fill in the new ranges
h.size   = [h.size(1:index(1)-1) size h.size(min(length(h.size),index(end)+1):max(0,length(h.size)-length(size)+1))];
h.origin = [h.origin(1:index(1)-1) header.origin h.origin(min(length(h.origin),index(end)+1):max(0,length(h.origin)-length(header.origin) + 1))];
h.delta  = [h.delta(1:index(1)-1) header.delta h.delta(min(length(h.delta),index(end)+1):max(0,length(h.delta))-length(header.delta)+1)];
h.unit   = [h.unit(1:index(1)-1) header.unit h.unit(min(length(h.unit),index(end)+1):max(0,length(h.unit)-length(header.unit)+1))];
h.label  = [h.label(1:index(1)-1) header.label h.label(min(length(h.label),index(end)+1):max(length(h.label)-length(header.label)+1))];
h.dims = length(h.size);