function OUT = PTranspose(Data, dim1, dim2, dist)

if numlabs ~= 1
    error('This function is not meant to be called from an spmd block')
end
if ~isa(Data, 'distributed')
    error('Input must be distributed')
end
if ~isnumeric(dim1) || ~isnumeric(dim2) || length(dim1) ~= 1 ...
        || length(dim2)~=1 || dim1 > ndims(Data) || dim2 > ndims(Data)
    error('Transpose dimensions must be individual dimensions of Data')
end
if ~isnumeric(dist) || length(dist) ~= 1 || dist > ndims(Data)
    error('Distribution dimension must be an individual dimension of Data')
end

dims = size(Data);
temp = dims(dim1);
dims(dim1) = dims(dim2);
dims(dim2) = temp;

perm = 1:ndims(Data);
perm(dim1) = dim2;
perm(dim2) = dim1;

spmd
    codistr = getCodistributor(Data);
    dimdist = codistr.Dimension;
    if find( [dim1,dim2] == dimdist ) 
        if find([dim1,dim2] == dist)
            if dist == dimdist
                dist = perm(dist);
                dimdist = dist;
            else
                temp = dist;
                dist = dimdist;
                dimdist = temp;
            end
        else
            if dist ~= dimdist             
                codistr = codistributor1d( dist);
                Data = redistribute( Data, codistr );
                dimdist = dist;
            end
        end
    end
    d = permute(getLocalPart(Data), perm);
    codistr = codistributor1d( dimdist, ...
        codistributor1d.unsetPartition, dims);
    Data = codistributed.build(d, codistr);
    if dimdist ~= dist
        dist = perm(dist);
        codistr = codistributor1d(dist);
        Data = redistribute( Data, codistr);
    end
end
OUT = Data;


end