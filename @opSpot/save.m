function save( obj,fileName,overwrite,theFieldNames )
%SAVE Summary of this function goes here
%   Detailed explanation goes here

if(nargin == 2)
    overwrite = 0;
end

if (overwrite == 1)
    if (exist(fileName))
        delete(fileName);
    end
else
    assert(~exist(fileName),'file already exists')
end

save(fileName,'obj')
% %save(fileName,'-struct','obj');
% if nargin < 4
%     theFieldNames = fieldnames(obj);
% end
% 
% for ind = 1:length(theFieldNames)
%     if ind == 1
%         eval(strcat(theFieldNames{ind},'=obj.',theFieldNames{ind},';'));
%         save(fileName,theFieldNames{ind});
%     else
%         eval(strcat(theFieldNames{ind},'=obj.',theFieldNames{ind},';'));
%         save(fileName,'-append',theFieldNames{ind});
%     end
% end
end

