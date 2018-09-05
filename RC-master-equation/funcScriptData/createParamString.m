function [ paramString ] = createParamString( names, params )
%CREATEPARAMSTRING Creates file name

paramString = '';
for i=1:length(names)
    paramString = strcat(paramString, names{i}, num2str(params(i)));
end

end

