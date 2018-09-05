function [ densityHEOM ] = importHEOM( folder, paramString  )
%IMPORTHEOM Load HEOM data
%   Columns must be seperated with \t zenklu. First column should be time
%   (fs) other - real and imaginary density operator parts. Elements should
%   be in this order: 11, 12, 21, 22.

    HEOMdata = importdata( strcat(folder, paramString,'.txt'),'\t' );
    HEOMdata = HEOMdata.data;

    dim = sqrt( (size(HEOMdata, 2) - 1)/2 );
    
    totalSteps = length( HEOMdata );

    % Determines which HEODdata elements are real and imaginary
    map = 2:2:2*dim^2;
    map = reshape(map, [dim dim]);
    map = map.';
    
    % Takes HEOM data
    densityHEOM = zeros(totalSteps, dim, dim);
    for row=1:dim
        for col=1:dim
            densityHEOM(:,row,col) = HEOMdata(:, map(row,col)) ...
                + 1i*HEOMdata(:, map(row,col) + 1 ) ;
        end
    end

end

