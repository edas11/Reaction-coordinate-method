function [ xNew ] = eqPointDistance( x, y, limit )
%EQPOINTDISTANCE Reduces so that all new x and corresponding y points would
%have more or less the same distance to each other.

xNew(1) = x(1);
distances = sqrt(diff(x).^2 + diff(y).^2);

culmDistance = 0;
for i=2:length(x)
    culmDistance = culmDistance + distances(i-1);
    if culmDistance<limit
        continue;
    else
        xNew(end+1) = x(i);
        culmDistance = 0;
    end
end

