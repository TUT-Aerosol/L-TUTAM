function out = rajat( lowerLimit,in,upperLimit )
%RAJAT limits 'in' between 'lowerLimit' and 'upperLimit'

out = min(max(lowerLimit,in),upperLimit);

end

