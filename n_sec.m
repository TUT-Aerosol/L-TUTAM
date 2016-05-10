function nSec = n_sec( malli )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
switch(malli(1:2))
    case 'PL'
        nSec=6;
    case 'LN'
        nSec=6;
    otherwise
        nSec=str2double(malli(3:end));
end

end

