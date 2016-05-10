function [ velo ] = velo( mp,temp )

    velo = sqrt(8*1.381e-23*temp/pi./mp);
    
end

