function [LON]=CorrectLongitude(LON)
    % longitude(-180:180) concept is to be used!
    if max(LON(:))>180
        lontrans=true;
    else
        lontrans=false;
    end
    if lontrans
        LON(LON>180)=LON(LON>180)-360;
    end
end