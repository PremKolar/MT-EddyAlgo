function [lon]=CorrectLongitude(lon)
    % longitude(-180:180) concept is to be used!
    if max(lon(:))>180
        lontrans=true;
    else
        lontrans=false;
    end
    if lontrans
        lon(lon>180)=lon(lon>180)-360;
    end
end
