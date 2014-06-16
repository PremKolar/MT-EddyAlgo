function [F]=GetFields(file,keys)
    F=struct;    
    for field=fieldnames(keys)';ff=field{1};
        if strcmpi(ff,'lon')
            F.(ff) = CorrectLongitude(squeeze(nc_varget(file,keys.(ff))));
        else
            F.(ff) = squeeze(nc_varget(file,keys.(ff)));
        end
    end    
    if min(size(F.lon))==1
        F.lon=repmat(standVectorUp(F.lon)',length(F.lat),1);
        F.lat=repmat(standVectorUp(F.lat) ,1,length(F.lon));
    end
end


