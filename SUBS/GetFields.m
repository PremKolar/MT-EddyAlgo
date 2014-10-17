function [F,unreadable]=GetFields(file,keys)
    F=struct;
    unreadable.is=false;
    for field=fieldnames(keys)';ff=field{1};
        %%
        if isempty(ff),continue;end
        %%
        try
            F.(ff)=tryCase(ff,file,keys);
        catch uc
            unreadable=cathCase(uc);
            return
        end
    end
    %%
    if isfield(F,'lon')
        if numel(F.lon)==length(F.lon)
            [F.lon,F.lat]=meshgrid(F.lon,F.lat);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fff = tryCase(ff,file,keys)
    sqDouNCv=@(F,k,ff) squeeze(double(nc_varget(F,k.(ff))));
    Fff = sqDouNCv(file,keys,ff);
    if strcmpi(ff,'lon')
        Fff =  wrapTo360(Fff);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function unreadable=cathCase(uc)
    unreadable.is=true;
    unreadable.catch=uc;
    disp('skipping'); disp(uc);
    disp(uc.message);
    disp(uc.getReport);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%