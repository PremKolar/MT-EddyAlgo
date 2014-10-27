function avi2pop
    init_threads(12)
    ncfile.avi='../avidim.nc';
    ncfile.pop='../popdim.nc';
    path.in='/scratch/uni/ifmto/u241194/DAILY/EULERIAN/SSH/';
    path.out='/scratch/uni/ifmto/u300065/FINAL/POP2AVIssh/';
    files.in=dir([path.in 'SSH_GLB_t.t0.1_42l_CORE.*.nc']);
    %%
    [lat,lon]=inits(ncfile);
    %%
    [idx]=CrossRef2Closest(lon,lat);
    %%
    %     [Ya,Xa]=size(lon.avi);
    ref=makeXref(idx);
    %%
    remapSSH(ref,files,path,lon)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function remapSSH(ref,files,path,lon)
    sshAvi=nan*lon.avi;
    Tf=disp_progress('init','files');
    FF=thread_distro(12,numel(files.in));
    spmd
        ff=FF(labindex,1):FF(labindex,2) ;
        for f=ff
            Tf=disp_progress('c',Tf,numel(ff));
            try
                overFiles(files.in(f).name,path,ref,sshAvi);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function overFiles(file,path,ref,sshPOP2AVI)
    outfile=[path.out file(1:end-3) '.mat']; 
    
    ssh=nc_varget([path.in file],'SSH');
    Tp=disp_progress('init','averaging');
    for ii=1:numel(ref)
        Tp=disp_progress('c',Tp,numel(ref),10);
        sshPOP2AVI(ii)=nanmean(ssh(ref{ii}));
    end
    save(outfile,'sshPOP2AVI')
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lat,lon]=inits(ncfile)
    %%
    lat.avi=nc_varget(ncfile.avi,'lat');
    lon.avi=nc_varget(ncfile.avi,'lon');
    [lon.avi,lat.avi]=meshgrid(lon.avi,lat.avi);
    %%
    lat.pop=nc_varget(ncfile.pop,'U_LAT_2D');
    lon.pop=nc_varget(ncfile.pop,'U_LON_2D');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ref]=makeXref(idx)
    if ~exist('ref.mat','file')
        ref=makeXrefCalc(idx);
    else
        load ref;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ref=makeXrefCalc(idx)
    ref=repmat({[]},Ya*Xa,1); % crossref cell arr
    [idxS,idxO]=sort(idx); %
    UN=unique(idxS); %
    [~,bin]=histc(idxS,UN); %
    fbin=find(diff(bin)); % find bin bndries
    lims =[ [1; fbin+1] [fbin; numel(bin)] ]; % from:till lims (with respect to sorted idx vec from pop)
    Tp=disp_progress('init','gffcztc');
    for ii=1:numel(UN)
        Tp=disp_progress('c',Tp,numel(UN),100);
        popx=idxO(lims(ii,1):lims(ii,2));
        avix=UN(ii);
        ref{avix}=popx;
    end
    save avi2pop ref
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [idx]=CrossRef2Closest(lon,lat)
    azi=deg2rad(lon.pop(:));
    elev=deg2rad(lat.pop(:));
    [x,y,z] = sph2cart(azi,elev,1);
    qazi= deg2rad(lon.avi(:));
    qelev= deg2rad(lat.avi(:));
    [qx,qy,qz] = sph2cart(qazi,qelev,1);
    inxyz=[x,y,z];
    outxyz=[qx,qy,qz];
    JJ=thread_distro(12,numel(azi));
    if exist('idx.mat','file')
        load idx
    else
        spmd
            labindex
            myII=JJ(labindex,1):JJ(labindex,2);
            idx = dsearchn(outxyz,inxyz(myII,:));
            idx = gcat(idx,1,1);
        end
        idx=idx{1};
        save idx idx
    end
end
