%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 15-Jul-2014 13:52:44
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compPlots
    dirs={'iq2';'iq6';'iq6dl';'iq6dlMr';'iq6dlMrIc';'iq8';'iq4';'ch'};
    dirs={'iq2'};
    D=inIt(dirs);
    
    
    %     L=numel(dirs)
    %     spmd(L)
    %         ii=labindex
    %           longesttracks(D.basedir,D.out(ii),D.thresh.ampArea,ii);
    %     end
    
    for ii=1:numel(D.out)
        %          longesttracks(D.basedir,D.out(ii),D.thresh.ampArea,ii);
        LTscatter(D.basedir,D.out(ii),D.thresh.ampArea,ii)
    end
    %    printoutsLT(D);
    
    
    
    
    %%
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LTscatter(basedir,dout,threshampArea,ii)
    
    
    
    all=collectScatD(basedir,dout);
    
   
   
    scatter(cat(2,all.aol),cat(2,all.lat))
    scatter(cat(2,all.iq),cat(2,all.lat))
    scatter(cat(2,all.iq),cat(2,all.lat),10,cat(2,all.aol))
   colorbar
   caxis([0 6])
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function all=collectScatD(basedir,dout)
    fs=dir([basedir 'data' dout.name '/TRACKS/*.mat']);
   
    ff=1
    while ff<numel(fs)
        try
            file=[basedir 'data' dout.name '/TRACKS/' fs(ff).name];
            MF=matfile(file);
            trck=MF.trck;
            all(ff).lat=extractdeepfield(trck,'geo.lat');
            all(ff).aol=extractdeepfield(trck,'area.RadiusOverRossbyL');
            all(ff).iq=extractdeepfield(trck,'isoper');
            all(ff).peak2contour=extractdeepfield(trck,'peak.to_contour');
            all(ff).peak2mean=extractdeepfield(trck,'peak.to_mean');
            all(ff).peak2ellip=extractdeepfield(trck,'peak.to_ellipse');
        catch
            continue
        end
        ff=ff+1
    end
   
    
    
      ff=1
    while ff<numel(fs)
        try
            file=[basedir 'data' dout.name '/TRACKS/' fs(ff).name];
            MF=matfile(file);
            trck=MF.trck;
            all(ff).peak2contour=extractdeepfield(trck,'peak.to_contour');
            all(ff).peak2mean=extractdeepfield(trck,'peak.to_mean');
            all(ff).peak2ellip=extractdeepfield(trck,'peak.to_ellipse');
        catch
            continue
        end
        ff=ff+1
    end
    
    
    
    
    
end
function scatplots(LTdata)
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D=inIt(dirs)
    addpath(genpath('./'))
    flsh= @(x) deal(x{:});
    dbstop if error
    dbstop if warning
    D=INPUT;
    D.threads.num=init_threads(12);
    D.here=pwd;
    D.basedir=['/scratch/uni/ifmto/u300065/FINAL/smAor/'];
    D.out(numel(dirs))=struct;
    [D.out(:).name]=flsh(dirs);
    [D.out(:).LTfile]=  flsh(cellfun(@(c) [D.basedir 'LT_' c '.mat'],dirs,'uniformoutput',false));
    %     [D.out(:).file]=  flsh(cellfun(@(c) [D.basedir 'tracks_' c '.mat'],dirs,'uniformoutput',false));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printoutsLT(D)
    system(sprintf('pdfjam -o tmp.pdf crpd*pdf'))
    outtit=[cat(2,D.out(:).name),'.pdf'];
    system(['pdfcrop  --margins "1 3 1 1" tmp.pdf ' outtit])
    delete('crpd*pdf')
    delete('tmp.pdf')
end
