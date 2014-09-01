%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 15-Jul-2014 13:52:44
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function compPlots
    dirs={'iq2';'iq6';'iq6dl';'iq6dlMr';'iq6dlMrIc';'iq8';'iq4';'ch'};    
    D=inIt(dirs);
   
    
%     L=numel(dirs)
%     spmd(L)
%         ii=labindex
%           longesttracks(D.basedir,D.out(ii),D.thresh.ampArea,ii);             
%     end
    
     for ii=1:numel(D.out)          
        longesttracks(D.basedir,D.out(ii),D.thresh.ampArea,ii);   
     end
   printoutsLT(D);
   
     
     
    
    %%
     LTscatter('longesttracks.mat')
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LTscatter(matname) 
    load(matname)
    
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
    D.basedir=['/scratch/uni/ifmto/u300065/FINAL/aorStuff/'];
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
