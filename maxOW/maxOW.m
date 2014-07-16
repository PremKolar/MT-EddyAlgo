%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 15-Jul-2014 13:52:44
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maxOW
	%% set up
	[DD]=maxOWsetUp;
	%% spmd
	maxOWmain(DD);
	%% make netcdf
	maxOWwriteNCfile(DD);
    %% 
    OWprocess(DD)
	%% update DD
	save_info(DD);
end

function OWprocess(DD)   
  load([DD.path.Rossby.name 'filenames.mat'])
  load(files.OkuboWeiss,'out')
  OW=out;
  OW(abs(OW)>1)=nan;
  [T,Z,Y,X]=size(OW);
  for t=1:T
     [owm,zi]=nanmin(squeeze(OW(t,:,:,:) ),[], 1);
     owm=squeeze(owm);
     zi=squeeze(zi);
     figure
    shearDom=owm>=0;
     owm(shearDom)=nan;
     zi(shearDom)=nan;
    ppc(log10(-owm))
    caxis([-7 -3])
    figure
    ppc(zi)
  end
  
end