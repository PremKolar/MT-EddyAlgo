%%
close all
clear all
cmc=jet(200)
cma=cmc(1:100,:)
cmb=cmc(101:end,:)

D=textscan(pwd,'%s','delimiter','/')
d=D{1}{end}

dd=getfield(dir(['../data' d '/EDDIES/*19941120*' ]),'name')
cc=getfield(dir(['../dataC/CUTS/*19941120*' ]),'name')
load(['../data' d '/EDDIES/' dd ])
load(['../dataC/CUTS/' cc ])
whos
%%
DD=load(['../data' d '/DD.mat' ])
%%
csteps=nanmin(grids.ssh(:)):.01:nanmax(grids.ssh(:))
contour(grids.ssh,csteps,'linewidth',1.5)
%%
save(['tmp-' d])
%%
axis on
grid on
set(gca,'xticklabel','','yticklabel','')
% colormap(hsv)
cl=[-.5 .7]
caxis(cl) 
doublemap([cl(1) 0 cl(2)],cma,flipud(cmb),[.1 .6 .1],20)
%%
hold on
for ee=1:numel(anticyclones)
   x=anticyclones(ee).coordinates.exact.x;
   y=anticyclones(ee).coordinates.exact.y;
   plot(x,y,'--','color','black','linewidth',2)
end
%%
hold on
for ee=1:numel(cyclones)
   x=cyclones(ee).coordinates.exact.x;
   y=cyclones(ee).coordinates.exact.y;
   plot(x,y,'red','linewidth',2)
end

%%

axis([120 320 325 475])

tit=[sprintf('IQ: %2.1f', DD.thresh.shape.iq) ' - ' ...
     sprintf('ID: %2.1f', DD.thresh.IdentityCheck) ' - ' ...
     sprintf('r/Lr: %2d', DD.thresh.maxRadiusOverRossbyL) ]
title([tit])
colorbar('location', 'southOutside')
colorbar off
%%
savefig('./',72,800,600, d)




