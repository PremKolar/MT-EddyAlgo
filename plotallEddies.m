a = load('../dataaviI/EDDIES/EDDIE_19940105_-80-+80_000-360.mat');
b = load('../dataaviI/CUTS/CUT_19940105_-80-+80_000-360.mat');
%%
LA = b.fields.lat;
LO = wrapTo180(b.fields.lon);
%%
close all
hold on
% axesm('aitoff','MapLatLimit',[-80 80],'MapLonlimit',[-180 180],'frame','off','grid', 'off')
axesm('aitoff','frame','on','grid', 'off')
load coast
plotm(lat,long)
%  axesm('stereo','origin',[-90 0 0],'MapLatLimit',[-90 -10],'frame','on','flinewidth',1,'grid', 'on')

for cc = 1:1:numel(a.AntiCycs)
    fprintf('%d promil\n',round(1000*cc/numel(a.AntiCycs)))
    coo = a.AntiCycs(cc).coor.exact;   
    la = interp2(LA,coo.x,coo.y);
    lo = interp2(LO,coo.x,coo.y);
    flag = diff(la)>10;
    la([false flag']) = nan;
    lo([false flag']) = nan;
    plotm(la,lo,'r')
end
for cc = 1:1:numel(a.Cycs)
    fprintf('%d promil\n',round(1000*cc/numel(a.Cycs)))
    coo = a.Cycs(cc).coor.exact;   
    la = interp2(LA,coo.x,coo.y);
    lo = interp2(LO,coo.x,coo.y);
    flag = diff(la)>10;
    la([false flag']) = nan;
    lo([false flag']) = nan;
    plotm(interp2(LA,coo.x,coo.y),interp2(LO,coo.x,coo.y),'g')
end
axis tight
