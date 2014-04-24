%close all
clear all
addpath(genpath('/home/zmaw/u300065/MASTER_THESIS/CODE/SUBS'))
addpath(genpath('/home/zmaw/u300065/MASTER_THESIS/mexcdf'))


%% 
south=50;
north=58;
west=1;
east=12;

from='19940430';
till='19950102';

path_in=('/scratch/uni/ifmto/u241194/DAILY/EULERIAN/SSH/');




%% get geo stuff
radius_earth=6371000;
timestep=0;
EddiesAC=[];

fs=dir([path_in, '*.nc']);
file=[path_in, fs(round(length(fs)/2)).name];
[SSH_c,LON_c,LAT_c,DY_c,DX_c,YDim,XDim,ystart,yend,xstart,xend,F]=get_geocoor(file, south, north, west, east, radius_earth);

LON = ncread(file,'U_LON_2D')';
LAT = ncread(file,'U_LAT_2D')';
SSH  = ncread(file,'SSH')';
 
    SSH_min=min(SSH_c(:));
    SSH_max=max(SSH_c(:));


% Set up axes
axesm ('globe','Grid', 'on');
view(60,60)
axis off



% Display coastline vectors
load coast
plotm(lat, long,'r','linewidth',3)

scalefac=20;

hold on
plotm(LAT_c(1:scalefac:end,1:scalefac:end),LON_c(1:scalefac:end,1:scalefac:end),'-','markersize',0.0002,'color','blue')
plotm(LAT_c(1:scalefac:end,1:scalefac:end)',LON_c(1:scalefac:end,1:scalefac:end)','-','markersize',0.0002,'color','blue')

% surfm(LAT(1:scalefac:end,1:scalefac:end),LON(1:scalefac:end,1:scalefac:end),SSH(1:scalefac:end,1:scalefac:end))
% 
%  LAT_w=LAT(:,end-10:end) ;
%  LAT_e=LAT(:,1:10) ;
%  LON_e=LON(:,1:10) ;
%  LON_w=LON(:,end-10:end) ; 
%  plotm(LAT_w,LON_w,'color','blue')
%  plotm(LAT_e,LON_e,'color','red')
%   plotm(LAT_w',LON_w','color','blue')
%  plotm(LAT_e',LON_e','color','red')
 
 scalefac=10;
 SSH(SSH==0)=nan;
%pcolor(SSH(1:scalefac:end,1:scalefac:end));
% pcolorm(LAT(1:scalefac:end,1:scalefac:end),LON(1:scalefac:end,1:scalefac:end),SSH(1:scalefac:end,1:scalefac:end));
%  
% pcolor(LON(1:scalefac:end,1:scalefac:end)); shading interp; colorbar
% 
% 
% lonmin=min(LON(:))
% lonmax=max(LON(:))






figure
projec='globe';

axesm('MapProjection',projec,'MapLatLimit',[south north],'MapLonLimit',[west east])
col={'red','blue'}
for s=1:2
    plot_tracks(gcf,POS{s},col{s})
end





GRD=nan(180,360);
[A,O]=size(GRD);
EDDY_SPEEDS.ANTICYC.MEAN_SPEED=cell(size(GRD));
EDDY_SPEEDS.CYC.MEAN_SPEED=cell(size(GRD));
EDDY_SPEEDS.GRID.LAT=(-89:1:90);
EDDY_SPEEDS.GRID.LON=(-179:1:180);

[VELOCITIES_mean,VELOCITIES_median,vel_daily_total,LIFESPAN,LAT_mean,LON_mean]=get_eddy_statistics(POSITIONS,LAT,LON,radius_earth,TIMES);


for s=1:2
    LA=round(LAT_mean{s});
    LO=round(LON_mean{s});
    VE=VELOCITIES_mean{s};
    
    %
    %
    % for 1:length(VE)
    %
    % end
    
    
    
end









for la= EDDY_SPEEDS.GRID.LAT
    for lo= EDDY_SPEEDS.GRID.LON
        
        jj=find(LAT_mean{s}==la & LON_mean{s}==lo);
        
        if ~isempty(jj)
            sdfg
        end
        
        
    end
end













a=full(AREAS{1}(2:end,:));
a(a==0)=nan;
radius_thresh=sqrt(area_threshold/pi);

b=full(AREAS{2}(2:end,:));
b(b==0)=nan;


aa=[sqrt(a(:)/pi)]';

[hh,aa]=hist(a(:),((radius_thresh:5:100)*pi).^2)
[hb,bb]=hist(b(:),((radius_thresh:5:100)*pi).^2)
bar([hh;hb]')
L=get(gca,'xtick')
set(gca,'xtick',[]);
xticks= (linspace(1,sqrt(20),10)).^2;
xlabs=num2str([round(sqrt(aa(round(xticks))/pi))]');
set(gca,'xtick',xticks);
set(gca,'xtickLabel',xlabs);

bar([hh;hb])
colorbar
L=get(gca,'xtick')
set(gca,'xticklabel',['anti-cyc';'cyclonic']);
xticks= (linspace(1,sqrt(20),10)).^2;
xlabs=num2str([round(sqrt(aa(round(xticks))/pi))]');
colorbar('Ytick',xticks,'yticklabel',xlabs)
title('eddy scale, radius [m]')
g=(gcf)

savefig('eddy_scale_histo_indian.pdf',g,'pdf')







lvec=(round(south):1:round(north));

VEL_mean_binned=cell(1,2);
VEL_std_binned=cell(1,2);
for s=1:2
    L=round(LAT_mean{s});
    V=VELOCITIES_mean{s};
    V(V==inf)=nan;
    Vb=nan(size(lvec));
    Vs=nan(size(lvec));
    k=0;
    for lat=lvec(1:end)
        k=k+1;
        l=find(L==lat )
        Vb(k)=nanmean(V(l));
        Vs(k)=nanstd(V(l));
    end
    
    figure
    
    hold on
    plot(LAT_mean{s},VELOCITIES_mean{s},'.','color','black','markersize',1)
    VEL_mean_binned{s}=Vb;
    VEL_std_binned{s}=Vs;
    plot(lvec,Vb,lvec,Vb+Vs,'r',lvec,Vb-Vs,'r')
    switch s
        case 1
            title('mean anti-cyc. eddy propagation speed')
        case 2
            title('mean cyclonal eddy propagation speed')
    end
    xlabel('latitude')
    ylabel('[km/day]')
    axis([-70 20 0 30])
    % savefig(['eddy_speeds_indian_',num2str(s),'.pdf'],gcf,'pdf');
end






figure(20)
hist(LIFESPAN{1},(1:10:365))

hist(vel_daily_total,1000)
















% get last SSH

SSH  = ncread(RESULTS.fileinfo.file,'SSH')'; SSH = SSH(y_start:y_end,x_start:x_end);
figure(20)
% resize for plot


pixels=335*503
horvar_ratio=XDim/YDim
YDim_new=sqrt(pixels/horvar_ratio);
XDim_new=pixels/YDim_new;
Y_inc=round(YDim/YDim_new);
X_inc=round(XDim/XDim_new);
SSH_plot=SSH(1:Y_inc:end,1:X_inc:end);
LAT_plot=LAT(1:Y_inc:end,1:X_inc:end);
LON_plot=LON(1:Y_inc:end,1:X_inc:end);

DX_plot=DX(1:Y_inc:end,1:X_inc:end);
DY_plot=DY(1:Y_inc:end,1:X_inc:end);




south=min(LAT(:));
north=max(LAT(:));

west=min(LON(:));
east=max(LON(:));
ppp=figure(111);

projec='blub';

axesm('MapProjection',projec,'MapLatLimit',[south north],'MapLonLimit',[west east])
contourfm(LAT,LON,SSH,(min(SSH(:)):1:max(SSH(:))));
%contourfm_ssh(LAT_plot,LON_plot,SSH_plot,south,north,west,east,111,projec)


figure
load coast
hold on
axesm('MapProjection',projec,'MapLatLimit',[south north],'MapLonLimit',[west east])
plotm(lat,long)

hold on
for s=1:2
    %       axesm('MapProjection',projec,'MapLatLimit',...
    %          [south north],'MapLonLimit',[west east]...
    %           ,'Grid','on')
    %grid on
    sense=senses(s);
    %     %
    switch sense
        case 1
            co='red';
        case -1
            co='green';
    end
    %
    %
    POS=POSITIONS{s};
    
    for n=1:length(POS(1,:))
        t = find(POS(:,n),1,'last');
        plotm(LAT(POS(2:t,n)),LON(POS(2:t,n)),'linewidth',1.5,'color','black')
    end
    
    
    
    %     length_total=length(Rims_all{s});
    %
    %     k=1;
    %     k_total=k;
    %     while k<length_total
    %         length_eddy=Rims_all{s}(k);
    %         eddy_start=k+1;
    %         eddy_end  =k+length_eddy;
    %         hold on;
    %         plotm(LAT(Rims_all{s}(eddy_start:eddy_end)),LON(Rims_all{s}(eddy_start:eddy_end)),'-','linewidth',2,'color',co)
    %         k=k+length_eddy+1;
    %     end
    %
    %     P=POSITIONS{s};
    %     Pe=Peaks{s};
    %     for n=1:length(P(1,:))
    %         t = find(P(:,n),1,'last');
    %         if t==size(P,1)
    %             x1d=full(P(t,n));
    %             hold on
    %             plotm(LAT(x1d),LON(x1d),'b*','markersize',7);
    %         end
    %     end
    %
    
    
    %     figure(10)
    %     saveas(gcf, 'MAP_test.eps','eps');
    %
    %
    %       plot_tracks(111,POSITIONS{s})
    %         plot_eddy_perimeters(111,Rims_all{s},sense)
    %    plot_eddy_peaks(POSITIONS{s},LAT_plot,LON_plot,1)
    %      tag_eddy(POSITIONS{s},111,sense)
    %     showaxes
    %
end


axis equal

 
%%
path_out='~/MASTER_THESIS/RESULTS/SPEEDS/'
xdim=Xdim;
ydim=Ydim;
titles={'anti_cyc','cyc'};
for s=1:2
    V= VB{s};
    figure(s)
    figure('renderer','zbuffer')
    pcolor(GLO,GLA,V);shading flat;colorbar
    
    caxis([0 10])
    title(titles{s})
    xlabel('[m/s]');
    
    
    figpos=getpixelposition(gcf);
    resolution=get(0,'ScreenPixelsPerInch');
    rez=resolution;
    rez=300;
    
    
    set(gcf,'position',[0 (s-1)*ydim xdim ydim]);
    set(gcf,'paperunits','inch','papersize',[xdim ydim],'paperposition',[0 0 [xdim ydim]]);
    name=titles{s}
    
    
    savefig(gcf)
    %  print(gcf,fullfile(path_out,name),'-dpsc',['-r',num2str(rez)],'-zbuffer') %save file
    
end

%%

save(['/scratch/uni/ifmto/u300065/RESULTS/EDDY_SPEEDS/',dir_in]);


%%
return
matlabpool CLOSE
% Set up axes
axesm ('mercator','Grid', 'on');
%view(60,60)
%axis off

% Display coastline vectors
load coast
plotm(lat, long,'r','linewidth',1)
hold on
plot3m(LAT_stepped{1},LON_stepped{1},VEL_stepped{1},'.')

pcolorm(EDDY_SPEEDS.grid.lat,EDDY_SPEEDS.grid.lon,  EDDY_SPEEDS.anticyclonic.speeds_mean)



saveas(gcf,'eddy_speeds','png')
shading flat
savefig('speeds','eps','-r500')
xdim=1600*3
ydim=1200*3
f=gcf; %f is the handle of the figure you want to export
figpos=getpixelposition(f);
resolution=get(0,'ScreenPixelsPerInch');
rez=resolution;
rez=300
set(f,'paperunits','inch','papersize',[xdim ydim]/rez,'paperposition',[0 0 [xdim ydim]/rez])

path_out=['./']; %the folder where you want to put the file
if ~exist(path_out)
    mkdir(path_out)
end
name=['speeds_',dir_in]; %what you want the file to be called
print(f,fullfile(path_out,name),'-dpsc2',['-r',num2str(rez)],'-zbuffer') %save file
