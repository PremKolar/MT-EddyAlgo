function sub09_trackstuff
	load S09main II DD T
	%%
	try
		TR=getTR(DD);
	catch me
		disp(me.message)
		sub09_trackinit;
		TR=getTR(DD) ;
	end
	%%
	senses=DD.FieldKeys.senses;
	catsen= @(f) [TR.(senses{1}).(f); TR.(senses{2}).(f) ];
	S.t2l=@(t) round(linspace(t(1),t(2),t(3)));
	%%
	rad=round(catsen('rad')/1000);
	vel=catsen('vel')*100;
	age=catsen('age');
	lat=catsen('lat');
	lon=catsen('lon'); %#ok<NASGU>
	%%
	S.rightyscalenum=5;
	age(end+1:end+S.rightyscalenum)=max(age)-0;
	lat(end+1:end+S.rightyscalenum)=S.t2l([min(lat) max(lat) S.rightyscalenum]);
	rad(end+1:end+S.rightyscalenum)=S.t2l([min(rad) max(rad) S.rightyscalenum]);
	vel(end+1:end+S.rightyscalenum)=10;
	%%
	[~,sml2lrg] = sort(rad)  ;
	S.age=age(fliplr(sml2lrg));
	S.lat=lat(fliplr(sml2lrg));
	S.rad=rad(fliplr(sml2lrg));
	S.vel=vel(fliplr(sml2lrg));
	%%
	zerage = S.age<=0  ;
	velHigh= S.vel>20 | S.vel <-30;
	radnill = isnan(S.rad) | S.rad==0;
	killTag= zerage | velHigh | radnill ;
	S.age(killTag)=[];
	S.lat(killTag)=[];
	S.rad(killTag)=[];
	S.vel(killTag)=[];
	%%
	scattStuff(S,T,DD,II)
	%%
	velZonmeans(S,DD,II,T)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function velZonmeans(S,DD,II,T) %#ok<INUSD>
	close all
	LA     = round(S.lat);
	LAuniq = unique(LA)';
	vvM=nan(size(LAuniq));
	vvS=nan(size(LAuniq));
	for cc=1:(numel(LAuniq))
		vvM(cc)=mean(S.vel(LA==LAuniq(cc)));
		vvS(cc)=std(S.vel(LA==LAuniq(cc)));
	end
	vvM(abs(LAuniq)<5)=nan;
	vvS(abs(LAuniq)<5)=nan;
	%%
	h.own=ownPlot(DD,II,LAuniq,vvM,vvS);
	savefig(DD.path.plots,100,800,800,['S-velZonmean'],'dpdf');
	%%
	chelt = imread('/scratch/uni/ifmto/u300065/presMT/FIGS/chelt11Ucomp.jpg');
	chelt= chelt(135:3595,415:3790,:);
	h.ch=chOverLay(S,DD,chelt,LAuniq,vvM);
	savefig(DD.path.plots,100,800,800,['S-velZonmean4chelt11comp'],'dpdf');
	
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h=ownPlot(DD,II,LAuniq,vvM,vvS)
	%%
	clf
	lw=2;
	pp(1)=plot(II.la(:,1),-II.maps.zonMean.Rossby.small.phaseSpeed*100); 	hold on
	pp(2)=plot(LAuniq,-vvM,'r');
	pp(4)=plot(LAuniq,-vvM+vvS,'y');
	pp(5)=plot(LAuniq,-vvM-vvS,'y');
	pp(3)=plot([DD.map.out.south DD.map.out.north], [0 0],'b--');
	axis([-70 70 -5 20])
	set(pp(1:3),'linewidth',lw)
	leg=legend('Rossby-wave phase-speed',2,'all eddies',2,'std');
	set( get(leg,'children'),'linewidth',lw)
	ylabel('[cm/s]')
	xlabel('[latitude]')
	title(['westward propagation [cm/s]'])
	set(get(gcf,'children'),'linewidth',lw)
	grid minor
	h=gcf;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h=chOverLay(S,DD,chelt,LAuniq,vvM)
	[Y,X,Z]=size(chelt);
	ch=reshape(flipud(reshape(chelt,Y,[])),[Y,X,Z]);
	ch=permute(reshape(fliplr(reshape(permute(ch,[3,1,2]),1,[])),[Z,Y,X]),[2,3,1]);
	ch=permute(reshape(flipud(reshape(permute(ch,[2,1,3]),[X,Y*Z])),[X,Y,Z]),[2,1,3]);
	%%
	vvm=vvM+15;
	lau=LAuniq;
	kill=isnan(lau) | isnan(vvM) | abs(lau)<10 | abs(lau)>50;
	lau(kill)=nan;
	vvm(kill)=nan;
	
	%%
	clf
	imagesc(linspace(-50,50,X),linspace(-5,20,Y),ch)
	hold on
	mid=floor(numel(lau)/2);
	x.a=lau(1:mid);
	x.b=lau(mid+2:end);
	y.a=vvm(1:mid);
	y.b=vvm(mid+2:end);
	plot(x.a,spline(x.a,y.a,x.a),'g-','linewidth',1);
	plot(x.b,spline(x.b,y.b,x.b),'g-','linewidth',1);
	plot(lau,vvm,'g.','markersize',8);
	set(gca, 'yticklabel', flipud(get(gca,'yticklabel')));
	axis tight
	grid minor
	ylabel('[cm/s]')
	xlabel('[latitude]')
	title(['westward propagation [cm/s]'])
	text(-45,-2,['R/L: ' num2str(DD.thresh.maxRadiusOverRossbyL)])
	text(-45,-1,['iq: ' num2str(DD.thresh.shape.iq)])
	text(-45,0,['life: ' num2str(DD.thresh.life)])
	text(-45,1,['id: ' sprintf('%02d%%',-100+100*DD.thresh.IdentityCheck)])
	text(-45,2,['vrtcs: ' sprintf('%2d',DD.thresh.corners.min)])
	text(-45,3,['days: ' sprintf('%5d',DD.time.span)])
	text(-45,4,['dt: ' sprintf('%2d',DD.time.delta_t)])
	text(-45,5,['meanAge: ' sprintf('%5d',round(mean(S.age)))])
	text(-45,6,['data: ' sprintf('%8d',numel(S.age))])
	h=gcf;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scattStuff(S,T,DD,II)
	age=S.age;
	lat=S.lat;
	vel=S.vel;
	rad=round(S.rad);
	inc=1;
	%%
	oie=@(inc,x) x(1:inc:end);
	incscatter=@(inc,a,b,c,d) scatter(oie(inc,a),oie(inc,b),oie(inc,c),oie(inc,d));
	incscatter(inc,age,lat,rad,vel);
	grid on
	axis tight
	set(gca,'XAxisLocation','bottom')
	cb=colorbar;
	cb1 = findobj(gcf,'Type','axes','Tag','Colorbar');
	cbIm = findobj(cb1,'Type','image');
	alpha(cbIm,0.5)
	set(cb,'location','north','xtick',(S.t2l(T.vel)),'xlim',T.vel([1 2]))
	doublemap([T.vel(1) 0 T.vel(2)],II.aut,II.win,[.9 1 .9],20)
	h1=gca;
	h1pos = get(h1,'Position'); % store position of first axes
	h2 = axes('Position',h1pos,...
		'XAxisLocation','top',...
		'YAxisLocation','right',...
		'Color','none');
	set(h2, ...
		'ytick',linspace(0,1,S.rightyscalenum),...
		'xtick',[],...
		'yticklabel',(S.t2l([min(rad) max(rad) S.rightyscalenum])))
	set(h1, ...
		'ytick',[-50 -30 -10 0 10 30 50],...
		'xtick',S.t2l(T.age))
	ylabel(h2,'radius [km]')
	ylabel(h1,'lat  [{\circ}]')
	xlabel(h1,'age [d]')
	xlabel(h2,'zon. vel.  [cm/s] - eddies beyond scale dismissed!')
	set(get(gcf,'children'),'clipping','off')
	%%
	savefig(DD.path.plots,T.rez,T.width,T.height,['sct-ageLatRadU'],'dpdf');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TR=getTR(DD)
	xlt=@(sen,f) extractfield(load(['TR-' sen '-' f '.mat']),'tmp');
	F={'rad','age','lat','lon'};
	g=@(c) cat(1,c{:});
	for ss=1:2
		for fi=1:numel(F);f=F{fi};
			sen=DD.FieldKeys.senses{ss};
			TR.(sen).(f)=((xlt(sen,f))');
		end
		f='vel';
		TR.(sen).(f)=g(g(xlt(sen,f)));
	end
end
