function U=input_vars
	%% threads
	U.threads.num=12;
	%% time
 	U.time.from.str='19941030';
 	U.time.till.str='19951030';
 	U.time.delta_t=1; % [days]!
	%% dirs
	U.path.TempSalt.name='../TempSalt/';
 	U.path.raw.name='/scratch/uni/ifmto/u241194/DAILY/EULERIAN/SSH/';
  	U.path.root='../dataChelt2/';
    U.path.plots='../plotsChelt/';
 	%% thresholds
	U.contour.step=0.01; % [SI]
	U.thresh.ssh_filter_size=1;
	U.thresh.radius=0; % [SI]
	U.thresh.amp=0.01; % [SI]
	U.thresh.shape.iq=0.3; % isoperimetric quotient
	U.thresh.shape.chelt=0.2; % (diameter of circle with equal area)/(maximum distance between nodes) (if ~switch.IQ) 
	U.thresh.corners=6; % min number of data points for the perimeter of an eddy
	U.thresh.dist=.8*24*60^2; % max distance travelled per day
	U.thresh.life=10; % min num of living days for saving
	%% dims for map plots
	U.dim.X=20*1+1;
 	U.dim.Y=20*2+1;
   	U.dim.west=-80;
 	U.dim.east=-60;
 	U.dim.south=20;
 	U.dim.north=40;
	U.dim.NumOfDecimals=1;
	%% switches
	U.switchs.RossbyStuff=false;  % TODO
	U.switchs.IQ=false;	
	U.switchs.chelt=true;	
	
	%% technical params
	U.RossbyStuff.splits = 10; % number of chunks for brunt väis calculations
	%% fields that must end with .mean and .std - for output plot maps
	U.FieldKeys.MeanStdFields= { ...
		'age';
		'dist.traj.fromBirth';
		'dist.traj.tillDeath';
		'dist.zonal.fromBirth';
		'dist.zonal.tillDeath';
		'dist.merid.fromBirth';
		'dist.merid.tillDeath';
		'radius.mean';
		'radius.zonal';
		'radius.meridional';
		'vel.traj';
		'vel.zonal';
		'vel.merid';
		'amp.to_contour';
		'amp.to_ellipse';
		 'amp.to_mean.of_contour';
		};
	
	%% fields 4 colorcoded track plots
	U.FieldKeys.trackPlots= { ...
		'isoper';
		'radius.mean';
		'radius.meridional';
		'radius.zonal';
		%'radius.volume';
		'age';
		'peak.amp.to_contour';
      'peak.amp.to_mean.of_contour';
		'peak.amp.to_ellipse';
		};	
end
