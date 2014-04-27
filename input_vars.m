function U=input_vars
	%% threads
	U.threads.num=12;
	%% time
 	U.time.from.str='19940425';
 	U.time.till.str='19980131';
 	U.time.delta_t=1; % [days]!
	%% dirs
	U.path.root='../datatemp/';
	U.path.TempSalt.name='TempSalt/';
	U.path.raw.name='/scratch/uni/ifmto/u241194/DAILY/EULERIAN/SSH/';
 	U.pattern.in='CUT_yyyymmdd_SSSSsNNNNnWWWWwEEEEe.mat'; %shouldnt be here
	%% thresholds
	U.contour.step=0.01; % [SI]
	U.thresh.ssh_filter_size=1;
	U.thresh.radius=1e4; % [SI]
	U.thresh.amp=0.01; % [SI]
	U.thresh.shape.iq=0.35; % isoperimetric quotient
	U.thresh.shape.chelt=0.5; % (diameter of circle with equal area)/(maximum distance between nodes) (if ~switch.IQ) 
	U.thresh.corners=8; % min number of data points for the perimeter of an eddy
	U.thresh.dist=.3*24*60^2; % max distance travelled per day
	U.thresh.life=10; % min num of living days for saving
	%% dims for map plots
	U.dim.X=90*2+1;
	U.dim.Y=40*2+1;
	U.dim.west=-90;
	U.dim.east=0;
	U.dim.south=20;
	U.dim.north=60;
	U.dim.NumOfDecimals=1;
	%% switches
	U.switchs.RossbyStuff=false;
	U.switchs.IQ=true;	
	
	%% technical params
	U.RossbyStuff.splits = 10; % number of chunks for brunt väis calculations
	%% fields that must end with .mean and .std - for output plot maps
	U.FieldKeys.MeanStdFields= { ...
		'age';
		'sense';
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
