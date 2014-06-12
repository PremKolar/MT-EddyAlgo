function U=input_vars
	%% threads / debug
	U.threads.num=12;
	U.debugmode=0;
	%% time
%     U.time.from.str='19920114';  % min aviso
%     U.time.till.str='20130807'; % max aviso
% 	U.time.from.str='19940105';  % min pop
% 	U.time.till.str='20061231'; % max pop   
    U.time.from.str='19940105';  % min pop
	U.time.till.str='19960105'; % max pop   
    U.time.delta_t=7; % [days]!
	%% dirs
	U.path.OutDirBaseName='avitestTiny';
% 	U.path.TempSalt.name='../TempSalt/';
	U.path.TempSalt.name='/home/niko/ROMnew/TempSalt/';
	U.path.raw.name='/data/icdc/ocean/aviso_ssh/DATA/weekly/msla/';
	%% output MAP STUFF
	U.map.out.X=20*1+1;
	U.map.out.Y=30*1+1;
	U.map.out.west=170;
	U.map.out.east=180;
	U.map.out.south=-60;
	U.map.out.north=-40;
    %% input MAP STUFF
	U.map.in.west=U.map.out.west;
	U.map.in.east=U.map.out.east;
	U.map.in.south=U.map.out.south;
	U.map.in.north=U.map.out.north;
	U.map.in.time.delta_t = 1; % [days]
	U.map.in.ssh_unitFactor = 100; % eg 100 if SSH data in cm, 1/10 if in deka m etc..
    %% input patterns
    U.map.in.fname='SsaltoDuacs__merged_msla__AVISO__ref__0.333deg__yyyymmdd.nc';
	U.map.in.keys.lat='lat';
	U.map.in.keys.lon='lon';
	U.map.in.keys.ssh='msla'; 
    U.TS.keys.lat='U_LAT_2D';
    U.TS.keys.lon='U_LON_2D';  
    U.TS.keys.salt='SALT';
    U.TS.keys.temp='TEMP';  
    U.TS.keys.depth='depth_t';  
    %% thresholds
	U.contour.step=0.01; % [SI]
	U.thresh.ssh_filter_size=1;
	U.thresh.radius=0; % [SI]
	U.thresh.amp=0.01; % [SI]
	U.thresh.shape.iq=0.1; % isoperimetric quotient
	U.thresh.shape.chelt=0.3; % (diameter of circle with equal area)/(maximum distance between nodes) (if ~switch.IQ)
	U.thresh.corners=4; % min number of data points for the perimeter of an eddy
	U.thresh.dist=1*24*60^2; % max distance travelled per day
	U.thresh.life=3; % min num of living days for saving
	U.thresh.ampArea=[.25 2.5]; % allowable factor between old and new time step for amplitude and area (1/4 and 5/1 ??? chelton)
	%% switches
	U.switchs.RossbyStuff=true; % TODO
	U.switchs.IQ=false;
	U.switchs.chelt=false;
	U.switchs.distlimit=0;
	U.switchs.AmpAreaCheck=0;
	U.switchs.netUstuff=1;
	%% parameters
	U.parameters.rossbySpeedFactor=1.75; % only relevant if cheltons method is used. eddy translation speed assumed factor*rossbyWavePhaseSpeed for tracking projections
	U.parameters.meanU=100; % depth from which to take mean U
	U.parameters.minProjecDist=150e3; % minimum linear_eccentricity*2 of ellipse (see chelton 2011)
	U.parameters.trackingRef='CenterOfVolume'; % choices: 'centroid', 'CenterOfVolume', 'Peak'
	%% technical params
	U.RossbyStuff.splits =4; % number of chunks for brunt väis calculations
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
		'amp.to_mean';
		};
	
	%% fields 4 colorcoded track plots
	U.FieldKeys.trackPlots= { ...
		'isoper';
		'radius.mean';
		'radius.meridional';
		'radius.zonal';
		'age';
		'peak.amp.to_contour';
		'peak.amp.to_mean';
		'peak.amp.to_ellipse';
		};
	%% TODO
	U.FieldKeys.senses= { ...
		'AntiCycs';
		'Cycs';
		};
end
