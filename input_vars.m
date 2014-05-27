function U=input_vars
	%% threads
	U.threads.num=2;
	%% time
	U.time.from.str='19091231';
	U.time.till.str='19160627';
	U.time.delta_t=3; % [days]!
	%% dirs
 	U.path.root='../MAD/';
 	U.path.plots='../plotsMAD/';
	 U.path.TempSalt.name='/media/ROM/TempSalt/';
 	U.path.raw.name=U.path.root;
	%% dims for map plots
	U.dim.X=46*1+1;
	U.dim.Y=37*1+1;		
	U.dim.west=0;
	U.dim.east=46;
	U.dim.south=0;
	U.dim.north=37;
	U.dim.NumOfDecimals=1;
	%% thresholds
	U.contour.step=0.01; % [SI]
	U.thresh.radius=1e4; % [SI]
	U.thresh.amp=0.01; % [SI]
	U.thresh.shape.iq=0.3; % isoperimetric quotient
	U.thresh.shape.chelt=0.5; % (diameter of circle with equal area)/(maximum distance between nodes) (if ~switch.IQ) 
	U.thresh.corners=6; % min number of data points for the perimeter of an eddy
	U.thresh.dist=.3*24*60^2; % max distance travelled per day
	U.thresh.life=9; % min num of living days for saving
	U.thresh.ampArea=[.25 2.5]; % allowable factor between old and new time step for amplitude and area (1/4 and 5/1 ??? chelton)
	
	
	
	
	%% best ignore all from here
	%###########################################################################
	%###########################################################################	
	%% thresholds
	U.thresh.ssh_filter_size=1;
	U.thresh.ampArea=[.25 2.5]; % allowable factor between old and new time step for amplitude and area (1/4 and 5/1 ??? chelton)	
	%% switches
	U.switchs.IQ=true;	
	U.switchs.RossbyStuff=false; 
	U.switchs.chelt=false;		
	 U.switchs.distlimit=false;
	 U.switchs.AmpAreaCheck=false;
	 %% parameters
	 U.parameters.rossbySpeedFactor=1.75; % only relevant if cheltons method is used. eddy translation speed assumed factor*rossbyWavePhaseSpeed for tracking projections
    U.parameters.depthRossby=100; % depth from which to take rossby phase speed and radius
	 U.parameters.minProjecDist=150e3; % minimum  linear_eccentricity*2 of ellipse (see chelton 2011)
	 U.parameters.trackingRef='CenterOfVolume'; % choices: 'centroid', 'CenterOfVolume', 'Peak'
    %% technical params
    U.RossbyStuff.splits = 12; % number of chunks for brunt väis calculations
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
        %'radius.volume';
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
