function U=input_vars
  %% threads / debug
    U.threads.num=12;
    U.debugmode=1;
    %% time
	U.time.from.str='19091231';
	U.time.till.str='19110627';
	U.time.delta_t=3; % [days]!
	U.path.OutDirBaseName='MAD';
    U.path.TempSalt.name='../TempSalt/';
    U.path.raw.name='../';
   %% output MAP STUFF   
	U.map.out.X=46*1+1;
	U.map.out.Y=37*1+1;		
	U.map.out.west=0;
	U.map.out.east=46;
	U.map.out.south=0;
	U.map.out.north=37;
	%% input MAP STUFF   
    U.map.in.west=U.map.out.west;
    U.map.in.east=U.map.out.east;
    U.map.in.south=U.map.out.south;
    U.map.in.north=U.map.out.north;
    U.map.in.time.delta_t = 1; % [days]
    U.map.in.SSH_unitFactor = 100; % eg 100 if SSH data in cm, 1/10 if in deka m etc..
    U.map.in.pattern.fname='psvar.cdf';
    U.map.in.pattern.lat='U_LAT_2D';
    U.map.in.pattern.lon='U_LON_2D';
    U.map.in.pattern.ssh='SSH';
    %% thresholds
    U.contour.step=0.01; % [SI]
    U.thresh.ssh_filter_size=1;
    U.thresh.radius=0; % [SI]
    U.thresh.amp=0.01; % [SI]
    U.thresh.shape.iq=0.3; % isoperimetric quotient
    U.thresh.shape.chelt=0.3; % (diameter of circle with equal area)/(maximum distance between nodes) (if ~switch.IQ)
    U.thresh.corners=6; % min number of data points for the perimeter of an eddy
    U.thresh.dist=5*24*60^2; % max distance travelled per day
    U.thresh.life=3; % min num of living days for saving
    U.thresh.ampArea=[.25 2.5]; % allowable factor between old and new time step for amplitude and area (1/4 and 5/1 ??? chelton)
    %% switches
    U.switchs.RossbyStuff=false; % TODO
    U.switchs.IQ=true;
    U.switchs.chelt=false;
    U.switchs.distlimit=false;
    U.switchs.AmpAreaCheck=false;
    %% parameters
    U.parameters.rossbySpeedFactor=1.75; % only relevant if cheltons method is used. eddy translation speed assumed factor*rossbyWavePhaseSpeed for tracking projections
    U.parameters.meanU=100; % depth from which to take mean U
    U.parameters.minProjecDist=150e3; % minimum linear_eccentricity*2 of ellipse (see chelton 2011)
    U.parameters.trackingRef='CenterOfVolume'; % choices: 'centroid', 'CenterOfVolume', 'Peak'
    %% technical params
    U.RossbyStuff.splits =36; % number of chunks for brunt väis calculations
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
        'peak.amp.to_mean';
        %         'peak.amp.to_mean.of_contour';
        'peak.amp.to_ellipse';
        };
    %% TODO
    U.FieldKeys.senses= { ...
        'AntiCycs';
        'Cycs';
        };
end
