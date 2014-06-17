function U=input_vars
	%% threads / debug
	U.threads.num=12;
	U.debugmode=false;
	%     U.debugmode=true;
	%% time
	U.time.delta_t=1; % [days]!
	U.time.from.str='19940102';
	% 	 U.time.from.str='19940425';
	%     U.time.till.str='19960730';
	U.time.till.str='20040102';
	%% dirs
	U.path.OutDirBaseName='mad';
	U.path.TempSalt.name='../TempSalt/';
	%     U.path.TempSalt.name='/media/ROM/TempSalt/';
	U.path.raw.name='../';
	%     U.path.raw.name='/media/ROM/SSH_POP/';
	%% output MAP STUFF
	U.map.out.X=46*1+1;
	U.map.out.Y=37*1+1;
	U.map.out.west=0;
	U.map.out.east=46;
	U.map.out.south=10;
	U.map.out.north=47;
	%% input MAP STUFF
	U.map.in.west=U.map.out.west;
	U.map.in.east=U.map.out.east;
	U.map.in.south=U.map.out.south;
	U.map.in.north=U.map.out.north;
	U.map.in.time.delta_t = 1; % [days]
	U.map.in.ssh_unitFactor = 100; % eg 100 if SSH data in cm, 1/10 if in deka m etc..
	U.map.in.fname='psvar.cdf';
	U.map.in.keys.lat='U_LAT_2D';
	U.map.in.keys.lon='U_LON_2D';
	U.map.in.keys.ssh='SSH';
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
	U.thresh.shape.iq=0.3; % isoperimetric quotient
	U.thresh.shape.chelt=0.3; % (diameter of circle with equal area)/(maximum distance between nodes) (if ~switch.IQ)
	U.thresh.corners=6; % min number of data points for the perimeter of an eddy
	U.thresh.dist=.2*24*60^2; % max distance travelled per day
	U.thresh.life=10; % min num of living days for saving
	U.thresh.ampArea=[.25 2.5]; % allowable factor between old and new time step for amplitude and area (1/4 and 5/1 ??? chelton)
	%% switches
	U.switchs.RossbyStuff=true;
	U.switchs.IQ=false;
	U.switchs.chelt=true;
	U.switchs.distlimit=true;
	U.switchs.AmpAreaCheck=true;
	U.switchs.netUstuff=true;
	%% parameters
	U.parameters.rossbySpeedFactor=1.75; % only relevant if cheltons method is used. eddy translation speed assumed factor*rossbyWavePhaseSpeed for tracking projections
	U.parameters.meanU=100; % depth from which to take mean U
	U.parameters.minProjecDist=150e3; % minimum linear_eccentricity*2 of ellipse (see chelton 2011)
	U.parameters.trackingRef='CenterOfVolume'; % choices: 'centroid', 'CenterOfVolume', 'Peak'
	%% technical params
	U.RossbyStuff.splits =12; % number of chunks for brunt väis calculations	
end
