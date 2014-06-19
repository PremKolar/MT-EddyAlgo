function DD=input_vars
	%% threads / debug
	DD.threads.num=12;
		DD.debugmode=false;
% 	DD.debugmode=true;
	%% time
	DD.time.delta_t=3; % [days]!
% 	DD.time.from.str='19091231';
	DD.time.from.str=datestr(now,'yyyymmdd');
	% 	 DD.time.from.str='19940425';
	%     DD.time.till.str='19960730';
	DD.time.till.str=datestr(now+424242,'yyyymmdd');
% 	DD.time.till.str='30000101';
	%% dirs
	DD.path.OutDirBaseName='madDemo';
		DD.path.TempSalt.name='../TempSalt/';
% 	    DD.path.TempSalt.name='~/ROMnew/TempSalt/';
	DD.path.raw.name='../demo/';
	%     DD.path.raw.name='/media/ROM/SSH_POP/';
	%% output MAP STUFF
	DD.map.out.X=46*1+1;
	DD.map.out.Y=37*1+1;
	DD.map.out.west=0;
	DD.map.out.east=46;
	DD.map.out.south=10;
	DD.map.out.north=47;
%     %% output MAP STUFF
% 	DD.map.out.X=10*1+1;
% 	DD.map.out.Y=10*1+1;
% 	DD.map.out.west=20;
% 	DD.map.out.east=30;
% 	DD.map.out.south=10;
% 	DD.map.out.north=20;
	%% input MAP STUFF
	DD.map.in.west=DD.map.out.west;
	DD.map.in.east=DD.map.out.east;
	DD.map.in.south=DD.map.out.south;
	DD.map.in.north=DD.map.out.north;
	DD.map.in.time.delta_t = 1; % [days]
	DD.map.in.ssh_unitFactor = 10; % eg 100 if SSH data in cm, 1/10 if in deka m etc..
	DD.map.in.fname='RAWyyyymmdd.nc';
	DD.map.in.cdfName='psvar.cdf';
	DD.map.in.keys.lat='U_LAT_2D';
	DD.map.in.keys.lon='U_LON_2D';
	DD.map.in.keys.ssh='SSH';
	DD.TS.keys.lat='U_LAT_2D';
	DD.TS.keys.lon='U_LON_2D';
	DD.TS.keys.salt='SALT';
	DD.TS.keys.temp='TEMP';
	DD.TS.keys.depth='depth_t';
	%% thresholds
	DD.contour.step=0.01; % [SI]
	DD.thresh.radius=0; % [SI]
	DD.thresh.amp=0.01; % [SI]
	DD.thresh.shape.iq=0.1; % isoperimetric quotient
	DD.thresh.shape.chelt=0.3; % (diameter of circle with equal area)/(maximum distance between nodes) (if ~switch.IQ)
	DD.thresh.corners=6; % min number of data points for the perimeter of an eddy
	DD.thresh.dist=.5*24*60^2; % max distance travelled per day
	DD.thresh.life=9; % min num of living days for saving TODO check timestep or day
	DD.thresh.ampArea=[.25 2.5]; % allowable factor between old and new time step for amplitude and area (1/4 and 5/1 ??? chelton)
	%% switches
	DD.switchs.RossbyStuff=false;
	DD.switchs.IQ=true;
	DD.switchs.chelt=false;
	DD.switchs.distlimit=false;
	DD.switchs.AmpAreaCheck=false;
	DD.switchs.netUstuff=false;
	%% parameters
	DD.parameters.rossbySpeedFactor=1.75; % only relevant if cheltons method is used. eddy translation speed assumed factor*rossbyWavePhaseSpeed for tracking projections
	DD.parameters.meanU=100; % depth from which to take mean U
	DD.parameters.minProjecDist=150e3; % minimum linear_eccentricity*2 of ellipse (see chelton 2011)
	DD.parameters.trackingRef='CenterOfVolume'; % choices: 'centroid', 'CenterOfVolume', 'Peak'
	%% technical params
	DD.RossbyStuff.splits =12; % number of chunks for brunt väis calculations
	%% only relevant for S000
	DD.parameters.boxlims.south=10;
	DD.parameters.boxlims.west=0;
end

















