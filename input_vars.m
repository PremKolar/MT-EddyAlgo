function DD=input_vars
	%% threads / debug
	DD.threads.num=12;
   	DD.debugmode=false;
% 	DD.debugmode=true;
	%% time
	DD.time.delta_t=3; % [days]!
% 	DD.time.from.str='19091231';
	DD.time.from.str=datestr(datenum('19840814','yyyymmdd'),'yyyymmdd');
	% 	 DD.time.from.str='19940425';
	%     DD.time.till.str='19960730';
	DD.time.till.str=datestr(datenum('19840814','yyyymmdd')+12*3000,'yyyymmdd');
% 	DD.time.till.str='30000101';
	%% dirs
	DD.path.OutDirBaseName='u1';
		DD.path.TempSalt.name='../TempSalt/';
% 	    DD.path.TempSalt.name='~/ROMnew/TempSalt/';
	DD.path.raw.name='../rawZ/';
	%     DD.path.raw.name='/media/ROM/SSH_POP/';
    %% output MAP STUFF
	DD.map.out.X=54*1+1;
	DD.map.out.Y=47*1+1;
	DD.map.out.west=0; 
	DD.map.out.east=54;
	DD.map.out.south=10;
	DD.map.out.north=47;
    %% input MAP STUFF
   	DD.map.in.west=-20; 
	DD.map.in.east=70;
	DD.map.in.south=10;
	DD.map.in.north=46.8;
	DD.map.in.ssh_unitFactor = 1000; % eg 100 if SSH data in cm, 1/10 if in deka m etc..
	DD.map.in.fname='RAWyyyymmdd.nc';	
	DD.map.in.keys.lat='U_LAT_2D';
	DD.map.in.keys.lon='U_LON_2D';
	DD.map.in.keys.ssh='SSHA';	
	DD.map.in.keys.x='XT';
	DD.map.in.keys.y='YT';
	DD.map.in.keys.z='ZT';
	DD.map.in.keys.time='TIME';
	DD.map.in.keys.U='U';
	DD.map.in.keys.V='V';
	DD.TS.keys.lat='U_LAT_2D';
	DD.TS.keys.lon='U_LON_2D';
	DD.TS.keys.salt='SALT';
	DD.TS.keys.temp='TEMP';
	DD.TS.keys.depth='depth_t';
	%% thresholds
	DD.contour.step=0.01; % [SI]
	DD.thresh.radius=0; % [SI]
	DD.thresh.amp=0.01; % [SI]
	DD.thresh.shape.iq=0.5; % isoperimetric quotient
	DD.thresh.shape.chelt=0.3; % (diameter of circle with equal area)/(maximum distance between nodes) (if ~switch.IQ)
	DD.thresh.corners=6; % min number of data points for the perimeter of an eddy
	DD.thresh.dist=.5*24*60^2; % max distance travelled per day
	DD.thresh.life=10; % min num of living days for saving
	DD.thresh.ampArea=[.25 2.5]; % allowable factor between old and new time step for amplitude and area (1/4 and 5/1 ??? chelton)
	%% switches
	DD.switchs.RossbyStuff=false;
	DD.switchs.IQ=true;
	DD.switchs.chelt=false;
	DD.switchs.distlimit=false;
	DD.switchs.AmpAreaCheck=false;
	DD.switchs.netUstuff=true;
	%% parameters
	DD.parameters.rossbySpeedFactor=1.75; % only relevant if cheltons method is used. eddy translation speed assumed factor*rossbyWavePhaseSpeed for tracking projections
	DD.parameters.meanU=100; % depth from which to take mean U
	DD.parameters.minProjecDist=150e3; % minimum linear_eccentricity*2 of ellipse (see chelton 2011)
	DD.parameters.trackingRef='CenterOfVolume'; % choices: 'centroid', 'CenterOfVolume', 'Peak'
   DD.parameters.Nknown=true; % Brunt-Väisälä f already in data
	DD.map.in.keys.N='N';
	%% technical params
	DD.RossbyStuff.splits =12; % number of chunks for brunt väis calculations
	%% only relevant for S000
   DD.parameters.SSHAdepth=-25;
	DD.parameters.boxlims.south=10;
	DD.parameters.boxlims.west=0;
    DD.parameters.overrideWindowType=true;
	DD.map.in.cdfName='new2.cdf';
end

















