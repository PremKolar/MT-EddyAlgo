function DD=input_vars
    %% threads / debug
    DD.threads.num=12;
    DD.debugmode=false;
      	    DD.debugmode=true;
    %% time
    DD.time.delta_t=3; % [days]!
    % 	DD.time.from.str='19940102'; % pop
    % % 	DD.time.till.str='20061231';
    %  	DD.time.till.str='19950101';
    DD.time.from.str=datestr(datenum('19840814','yyyymmdd'),'yyyymmdd');
    DD.time.till.str=datestr(datenum('19840814','yyyymmdd')+12*3*2,'yyyymmdd');
    %% dirs
    % 	DD.path.OutDirBaseName='AVso';
    [~,DD.path.OutDirBaseName]=fileparts(pwd);
    DD.path.TempSalt.name='../TempSalt/';
    DD.path.raw.name='../rawZ/';
    DD.path.UV.name='../rawZ/';
    %DD.path.TempSalt.name='/media/ROM/TempSalt/';
    %     DD.path.raw.name='/data/icdc/ocean/aviso_ssh/DATA/weekly/msla/';
    %     DD.path.UV.name='../TempSalt/';
    % 	DD.path.raw.name='/scratch/uni/ifmto/u241194/DAILY/EULERIAN/SSH/';
    %DD.path.raw.name='/media/ROM/SSH_POP/';
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
	DD.map.in.north=47;
%     %% output MAP STUFF
%     DD.map.out.X=20*1+1;
%     DD.map.out.Y=20* 1+1;
%     DD.map.out.west=-170;
%     DD.map.out.east=-150;
%     DD.map.out.south=-60;
%     DD.map.out.north=-40;
%     %% input MAP STUFF
%     DD.map.in.west=DD.map.out.west;
%     DD.map.in.east=DD.map.out.east;
%     DD.map.in.south=DD.map.out.south;
%     DD.map.in.north=DD.map.out.north;
   %    DD.map.in.ssh_unitFactor = 100; 
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
    % 	DD.map.in.fname='SSH_GLB_t.t0.1_42l_CORE.yyyymmdd.nc';
    % DD.map.in.fname='SsaltoDuacs__merged_msla__AVISO__ref__0.25deg__yyyymmdd.nc';
    %     DD.map.in.keys.lat='lat';
    % 	DD.map.in.keys.lon='lon';
    % 	DD.map.in.keys.ssh='msla';
    % 	DD.map.in.keys.time='time';
    % 	DD.map.in.keys.ssh='SSH';
    DD.TS.keys.lat='U_LAT_2D';
    DD.TS.keys.lon='U_LON_2D';
    DD.TS.keys.salt='SALT';
    DD.TS.keys.temp='TEMP';
    DD.TS.keys.depth='depth_t';
    %% thresholds
    DD.contour.step=0.01; % [SI]
    DD.thresh.ssh_filter_size=1;
    DD.thresh.radius=0; % [SI]
    DD.thresh.amp=0.01; % [SI]
    DD.thresh.shape.iq=0.3; % isoperimetric quotient
    DD.thresh.shape.chelt=0.3; % (diameter of circle with equal area)/(maximum distance between nodes) (if ~switch.IQ)
    DD.thresh.corners=8; % min number of data points for the perimeter of an eddy
    DD.thresh.dist=.1*24*60^2; % max distance travelled per day
    DD.thresh.life=3; % min num of living days for saving
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
    DD.parameters.meanU=-200; % depth from which to take mean U
    DD.parameters.meanUunit=1; % depth from which to take mean U
    DD.parameters.minProjecDist=150e3; % minimum linear_eccentricity*2 of ellipse (see chelton 2011)
DD.parameters.Gausswidth=1e5;   
 DD.parameters.trackingRef='CenterOfVolume'; % choices: 'centroid', 'CenterOfVolume', 'Peak'
    DD.parameters.Nknown=false; % Brunt-Väisälä f already in data
    DD.map.in.keys.N='N';
    %% technical params
    DD.RossbyStuff.splits =12; % number of chunks for brunt väis calculations
    DD.parameters.overrideWindowType=false;
    DD.parameters.overrideWindowType=true;
    %% relevant for S000
%    DD.parameters.overrideWindowType=false; % need to be true
    DD.parameters.overrideWindowType=true; % need to be true
    DD.parameters.SSHAdepth=-25;
    DD.parameters.boxlims.south=10;
    DD.parameters.boxlims.west=0;
    DD.map.in.cdfName='new2.cdf';
end










