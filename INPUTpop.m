function DD=INPUTpop
    %% time step
    DD.time.delta_t=1; % [days]!
    %% dirs
    [~,DD.path.OutDirBaseName]=fileparts(pwd);
    DD.path.TempSalt.name='../TempSaltUV/';
    DD.path.UV.name='../TempSaltUV/';
    DD.path.raw.name='/scratch/uni/ifmto/u241194/DAILY/EULERIAN/SSH/';
    %% map in keys
    DD.map.in.fname='SSH_GLB_t.t0.1_42l_CORE.yyyymmdd.nc';
    DD.map.in.keys.lat='U_LAT_2D';
    DD.map.in.keys.lon='U_LON_2D';
    DD.map.in.keys.ssh='SSH';
    DD.map.in.keys.time='TIME';
    DD.map.in.keys.U='U';
    DD.map.in.keys.V='V';
    %{...
    DD.map.in.keys.x='XT';
    DD.map.in.keys.y='YT';
    DD.map.in.keys.z='ZT';
    DD.map.in.keys.N='N';
    DD.map.in.cdfName='new2.cdf';
    %...}
    %% temp salt keys
    DD.TS.keys.lat='U_LAT_2D';
    DD.TS.keys.lon='U_LON_2D';
    DD.TS.keys.salt='SALT';
    DD.TS.keys.temp='TEMP';
    DD.TS.keys.depth='depth_t';
    %% parameters
    DD.parameters.ssh_unitFactor = 100; % eg 100 if SSH data in cm, 1/10 if in deka m etc..
    DD.parameters.rossbySpeedFactor=1.75; % only relevant if cheltons method is used. eddy translation speed assumed factor*rossbyWavePhaseSpeed for tracking projections
    DD.parameters.meanU=100; % depth from which to take mean U
    DD.parameters.meanUunit=1; % depth from which to take mean U
    DD.parameters.minProjecDist=150e3; % minimum linear_eccentricity*2 of ellipse (see chelton 2011)
    DD.parameters.Gausswidth=1e5;
    DD.parameters.trackingRef='CenterOfVolume'; % choices: 'centroid', 'CenterOfVolume', 'Peak'
    DD.parameters.Nknown=false; % Brunt-Väisälä f already in data
    DD.parameters.forceZonalInf=false;
    DD.parameters.RossbySplits =12; % number of chunks for brunt väis calculations
    %{...
    DD.parameters.SSHAdepth=-25;
    %...}
     %%
     DD.switches.rehashDD=true;
    
    
    
    
    
    
    
    