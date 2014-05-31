function U=input_vars
%% threads
U.threads.num=12;
U.debugmode=0;
%% time
  U.time.from.str='19921014';
% U.time.till.str='20130807';
   U.time.till.str='19931212';
  U.time.delta_t=7; % [days]!
%% dirs
U.path.TempSalt.name='../TempSalt/';
  U.path.raw.name='/data/icdc/ocean/aviso_ssh/DATA/weekly/msla/';
   U.path.root='../dataTIN/';
    U.path.plots='../plotTIN/';
  %% thresholds
U.contour.step=0.01; % [SI]
U.thresh.ssh_filter_size=1;
U.thresh.radius=0; % [SI]
U.thresh.amp=0.01; % [SI]
U.thresh.shape.iq=0.3; % isoperimetric quotient
U.thresh.shape.chelt=0.3; % (diameter of circle with equal area)/(maximum distance between nodes) (if ~switch.IQ)
U.thresh.corners=6; % min number of data points for the perimeter of an eddy
U.thresh.dist=.5*24*60^2; % max distance travelled per day
U.thresh.life=20; % min num of living days for saving
 U.thresh.ampArea=[.25 2.5]; % allowable factor between old and new time step for amplitude and area (1/4 and 5/1 ??? chelton)
%% dims for map plots
U.dim.X=80*1+1;
  U.dim.Y=50*1+1;
  U.dim.west=-180;
  U.dim.east=180;
  U.dim.south=-50;
  U.dim.north=-10;
    %% switches
    U.switchs.RossbyStuff=true; % TODO
    U.switchs.IQ=false;
    U.switchs.chelt=true;
U.switchs.distlimit=true;
U.switchs.AmpAreaCheck=true;
    %% parameters
U.parameters.rossbySpeedFactor=1.75; % only relevant if cheltons method is used. eddy translation speed assumed factor*rossbyWavePhaseSpeed for tracking projections
    U.parameters.depthRossby=4242424242; % depth for which to integrate rossby phase speed and radius
    U.parameters.meanU=100; % depth from which to take mean U
U.parameters.minProjecDist=150e3; % minimum linear_eccentricity*2 of ellipse (see chelton 2011)
U.parameters.trackingRef='centroid'; % choices: 'centroid', 'CenterOfVolume', 'Peak'
    %% technical params
    U.RossbyStuff.splits = 24; % number of chunks for brunt väis calculations
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