function U=input_vars
%% threads
U.threads.num=12;
U.debugmode=1;
%% time
   U.time.from.str='19940102';
%    U.time.till.str='20061231';
  U.time.from.str='19940502';
  U.time.delta_t=1; % [days]!
%% dirs
 U.path.TempSalt.name='../TempSalt/';
% U.path.TempSalt.name='/media/ROM/TempSalt/';
U.path.raw.name='/scratch/uni/ifmto/u241194/DAILY/EULERIAN/SSH/';
% U.path.raw.name='/media/ROM/SSH_POP/';
   U.path.root='../dataWrld/';
    U.path.plots='../plotWrld/';
  %% thresholds
U.contour.step=0.1; % [SI]
U.thresh.ssh_filter_size=1;
U.thresh.radius=0; % [SI]
U.thresh.amp=0.1; % [SI]
U.thresh.shape.iq=0.3; % isoperimetric quotient
U.thresh.shape.chelt=0.3; % (diameter of circle with equal area)/(maximum distance between nodes) (if ~switch.IQ)
U.thresh.corners=6; % min number of data points for the perimeter of an eddy
U.thresh.dist=.5*24*60^2; % max distance travelled per day
U.thresh.life=20; % min num of living days for saving
 U.thresh.ampArea=[.25 2.5]; % allowable factor between old and new time step for amplitude and area (1/4 and 5/1 ??? chelton)
%% dims for map plots
U.dim.X=360*1+1;
  U.dim.Y=180*1+1;
  U.dim.west=-180;
  U.dim.east=-180;
  U.dim.south=-90;
  U.dim.north=90;
    %% switches
    U.switchs.RossbyStuff=true; % TODO
    U.switchs.IQ=true;
    U.switchs.chelt=false;
U.switchs.distlimit=false;
U.switchs.AmpAreaCheck=false;
    %% parameters
U.parameters.rossbySpeedFactor=1.75; % only relevant if cheltons method is used. eddy translation speed assumed factor*rossbyWavePhaseSpeed for tracking projections
    U.parameters.depthRossby=4242424242; % depth for which to integrate rossby phase speed and radius
    U.parameters.meanU=100; % depth from which to take mean U
U.parameters.minProjecDist=150e3; % minimum linear_eccentricity*2 of ellipse (see chelton 2011)
U.parameters.trackingRef='centroid'; % choices: 'centroid', 'CenterOfVolume', 'Peak'
    %% technical params
    U.RossbyStuff.splits = 4; % number of chunks for brunt väis calculations
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
