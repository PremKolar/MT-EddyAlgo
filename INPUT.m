% templates :
% 'udef' - user defined in INPUTuserDef.m
% 'pop' - template for POP SSH data
% 'aviso' - template for AVISO SSH data
% 'mad' - template for Madeleine's data
function DD=INPUT
    DD.template='aviso';
    %% threads / debug
    DD.threads.num=12;
    DD.debugmode=false;
%      DD.debugmode=true;
    DD.overwrite=false;
%        DD.overwrite=true;
    %% time
    DD.time.from.str='19940105';
%   DD.time.till.str='19991231';
     DD.time.till.str='20061231';
    %      threshlife=20*7
    threshlife=7*8;
    %% window on globe
    DD.map.in.west=-180;
    DD.map.in.east= 180;
    DD.map.in.south= -80;
    DD.map.in.north= 80;
    %% output map res
    DD.map.out.X=360*1+1; % TODO
    DD.map.out.Y=160*1+1;
    %% thresholds
    DD.contour.step=0.01; % [SI]
    DD.thresh.radius=0; % [SI]
    DD.thresh.maxRadiusOverRossbyL=4; %!
    DD.thresh.amp=0.01; % [SI]
    DD.thresh.shape.iq=0.55; % isoperimetric quotient
    DD.thresh.shape.chelt=0.2; % (diameter of circle with equal area)/(maximum distance between nodes) (if ~switch.IQ)
    DD.thresh.corners.min=8; % min number of data points for the perimeter of an eddy
    DD.thresh.corners.max=1e42; % dangerous..
    DD.thresh.life=threshlife; % min num of living days for saving
    DD.thresh.ampArea=[.25 2.5]; % allowable factor between old and new time step for amplitude and area (1/4 and 5/1 ??? chelton)
    DD.thresh.IdentityCheck=[2.5]; % 1: perfect fit, 2: 100% change ie factor 2 in either sigma or amp
    %% switches
    DD.switchs.IQ=1;
    DD.switchs.chelt=0;
    DD.switchs.RossbyStuff=1;
    DD.switchs.distlimit=1;
    DD.switchs.AmpAreaCheck=0;
    DD.switchs.netUstuff=0;
    DD.switchs.meanUviaOW=0;
    DD.switchs.IdentityCheck=1;
    DD.switchs.maxRadiusOverRossbyL=1;
    DD.switchs.spaciallyFilterSSH=0;
    DD.switchs.filterSSHinTime=1;
end

