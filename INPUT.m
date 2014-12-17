% templates :
% 'udef' - user defined in INPUTuserDef.m
% 'pop' - template for POP SSH data
% 'aviso' - template for AVISO SSH data
% 'mad' - template for Madeleine's data
% 'pop2avi' -
function DD=INPUT
    DD.template='aviso';
    %% threads / debug
    DD.threads.num=12;
    DD.debugmode=false;
%     DD.debugmode=true;
    DD.overwrite=false;
%         DD.overwrite=true;
    %% time
    DD.time.from.str  ='19940105'; %first pop/avi
    DD.time.till.str  ='20061227'; % last pop/avi
%     f='yyyymmdd';
%     dateplus=@(D,a,f) datestr(datenum(D,f)+a,f);
%     DD.time.till.str  = dateplus(DD.time.from.str,1*365,f);
%     DD.time.till.str='19950105';
%     threshlife=20*7
DD.time.delta_t=7; % [days]!   
    threshlife=7*8; % TODO
    %% window on globe (0:360° system)
    DD.map.in.west= 0;
    DD.map.in.east= 360;
    DD.map.in.south= -80;
    DD.map.in.north= 80;
    %% thresholds
    DD.contour.step=0.01; % [SI]
    DD.thresh.radius=0; % [SI]
    DD.thresh.maxRadiusOverRossbyL=4; %[ ]
    DD.thresh.minRossbyRadius=20e3; %[SI]
    DD.thresh.amp=DD.contour.step; % [SI]
    DD.thresh.shape.iq=0.55; % isoperimetric quotient [ ]
    DD.thresh.corners.min=8; % min number of data points for the perimeter of an eddy[ ]
    DD.thresh.corners.max=1e42; % dangerous.. [ ]
    DD.thresh.life=threshlife; % min num of living days for saving [days]
    DD.thresh.ampArea=[.25 2.5]; % allowable factor between old and new time step for amplitude and area (1/4 and 5/2 ??? chelton)
    DD.thresh.IdentityCheck=[2]; % 1: perfect fit, 2: 100% change ie factor 2 in either sigma or amp
    DD.thresh.phase = 0.2; % max(abs(rossby phase speed)) [SI]
    %% switches
    DD.switchs.IQ=1;
    DD.switchs.chelt=0;
    DD.switchs.RossbyStuff=1;  % TODO no choice
    DD.switchs.distlimit=1;      % TODO no choice
    DD.switchs.AmpAreaCheck=0;
    DD.switchs.netUstuff=0;
    DD.switchs.meanUviaOW=0;
    DD.switchs.IdentityCheck=1;
    DD.switchs.maxRadiusOverRossbyL=1;  % TODO no choice  
    DD.switchs.spaciallyFilterSSH=0;  % TODO delete
    DD.switchs.filterSSHinTime=1;
    %%
    DD.parameters.fourierOrder=4;
end

