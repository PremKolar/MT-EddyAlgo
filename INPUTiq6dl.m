
function DD=INPUT
    DD.template='pop';
    %% threads / debug
    DD.threads.num=12;
    DD.debugmode=false;
%     DD.debugmode=true;
    DD.overwrite=true;
     %% time
    DD.time.from.str='19940102';
    %     DD.time.till.str='19940602';
    DD.time.till.str='19990102';
    %% window on globe
    DD.map.in.west=-70;
    DD.map.in.east= -40;
    DD.map.in.south= 30;
    DD.map.in.north= 40;
    %% output map res
    DD.map.out.X=30*1+1; % TODO
    DD.map.out.Y=10*1+1;
    %% thresholds
    DD.contour.step=0.01; % [SI]
    DD.thresh.radius=0; % [SI]
    DD.thresh.maxRadiusOverRossbyL=10; % [SI]   %% GOOD???
    DD.thresh.amp=0.01; % [SI]
    DD.thresh.shape.iq=0.6; % isoperimetric quotient
    DD.thresh.shape.chelt=0.2; % (diameter of circle with equal area)/(maximum distance between nodes) (if ~switch.IQ)
    DD.thresh.corners.min=12; % min number of data points for the perimeter of an eddy
    DD.thresh.corners.max=pi*2e6*1e-4; % at dx ~1e-4 -> skip eddies(radius> ~1000km) , just for performance
    DD.thresh.life=20*7; % min num of living days for saving
    DD.thresh.ampArea=[.25 2.5]; % allowable factor between old and new time step for amplitude and area (1/4 and 5/1 ??? chelton)
    DD.thresh.IdentityCheck=[1/2.5 2.5];
    %% switches
    DD.switchs.IQ=true;
    DD.switchs.chelt=false;
    DD.switchs.RossbyStuff=true;
    DD.switchs.distlimit=true;
    DD.switchs.AmpAreaCheck=false;
    DD.switchs.netUstuff=false;
    DD.switchs.meanUviaOW=false;
    DD.switchs.IdentityCheck=false;
    DD.switchs.maxRadiusOverRossbyL=false;
    DD.switchs.spaciallyFilterSSH=false;
    DD.switchs.filterSSHinTime=true;
end