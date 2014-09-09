
function DD=INPUT
    DD.template='pop';
    %% threads / debug
    DD.threads.num=8;
    DD.debugmode=false;
%         DD.debugmode=true;
     DD.overwrite=false;
%     DD.overwrite=true;
   
     %% time
    DD.time.from.str='19940102';
%     DD.time.till.str='19990102';
%      threshlife=20*7
    
    DD.time.till.str='20500101';    
    threshlife=8*3;
 
%     %% window on globe
%     DD.map.in.west=-70;
%     DD.map.in.east= -30;
%     DD.map.in.south= 5;
%     DD.map.in.north= 45;
%     %% output map res
%     DD.map.out.X=40*1+1; % TODO
%     DD.map.out.Y=40*1+1;   
 %% window on globe
    DD.map.in.west=-180;
    DD.map.in.east= 180;
    DD.map.in.south= -80;
    DD.map.in.north= 80;
    %% output map res
    DD.map.out.X=360*1+1; % TODO
    DD.map.out.Y=160*1+1;   %% thresholds
    DD.contour.step=0.01; % [SI]
    DD.thresh.radius=0; % [SI]
    DD.thresh.maxRadiusOverRossbyL=5; %down to 4 next time!!!!!!!!!!!!!!!!!!!!!!!% [SI]   %% GOOD???
    DD.thresh.amp=0.01; % [SI]
    DD.thresh.shape.iq=0.6; % isoperimetric quotient
    DD.thresh.shape.chelt=0.2; % (diameter of circle with equal area)/(maximum distance between nodes) (if ~switch.IQ)
    DD.thresh.corners.min=8; % min number of data points for the perimeter of an eddy
    DD.thresh.corners.max=pi*2e6*1e-4; % at dx ~1e-4 -> skip eddies(radius> ~1000km) , just for performance
    DD.thresh.life=threshlife; % min num of living days for saving
     DD.thresh.ampArea=[.25 2.5]; % allowable factor between old and new time step for amplitude and area (1/4 and 5/1 ??? chelton)
    DD.thresh.IdentityCheck=[2];
    %% switches
    DD.switchs.IQ=true;
    DD.switchs.chelt=false;
    DD.switchs.RossbyStuff=true;
    DD.switchs.distlimit=true;
    DD.switchs.AmpAreaCheck=false;
    DD.switchs.netUstuff=false;
    DD.switchs.meanUviaOW=false;
    DD.switchs.IdentityCheck=false;
    DD.switchs.maxRadiusOverRossbyL=true;
    DD.switchs.spaciallyFilterSSH=false;
    DD.switchs.filterSSHinTime=true;
end
