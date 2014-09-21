% templates :
% 'udef' - user defined in INPUTuserDef.m
% 'pop' - template for POP SSH data
% 'aviso' - template for AVISO SSH data
% 'mad' - template for Madeleine's data

	function DD=INPUT
	    DD.template='pop';
	    %% threads / debug
	    DD.threads.num=12;
	    DD.debugmode=false;
% 	    DD.debugmode=true;	
	    DD.overwrite=false;
% 	     DD.overwrite=true;
	    %% time
	
	    DD.time.from.str='19940202';
	    %DD.time.till.str=datestr(datenum(DD.time.from.str,'yyyymmdd')+10*3,'yyyymmdd');
	  DD.time.till.str='19940505';
	    threshlife=3*7;
	    %% window on globe
	    DD.map.in.west=-60;
	    DD.map.in.east= -50;
	    DD.map.in.south= 30;
	    DD.map.in.north= 40;
	    %% output map res
	    DD.map.out.X=10*1+1; % TODO
	    DD.map.out.Y=10*1+1;
	    %% thresholds
	    DD.contour.step=0.01; % [SI]
	    DD.thresh.radius=0; % [SI]
	    DD.thresh.maxRadiusOverRossbyL=4; %!
	    DD.thresh.amp=0.01; % [SI]
	    DD.thresh.shape.iq=0.3; % isoperimetric quotient
	    DD.thresh.shape.chelt=0.2; % (diameter of circle with equal area)/(maximum distance between nodes) (if ~switch.IQ)
	    DD.thresh.corners.min=12; % min number of data points for the perimeter of an eddy
	     %% DANGEROUS !
	    DD.thresh.corners.max=1e42; % at dx ~1e-4 -> skip eddies(radius> ~1000km) , just for performance
	    DD.thresh.life=threshlife; % min num of living days for saving
	    DD.thresh.ampArea=[.25 2.5]; % allowable factor between old and new time step for amplitude and area (1/4 and 5/1 ??? chelton)
	    DD.thresh.IdentityCheck=[2];
	    %% switches
	    DD.switchs.RossbyStuff=1;
	    DD.switchs.distlimit=1;
	    DD.switchs.netUstuff=0;
	    DD.switchs.meanUviaOW=0;
	    DD.switchs.spaciallyFilterSSH=0;
	    DD.switchs.filterSSHinTime=1;
	    %%
	    DD.switchs.IQ=1;
	    DD.switchs.IdentityCheck=1;
	    DD.switchs.maxRadiusOverRossbyL=1;
	    %%
	    DD.switchs.chelt=0;
	    DD.switchs.AmpAreaCheck=0;
	end
	
