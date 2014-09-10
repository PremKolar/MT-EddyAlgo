
function DD=INPUT
    DD.template='pop';
    %% threads / debug
    DD.threads.num=12;
    DD.debugmode=false;
    %     DD.debugmode=true;
     DD.overwrite=false;
%     DD.overwrite=true;
   
     %% time
    DD.time.from.str='19940102';
    DD.time.till.str='19990102';
%      threshlife=20*7
    
%     DD.time.till.str='19950702';    
    threshlife=10*7;
 
%     %% window on globe
%     DD.map.in.west=-70;
%     DD.map.in.east= -30;
%     DD.map.in.south= 5;
%     DD.map.in.north= 45;
%     %% output map res
%     DD.map.out.X=40*1+1; % TODO
%     DD.map.out.Y=40*1+1;   
 %% window on globe
    DD.map.in.west=-70;
    DD.map.in.east= -20;
    DD.map.in.south= 0;
    DD.map.in.north= 60;
    %% output map res
    DD.map.out.X=50*1+1; % TODO
    DD.map.out.Y=60*1+1;   