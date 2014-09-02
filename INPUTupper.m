
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
%     DD.time.till.str='19990102';
%      threshlife=20*7
    
    DD.time.till.str='19940702';    
    threshlife=5*7;
 
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
    DD.map.in.east= -40;
    DD.map.in.south= 10;
    DD.map.in.north= 40;
    %% output map res
    DD.map.out.X=30*1+1; % TODO
    DD.map.out.Y=30*1+1;   