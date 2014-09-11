%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 04-Apr-2014 16:53:06
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% walks through all the contours and decides whether they qualify
function s05plotsAnima
    %% init
    DD=initialise('conts',mfilename);
    DD.threads.num=init_threads(DD.threads.num);
    %% spmd
    main(DD);
    %% update infofile
    conclude(DD);
    system(['mencoder "mf://*.jpeg" -mf fps=20  -o flat.avi -ovc lavc -lavcopts  vcodec=ljpeg'])
    system(['rm ./*.jpeg'])
    system(['mplayer flat.avi'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main(DD)
    if DD.debugmode
        spmd_body(DD)
    else
        spmd(DD.threads.num)
            spmd_body(DD)
            disp_progress('conclude');
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spmd_body(DD)
    [JJ]=SetThreadVar(DD);
    Td=disp_progress('init','making jpegs for movie');
    for jj=1:numel(JJ)
        work_day(DD,JJ(jj));
        Td=disp_progress('disp',Td,numel(JJ),numel(JJ));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [EE,skip]=work_day(DD,JJ)
    %% check for exisiting data
    skip=false;
    EE.filename.cont=JJ.files;
    EE.filename.cut=[DD.path.cuts.name, DD.pattern.prefix.cuts, JJ.protos];
    EE.filename.self=[DD.path.eddies.name, DD.pattern.prefix.eddies ,JJ.protos];
    %          if exist(EE.filename.self,'file'), skip=true; return; end
    %% plot
    makejpegs(EE,JJ.daynums);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makejpegs(EE,dayn)
    load(EE.filename.self);
    load(EE.filename.cut);
    figure(labindex);clf;
    %     CM=flipud(hot(100));
    contour(grids.ssh,-.6:.03:.6,'linewidth',2)
    for kk=1:numel(anticyclones)
        x=anticyclones(kk).coordinates.exact.x;
        y=anticyclones(kk).coordinates.exact.y;
        iq=anticyclones(kk).isoper-.55;
        iq=iq/0.45*100;
        iq(iq<1)=1;iq(iq>100)=100;
        iq=round(iq);
        hold on
        plot(x,y,'color','black','linewidth',ceil(iq*5/100))
    end
    
    for kk=1:numel(cyclones)
        x=cyclones(kk).coordinates.exact.x;
        y=cyclones(kk).coordinates.exact.y;
        iq=cyclones(kk).isoper-.55;
        iq=iq/0.45*100;
        iq(iq<1)=1;iq(iq>100)=100;
        iq=round(iq);
        hold on
        plot(x,y,'--','color','black','linewidth',ceil(iq*5/100))
    end
    axis tight off;
    title(num2str(dayn));
    savefig2png4mov('./',90,1024,768,datestr(dayn,'yymmdd'));
end

