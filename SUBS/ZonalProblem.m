function [OUT]=ZonalProblem(IN,window)
    %% shorthands
    Wflag=window.flag;
    Wlin=window.limits;
    %% full globe?
    full_globe.x=false;
    full_globe.y=false;
    if (Wlin.east-Wlin.west+1==size(Wflag,2))
        full_globe.x=true;
        if (Wlin.north-Wlin.south+1==size(Wflag,1))
            full_globe.y=true;
        end
        OUT=AppendIfFullZonal(IN,window);% longitude edge crossing has to be addressed
    end
    %% seam crossing?
    seam=false;
    if ~full_globe.x
        if (Wlin.west>Wlin.east) % ie not full globe but both seam ends are within desired window
            seam=true; % piece crosses long seam
            [OUT,window]=SeamCross(IN,window);
        else % desired piece is within global fields, not need for stitching
            OUT=AllGood(IN,Wlin);
        end
    end
    %% append params
    OUT.window	 =window;
    OUT.params.full_globe =full_globe;
    OUT.params.seam   =seam;
end
function [OUT]=AppendIfFullZonal(IN,window)
    %% append 1/10 of map to include eddies on seam
    % S04_track_eddies is able to avoid counting 1 eddy twice
    ss=window.limits.south;
    nn=window.limits.north;
    xadd=round(window.size.X/10);
    fields=fieldnames(IN);
    for field=fields';ff=field{1};
        %% init
        OUT.grids.(ff)=IN.(ff)(ss:nn,:);
        %% append
        OUT.grids.(ff)=OUT.grids.(ff)(:,[1:end, 1:xadd]);
    end
end
function [OUT,window]=SeamCross(IN,window)
    Wflag=window.flag;
    Wlin=window.limits;
    %% find new west and east
    easti =find(sum(double(Wflag))==0,1,'first');
    westi =find(sum(double(Wflag))==0,1,'last');
    southi=Wlin.south;
    northi=Wlin.north;
    %% reset east and west
    window.limits.west=westi;
    window.limits.east=easti;
    %% stitch 2 pieces 2g4
    fields=fieldnames(IN);
    for field=fields';ff=field{1};
        OUT.grids.(ff) =[IN.(ff)(southi:northi,westi:end) IN.(ff)(southi:northi,1:easti)];
    end
end
function OUT=AllGood(IN,Wlin)
    %% cut piece
    OUT.window.limits=Wlin;
    fields=fieldnames(IN);
    for field=fields';ff=field{1};
        OUT.grids.(ff) =IN.(ff)(Wlin.south:Wlin.north,Wlin.west:Wlin.east);
    end
end

