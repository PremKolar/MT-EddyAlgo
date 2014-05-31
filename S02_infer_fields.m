%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 04-Apr-2014 16:53:06
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates geostrophic data from SSH
function S02_infer_fields
    %% init
    DD=initialise('conts');
    %% read input file
    cut1=load([DD.path.cuts.name DD.path.cuts.files(1).name]);
    DD.coriolis=coriolisStuff(cut1.grids);
    DD.threads.num=init_threads(DD.threads.num);
    RS=getRossbyStuff(DD);
    %% spmd
    main(DD,RS)
    %% save info file
    save_info(DD)
end






function main(DD,RS)
    if DD.debugmode
        spmd_body(DD,RS);
    else
        spmd(DD.threads.num)
            spmd_body(DD,RS);
        end
    end
end

function RS=getRossbyStuff(DD)
    if DD.switchs.RossbyStuff
        file=[DD.path.Rossby.name DD.path.Rossby.files.name];
        RS.c=nc_varget(file,'RossbyPhaseSpeed');
       RS.Lr=nc_varget(file,'RossbyPhaseSpeed');    
    else
        RS=[];
    end
end


function spmd_body(DD,RS)
    %% loop
    JJ=DD.threads.lims(labindex,1):DD.threads.lims(labindex,2);
    T=disp_progress('init','infering fields');
    for jj=JJ
        T=disp_progress('disp',T,numel(JJ),100);
        %% load
        cut=LoadCut(jj,DD);      
   coriolis=coriolisStuff(cut.grids);
        %% calc
        grids=geostrophy(cut.grids,coriolis,RS);
        
        %% write
        write_fields(DD,jj,'cuts',grids);
        write_fields(DD,jj,'conts',grids);
    end
end
function gr=geostrophy(gr,corio,RS)
    %% SSH gradient
    [gr.SSHgrad_x,gr.SSHgrad_y]=dSSHdxi(gr.SSH,gr.DX,gr.DY);
    %% velocities
    gr.U=-corio.GOverF.*gr.SSHgrad_y;
    gr.V= corio.GOverF.*gr.SSHgrad_x;
    gr.absUV=hypot(abs(gr.U),abs(gr.V));
    %% deformation
    def=deformation(gr);
    gr.vorticity      = def.dVdx - def.dUdy;
    gr.divergence= def.dUdx + def.dVdy;
    gr.stretch   = def.dUdx - def.dVdy;
    gr.shear     = def.dVdx + def.dUdy;
    %% okubo weiss
    gr.OW=.5*(-gr.vorticity.*2+gr.divergence.*2+gr.stretch.*2+gr.shear.*2);
    %% assuming Ro=1
    gr.L=gr.absUV./corio.f;
    kinVis=1e-6;
    gr.Re=gr.absUV.*gr.L/kinVis;
    gr.Ro=ones(size(gr.L));
    gr.Rrhines=earthRadius./gr.L;
    gr.Lrhines=sqrt(gr.absUV./corio.beta);
    gr.L_R=abs(RS.c./corio.f);
    gr.Bu=(gr.L_R./gr.L).^2;
end
function def=deformation(grids)
    %% calc U gradients
    def.dUdy=[ nan(1,size(grids.U,2));diff(grids.U,1,1)] ./ grids.DY;
    def.dVdx=[ nan(size(grids.V,1),1), diff(grids.V,1,2)] ./ grids.DX;
    def.dVdy=[ nan(1,size(grids.V,2));diff(grids.V,1,1)] ./ grids.DY;
    def.dUdx=[ nan(size(grids.U,1),1), diff(grids.U,1,2)] ./ grids.DX;
end
function cut=LoadCut(jj,dd)
    file_f=dd.path.cuts.files(jj).name;
    file_b=dd.path.cuts.name;
    file =[file_b, file_f ];
    cut=load(file);
end
function [dSSHdx,dSSHdy]=dSSHdxi(SSH,DX,DY)
    %% calc SSH gradients
    dSSHdx=[diff(SSH,1,2), nan(size(SSH,1),1)] ./ DX;
    dSSHdy=[diff(SSH,1,1); nan(1,size(SSH,2))] ./ DY;
end
function out=coriolisStuff(fields)
    %% omega
    out.Omega=angularFreqEarth;
    %% f
    out.f=2*out.Omega*sind(fields.LAT);
    %% beta
    out.beta=2*out.Omega/earthRadius*cosd(fields.LAT);
    %% gravity
    out.g=sw_g(fields.LAT,zeros(size(fields.LAT)));
    %% g/f
    out.GOverF=out.g./out.f;
end
function om=angularFreqEarth
    T=day_sid;
    om=2.0*pi/T;
    function d=day_sid
        d=23.9344696*60*60; % wikipedia
    end
end
