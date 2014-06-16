%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 04-Apr-2014 16:53:06
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates geostrophic data from ssh
function S02_infer_fields
    %% init
    DD=initialise('cuts');
    %% read input file
   cut1=load( DD.checks.passed(1).filenames);
    DD.coriolis=coriolisStuff(cut1.grids);
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
        RS.Lr=nc_varget(file,'RossbyRadius');
    else
        RS=[];
    end
end

function spmd_body(DD,RS)
    %% distro chunks to threads
    [JJ]=SetThreadVar(DD);
    T=disp_progress('init','infering fields');
    for jj=1:numel(JJ)
        T=disp_progress('disp',T,numel(JJ),100);
        %% load
        cut=load(JJ(jj).files);
        coriolis=coriolisStuff(cut.grids);
        %% calc
        grids=geostrophy(cut.grids,coriolis,RS); %#ok<NASGU>
        %% write
        save(JJ(jj).files,'grids','-append');
    end
end
function gr=geostrophy(gr,corio,RS)
    %% ssh gradient
    [gr.sshgrad_x,gr.sshgrad_y]=dsshdxi(gr.ssh,gr.DX,gr.DY);
    %% velocities
    gr.U=-corio.GOverF.*gr.sshgrad_y;
    gr.V= corio.GOverF.*gr.sshgrad_x;
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
    if ~isempty(RS)
        gr.L=gr.absUV./corio.f;
        kinVis=1e-6;
        gr.Re=gr.absUV.*gr.L/kinVis;
        gr.Ro=ones(size(gr.L));
        gr.Rrhines=earthRadius./gr.L;
        gr.Lrhines=sqrt(gr.absUV./corio.beta);
        gr.L_R=abs(RS.c./corio.f);
        gr.Bu=(gr.L_R./gr.L).^2;
    end
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
function [dsshdx,dsshdy]=dsshdxi(ssh,DX,DY)
    %% calc ssh gradients
    dsshdx=[diff(ssh,1,2), nan(size(ssh,1),1)] ./ DX;
    dsshdy=[diff(ssh,1,1); nan(1,size(ssh,2))] ./ DY;
end
function out=coriolisStuff(fields)
    %% omega
    out.Omega=angularFreqEarth;
    %% f
    out.f=2*out.Omega*sind(fields.lat);
    %% beta
    out.beta=2*out.Omega/earthRadius*cosd(fields.lat);
    %% gravity
    out.g=sw_g(fields.lat,zeros(size(fields.lat)));
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
