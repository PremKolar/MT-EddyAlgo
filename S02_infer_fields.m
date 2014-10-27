%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 04 - Apr - 2014 16:53:06
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates geostrophic data from ssh
function S02_infer_fields
    %% init
    DD = initialise('cuts',mfilename);
    %% read input file
    cut1 = load( DD.checks.passed(1).filenames);
    DD.coriolis = coriolisStuff(cut1.grids.lat);
    RS = getRossbyStuff(DD);
    %% spmd
    main(DD,RS)
    %% save info file
    conclude(DD)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main(DD,RS)
    %% infer mean ssh
    if ~exist([DD.path.root, 'meanSSH.mat'],'file')
        if DD.debugmode
            [JJ] = SetThreadVar(DD);
            spmd_meanSsh(DD,JJ);
        else
            spmd(DD.threads.num)
                [JJ] = SetThreadVar(DD);
                spmd_meanSsh(DD,JJ);
            end
        end
        MeanSsh = saveMean(DD);
    else
        load([DD.path.root, 'meanSSH.mat']);
    end
    %% calc fields
    if DD.debugmode
        [JJ] = SetThreadVar(DD);
        spmd_fields(DD,RS,JJ,MeanSsh);
    else
        spmd(DD.threads.num)
            [JJ] = SetThreadVar(DD);
            spmd_fields(DD,RS,JJ,MeanSsh);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MeanSsh = saveMean(DD)
    MeanSsh = nan(DD.map.window.sizePlus.y * DD.map.window.sizePlus.x,1);
    Meancount = 0;
    for ll = 1:DD.threads.num
        cur = load(sprintf('meanTmp%03d.mat',ll));
        MeanSsh = nansum([MeanSsh cur.Mean.SshSum],2);
        Meancount = Meancount + cur.Mean.count;
        system(sprintf('rm meanTmp%03d.mat',ll));
    end
    MeanSsh = reshape(MeanSsh,[DD.map.window.sizePlus.y, DD.map.window.sizePlus.x])/Meancount;
    save([DD.path.root, 'meanSSH.mat'],'MeanSsh')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RS = getRossbyStuff(DD)
    if DD.switchs.RossbyStuff
        RS.Lr = getfield(load([DD.path.Rossby.name DD.FieldKeys.Rossby{1} '.mat']),'data');
        RS.c = getfield(load([DD.path.Rossby.name DD.FieldKeys.Rossby{2}  '.mat']),'data');
    else
        RS = [];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spmd_meanSsh(DD,JJ)
    T = disp_progress('init','infering mean ssh');
    Mean.SshSum = nan(DD.map.window.sizePlus.y*DD.map.window.sizePlus.x,1);
    for jj = 1:numel(JJ)
        T = disp_progress('disp',T,numel(JJ),100);
        %% load
        ssh = extractdeepfield(load(JJ(jj).files),'grids.ssh')';
        %% mean ssh
        Mean.SshSum = nansum([Mean.SshSum, ssh],2);
    end
    Mean.count = numel(JJ);
    save(sprintf('meanTmp%03d.mat',labindex),'Mean');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spmd_fields(DD,RS,JJ,MeanSsh)
    T = disp_progress('init','infering fields');
    for jj = 1:numel(JJ)
        T = disp_progress('disp',T,numel(JJ),100);
        %% skip
        try
            alreadyFltrd = isfield(JJ(jj).files,'filtered');
        catch me
            disp(me.message)
            disp(['removing corrupt CUT - run relevant steps prior with DD.overwrite = false again!'])
            system(['rm ' JJ(jj).files])
        end
        
        if alreadyFltrd && ~DD.overwrite, dispM('skipping');continue; end
        cut = load(JJ(jj).files);
          %% filter
        if DD.switchs.filterSSHinTime
            %% not yet built
            if ~isfield(cut.grids,'sshRaw')
                cut.grids.sshRaw = cut.grids.ssh;
            end
            %% filter
            cut.grids.ssh = cut.grids.sshRaw - MeanSsh;
        end
        %%
        coriolis = coriolisStuff(cut.grids.lat);
        %% calc
        grids = geostrophy(cut.grids,coriolis,RS); %#ok<NASGU>
        %% write
        cut.filtered = true;
        save(JJ(jj).files,'grids','-append');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gr = geostrophy(gr,corio,RS)
    %% ssh gradient
    [gr.sshgrad_x,gr.sshgrad_y] = dsshdxi(gr.ssh,gr.dx,gr.dy);
    %% velocities
    gr.U = - corio.GOverF.*gr.sshgrad_y;
    gr.V = corio.GOverF.*gr.sshgrad_x;
    gr.absUV = hypot(abs(gr.U),abs(gr.V));
    %% 2d - deformation
    def = deformation(gr);
    gr.vorticity = def.dVdx - def.dUdy;
    gr.divergence = def.dUdx + def.dVdy;
    gr.stretch = def.dUdx - def.dVdy;
    gr.shear = def.dVdx + def.dUdy;
    %% okubo weiss
    %     gr.OW = .5*( - gr.vorticity.*2 + gr.divergence.*2 + gr.stretch.*2 + gr.shear.*2);
    %% or in 2d
    gr.OW = 2*(def.dVdx.*def.dUdy + def.dUdx.^2);
    %% assuming Ro = 1
    if ~isempty(RS)
        gr.L = gr.absUV./corio.f;
        kinVis = 1e-6;
        gr.Re = gr.absUV.*gr.L/kinVis;
        gr.Ro = ones(size(gr.L));
        gr.Rrhines = earthRadius./gr.L;
        gr.Lrhines = sqrt(gr.absUV./corio.beta);
        gr.L_R = abs(RS.c./corio.f);
        gr.Bu = (gr.L_R./gr.L).^2;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function def = deformation(grids)
    %% calc U gradients
    dUdy = diff(grids.U,1,1);
    dUdx = diff(grids.U,1,2);
    dVdy = diff(grids.V,1,1);
    dVdx = diff(grids.V,1,2);
    def.dUdy = dUdy([1:end, end], :)  ./ grids.dy;
    def.dUdx = dUdx(:,[1:end, end] )  ./ grids.dx;
    def.dVdy = dVdy([1:end, end], :)  ./ grids.dy;
    def.dVdx = dVdx(:,[1:end, end] )  ./ grids.dx;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dsshdx,dsshdy] = dsshdxi(ssh,dx,dy)
    %% calc ssh gradients
    dsshdx = diff(ssh,1,2);
    dsshdy = diff(ssh,1,1);
    dsshdx = dsshdx(:,[1:end, end])./ dx;
    dsshdy = dsshdy([1:end, end],:)./ dy;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = coriolisStuff(lat)
    %% omega
    out.Omega = angularFreqEarth;
    %% f
    out.f = 2*out.Omega*sind(lat);
    %% beta
    out.beta = 2*out.Omega/earthRadius*cosd(lat);
    %% gravity
    out.g = sw_g(lat,zeros(size(lat)));
    %% g/f
    out.GOverF = out.g./out.f;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
