%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 04-Apr-2014 16:53:06
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculates geostrophic data from ssh
function S02_infer_fields
    %% init
    DD=initialise('cuts',mfilename);
    %% read input file
    cut1=load( DD.checks.passed(1).filenames);
    DD.coriolis=coriolisStuff(cut1.grids.lat);
    RS=getRossbyStuff(DD,cut1.grids);
    %% spmd
    main(DD,RS)
    %% save info file
    conclude(DD)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main(DD,RS)
    if DD.debugmode
        spmd_body(DD,RS);
    else
        spmd(DD.threads.num)
            spmd_body(DD,RS);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RS=getRossbyStuff(DD,gr)
    if DD.switchs.RossbyStuff
        RS.Lr=getfield(load([DD.path.Rossby.name 'RossbyRadius.mat']),'out');
        RS.c=getfield(load([DD.path.Rossby.name 'RossbyPhaseSpeed.mat']),'out');
        RS.Lr(RS.Lr<0 | RS.Lr > 100*nanmedian(abs(RS.Lr(:))))=nan;
        RS.c(RS.c<0 | RS.c > 100*nanmedian(abs(RS.c(:))))=nan;
        RS.LrInc.y=smooth2a(RS.Lr./gr.DY,10);
        RS.LrInc.x=smooth2a(RS.Lr./gr.DX,10);
        RS.LrInc.x=double(NeighbourValue(isnan(RS.LrInc.x),RS.LrInc.x));
        RS.LrInc.y=double(NeighbourValue(isnan(RS.LrInc.y),RS.LrInc.y));
    else
        RS=[];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        grids=geostrophy(cut.grids,coriolis,RS);
        if ~isfield(grids,'sshRaw')
            grids.sshRaw=grids.ssh;
            grids.ssh=filterStuff(cut.grids,RS);
        end
        %% write
        save(JJ(jj).files,'grids','-append');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sshHighPass=filterStuff(gr,RS)   
    %% get center, minor and major axis for ellipse      
    RossbyEqFlag=abs(gr.lat)<5  ;													%TODO {put before loop...
    semi.x=10*ceil(max(nanmedian(RS.LrInc.x(~RossbyEqFlag),2)));			%...
	     semi.y=10*ceil(max(nanmedian(RS.LrInc.y(~RossbyEqFlag),2)));		%.}
    [sshHighPass]=ellipseFltr(semi,gr.ssh);										
    %
    %     JET=repmat(jet,3,1);
    %     figure(1)
    %     b=gr.ssh - min(gr.ssh(:));
    %     %      surf(double(gr.lon),double(gr.lat), b*5)
    %     % X= gr.lon(1:2:end,1:2:end);
    %     % Y= gr.lat(1:2:end,1:2:end);
    %     % Z= b(1:2:end,1:2:end)*5
    %     %      surface(double(X),double(Y), Z,'EdgeColor',[.8 .8 .8],'FaceColor','none')
    %     % hold on
    %     contour3(double(gr.lon),double(gr.lat), b*5,(min(b(:)):0.01:max(b(:)))*5)
    %     axis tight equal;    view(3)
    %     % cb1=colorbar
    %     colormap(JET(round(.25/3*size(JET,1)):end,:))
    %     zt=get(gca,'ztick')
    %     set(gca,'zticklabel',num2str(zt'/5*100))
    %     % set(cb1,'ytick',[])
    %     xlabel('lon')
    %     ylabel('lat')
    %     zlabel('[cm]')
    %     tit='Non-filtered SSH'
    %     title(tit)
    %     savefig('./',200,8*200,6*200, space2underscore(tit))
    %     figure(2)
    %     a=sshHighPass - min(sshHighPass(:));
    %     %  Z= a(1:2:end,1:2:end)*5
    %     %      surface(double(X),double(Y), Z,'EdgeColor',[.8 .8 .8],'FaceColor','none')
    %     % hold on
    %     %     surf(double(gr.lon),double(gr.lat), a*5,'facecolor','black','facealpha',.5)
    %     contour3(double(gr.lon),double(gr.lat),a*5,(min(a(:)):0.01:max(a(:)))*5)
    %     axis tight equal;    view(3)
    %     % cb2=colorbar
    %     colormap(jet)
    %     zt=get(gca,'ztick')
    %     set(gca,'zticklabel',num2str(zt'/5*100))
    %     xlabel('lon')
    %     ylabel('lat')
    %     zlabel('[cm]')
    %     % set(cb2,'ytick',[])
    %     tit='High-pass filtered SSH'
    %     title(tit)
    %     savefig('./',200,8*200,6*200, space2underscore(tit))
    %     figure(3)
    %     c=b-a;
    %     %   Z= c(1:2:end,1:2:end)*5
    %     %      surface(double(X),double(Y), Z,'EdgeColor',[.8 .8 .8],'FaceColor','none')
    %     %     surf(double(gr.lon),double(gr.lat), c*5,'facecolor','black','facealpha',.5)
    %     contour3(double(gr.lon),double(gr.lat),c*5,(min(c(:)):0.01:max(c(:)))*5)
    %     axis tight equal;    view(3)
    %     colormap(JET(round(.25/3*size(JET,1)):end,:))
    %     zt=get(gca,'ztick')
    %     set(gca,'zticklabel',num2str(zt'/5*100))
    %     xlabel('lon')
    %     ylabel('lat')
    %     zlabel('[cm]')
    %     tit='Subtracted Low-Pass component'
    %     title(tit)
    %     savefig('./',200,8*200,6*200, space2underscore(tit))
    %
    %
    %%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function def=deformation(grids)
	%% calc U gradients	
	dUdy=diff(grids.U,1,1);
	dUdx=diff(grids.U,1,2);
	dVdy=diff(grids.V,1,1);
	dVdx=diff(grids.V,1,2);
	def.dUdy= dUdy([1:end, end], :)  ./ grids.DY;
	def.dUdx= dUdx(:,[1:end, end] )  ./ grids.DX;
	def.dVdy= dVdy([1:end, end], :)  ./ grids.DY;
	def.dVdx= dVdx(:,[1:end, end] )  ./ grids.DX;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dsshdx,dsshdy]=dsshdxi(ssh,DX,DY)
    %% calc ssh gradients
	 dsshdx=diff(ssh,1,2);
    dsshdy=diff(ssh,1,1);	 
	  dsshdx=dsshdx(:,[1:end, end])./ DX;
    dsshdy=dsshdy([1:end, end],:)./ DY;
%   
% 	 
% 	 dsshdx=[diff(ssh,1,2), nan(size(ssh,1),1)] ./ DX;
%     dsshdy=[diff(ssh,1,1); nan(1,size(ssh,2))] ./ DY;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=coriolisStuff(lat)
    %% omega
    out.Omega=angularFreqEarth;
    %% f
    out.f=2*out.Omega*sind(lat);
    %% beta
    out.beta=2*out.Omega/earthRadius*cosd(lat);
    %% gravity
    out.g=sw_g(lat,zeros(size(lat)));
    %% g/f
    out.GOverF=out.g./out.f;
end

