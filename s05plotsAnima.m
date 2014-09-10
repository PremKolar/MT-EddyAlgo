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
	%% get ssh data
	try
		cut=load(EE.filename.cut);
		%% get contours
		cont=load(EE.filename.cont);
	catch %#ok<CTCH>
		% TODO update dt!!!
		skip=true;
		%         return
	end
	%% put all eddies into a struct: ee(number of eddies).characteristica
	ee=eddies2struct(cont.all,DD.thresh.corners);
	%% remeber date
	[ee(:).daynum]=deal(JJ.daynums);
	%% avoid out of bounds integer coordinates close to boundaries
	[ee_clean,~]=CleanEDDies(ee,cut,DD.contour.step);
	%% find them
	EE=find_eddies(EE,ee_clean);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EE=find_eddies(EE,ee_clean)
	%
	makejpegs(EE,ee_clean);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function makejpegs(EE,ee)	
	load(EE.filename.self);
	load(EE.filename.cont);
	figure(labindex);clf;
	for kk=1:numel(ee)
		x=ee(kk).coordinates.exact.x;
		y=ee(kk).coordinates.exact.y;
		hold on;
		plot(x,y,'color',rainbow(1,1,1,kk,numel(ee)));
	end
	for kk=1:numel(anticyclones)
		x=anticyclones(kk).coordinates.exact.x;
		y=anticyclones(kk).coordinates.exact.y;
		hold on
		plot(x,y,'color','black','linewidth',2)
	end
	for kk=1:numel(cyclones)
		x=cyclones(kk).coordinates.exact.x;
		y=cyclones(kk).coordinates.exact.y;
		hold on
		plot(x,y,'--','color','black','linewidth',2)
	end
	axis tight;
	title(datestr(ee(1).daynum));
 	savefig2png4mov('./',100,800,600,datestr(ee(1).daynum,'yymmdd'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [EE]=eddies2struct(CC,thresh)
	EE=struct;
	ii=1;cc=0;
	while ii<size(CC,1);
		len=  CC(ii,2);% contourc saves the length of each contour before appending the next
		if len>=thresh
			cc=cc+1;
			EE(cc).level=CC(ii,1);
			EE(cc).circum.length= len;
			EE(cc).coordinates.exact.x=CC(1+ii:ii+EE(cc).circum.length,1);
			EE(cc).coordinates.exact.y=CC(1+ii:ii+EE(cc).circum.length,2);
			EE(cc).coordinates.int.x=int32(round(EE(cc).coordinates.exact.x));
			EE(cc).coordinates.int.y=int32(round(EE(cc).coordinates.exact.y));
		end
		ii=ii+len+1; % jump to next eddy for next iteration
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ee,cut]=CleanEDDies(ee,cut,contstep)  %#ok<INUSD>
	[cut.dim.Y,cut.dim.X]=size(cut.grids.ssh);
	for jj=1:numel(ee)
		x=ee(jj).coordinates.int.x;
		y=ee(jj).coordinates.int.y;
		%% the following also takes care of the overlap from S00 in the global case
		% x(x>cut.window.size.X)= x(x>cut.window.size.X)-cut.window.size.X ;
		
		x(x>cut.dim.X)=cut.dim.X;
		y(y>cut.dim.Y)=cut.dim.Y;
		x(x<1)=1;
		y(y<1)=1;
		ee(jj).coordinates.int.x=x;
		ee(jj).coordinates.int.y=y;
	end
end
