%%%%%%%%%
% Created: 08-Apr-2013 19:50:46
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CM=doublemap(cb,cm1,cm2,centercol,center)
	%% get colorbardata
	if nargin<5
        center=0;
    end
    zlim=(get(cb,'ylim'))-center;
	ztick=(get(cb,'ytick'));
	zticklabel=(get(cb,'yticklabel'));
	%% resample to fit ticks	
	if numel(cm1)~=numel(cm2)
		if numel(cm2)>numel(cm1)
			cm1=resample(cm1,size(cm2,1),size(cm1,1));
		else
			cm2=resample(cm2,size(cm1,1),size(cm2,1));
		end
    end
	cm1r=resample(cm1,round(1000*abs(zlim(1))),round(1000*abs(zlim(1)) * abs(zlim(2)/zlim(1))));
	cm2r=resample(cm2,round(1000*zlim(2)),round(1000*zlim(2)));
	CM=[cm1r;flipud(cm2r)];
	%% blend in the middle
	midfilt=linspace(-1,zlim(2)/abs(zlim(1)),length(CM));
	%     gp=repmat(gauspuls(midfilt,1,1,-1/5)',1,3);
	gp=repmat(gauspuls(midfilt,1.5,1.5,-1/5)',1,3);
	centercolvec=repmat(centercol,size(CM,1),1);
	CM=(1-gp).*CM + gp.*centercolvec;
	%% correct for round errors
	CM(CM<0)=0;
	CM(CM>1)=1;
	%% reset to old params
	colormap(CM);
	caxis(zlim);
	set(cb,'ytick',ztick);
	set(cb,'yticklabel',zticklabel);
end