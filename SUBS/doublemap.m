%%%%%%%%%
% Created: 08-Apr-2013 19:50:46
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CM=doublemap(abc,cm1,cm2,centercol,alpha)
	celldeal= @(x) x{:};
	tenfac= @(a,b) celldeal(num2cell([a,b] * 10^(max([ length(num2str(a)) length(num2str(b)) ])-2)));
	
	%% get colorbardata
	if nargin<=4
		alpha=1;
	end
	%% resample to fit ticks
	
	da=diff(abc([1 2]));
	db=diff(abc([2 3]));
	if any(abs(mod([da db],1)))
		[da,db]	= tenfac(da,db);
	end
	
	
	cm1=resample(cm1,da,db);
	cm1=resample(cm1,100,size(cm1,1));
	cm2=resample(cm2,100,size(cm2,1));
	CM=[cm1;flipud(cm2)];
	%% blend in the middle
	nrm=@(x) (x-min(x(:)))/max(x(:)-min(x(:)));
	gp=repmat(gausswin(size(CM,1),alpha),1,3);
	centercolvec=repmat(centercol,size(CM,1),1);
	CM=nrm((1-gp).*CM + gp.*centercolvec);
	%% correct for round errors
	CM(CM<0)=0;
	CM(CM>1)=1;
	%% reset to old params
	colormap(CM);
end
