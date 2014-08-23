function maxOWprocCalc
	load NC
	% 	load logOwMean
	main(NC)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=funcs
	f.locNC=@(in,codi) getLocalPart(codistributed(nc_varget(in,'OkuboWeiss'),codi));
	f.locCo=@(in,codi) getLocalPart(codistributed(in,codi));
	f.ncPut=@(n,f,data)  nc_varput(n.(f).fileName ,n.(f).varName,data);
	f.ncPutBig=@(n,f,data,t,s)  nc_varput(n.(f).fileName ,n.(f).varName,data,[t,0,0],[1 s.Y s.X]);
	f.ncPutYref=@(n,f,data,t,s)  nc_varput(n.(f).fileName ,n.(f).varName,data,[t,0,0],[1 s.Z s.X]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main(NC)
	codi=codistributor1d(3);
	f=funcs;
	Yref=500;
	%% get bathymetry
	[~,bathym]=max(~isnan(nc_varget(NC.files(1).full,'OkuboWeiss')));
	ndgridFromSize=@(in)  ndgrid(1:size(in,1),1:size(in,2));
	spmd
		deepest  =  f.locCo(bathym,codi)-1;
		deepest(deepest==0)=1;
		[Y,X]=ndgridFromSize(squeeze(deepest));
		deepestLin = sub2ind([NC.S.Z,NC.S.Y,NC.S.X], deepest(:), Y(:), X(:));
	end
	%%
	T=disp_progress('init','min OW''s')  ;
	for tt=1:NC.S.T
		T=disp_progress('show',T,NC.S.T);
		%% get min in z
		currFile=NC.files(tt).full;
		spmd
			mydata= 		log10OW(f.locNC(currFile,codi),nan);
			%% kill bottom layer
			mydata(deepestLin)=nan;
			[owMin,MinZi]=nanmax(mydata(:,:,:),[], 1);
			MinZi=gcat(squeeze(MinZi),2,1); 
			owMin=gcat(squeeze(owMin),2,1);
			owYref=gcat(squeeze(mydata(:,Yref,:)),2,1);
		end
		MinZi=MinZi{1};
		owMin=owMin{1};
		owYref=owYref{1};
		%% put to big files
		f.ncPutBig(NC.new,'minOWzi',MinZi,tt-1,NC.S);
		f.ncPutBig(NC.new,'minOW',owMin,tt-1,NC.S);
		f.ncPutYref(NC.new,'owYref',owYref,tt-1,NC.S);
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OW]=log10OW(OW,dummy)
	tag=isnan(OW) | isinf(OW) | OW>=0;
	OW(tag)=dummy;
	OW=log10(-OW);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%