function maxOWprocCalc
	load NC
	% 	load logOwMean
	NC.Yref=500;
	NC.yxref=[1570, 530];
	main(NC)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=funcs
	f.locNC=@(in,codi) getLocalPart(codistributed(nc_varget(in,'OkuboWeiss'),codi));
	f.locCo=@(in,codi) getLocalPart(codistributed(in,codi));
	f.ncPut=@(n,f,data)  nc_varput(n.(f).fileName ,n.(f).varName,data);
	f.ncPutBig=@(n,f,data,t,s)  nc_varput(n.(f).fileName ,n.(f).varName,data,[t,0,0],[1 s.Y s.X]);
	f.ncPutYref=@(n,f,data,t,s)  nc_varput(n.(f).fileName ,n.(f).varName,data,[t,0,0],[1 s.Z s.X]);
	f.ncPutXYref=@(n,f,data,t,s)  nc_varput(n.(f).fileName ,n.(f).varName,data,[t,0,0],[1 s.Z s.X]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [deepestLin,codi]=preploop(NC)
	codi=codistributor1d(3);
	f=funcs;
	%% get bathymetry
	[~,bathym]=max(~isnan(nc_varget(NC.files(1).full,'OkuboWeiss')));
	ndgridFromSize=@(in)  ndgrid(1:size(in,1),1:size(in,2));
	spmd
		deepest  =  f.locCo(bathym,codi)-1;
		deepest(deepest==0)=1;
		[Y,X]=ndgridFromSize(squeeze(deepest));
		deepestLin = sub2ind([NC.S.Z,NC.S.Y,NC.S.X], deepest(:), Y(:), X(:));
		%2nd deepest
		deepest  = deepest -1;
		deepest(deepest==0)=1;
		deepestLin = [reshape(deepestLin,1,[]) reshape(sub2ind([NC.S.Z,NC.S.Y,NC.S.X], deepest(:), Y(:), X(:)),1,[])];
		%3rd deepest
		deepest  = deepest -1;
		deepest(deepest==0)=1;
		deepestLin = [reshape(deepestLin,1,[]) reshape(sub2ind([NC.S.Z,NC.S.Y,NC.S.X], deepest(:), Y(:), X(:)),1,[])];
		%4th deepest
		deepest  = deepest -1;
		deepest(deepest==0)=1;
		deepestLin = [reshape(deepestLin,1,[]) reshape(sub2ind([NC.S.Z,NC.S.Y,NC.S.X], deepest(:), Y(:), X(:)),1,[])];
	end	%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calcMinZi(NC,tt)
	f=funcs;
	spmd
		mydata=log10OW(f.locNC(NC.currFile,NC.codi),nan);
		%% kill bottom layer
		mydata(NC.deepestLin)=nan;
		[owMin,MinZi]=nanmax(mydata(:,:,:),[], 1);
		MinZi=gcat(squeeze(MinZi),2,1);
		owMin=gcat(squeeze(owMin),2,1);
	end
	MinZi=MinZi{1};
	owMin=owMin{1};
	%% put to big files
	f.ncPutBig(NC.new,'minOWzi',MinZi,tt-1,NC.S);
	f.ncPutBig(NC.new,'minOW',owMin,tt-1,NC.S);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calcYref(NC,tt)
	f=funcs;
	spmd
		mydata= 		log10OW(f.locNC(NC.currFile,NC.codi),nan);
		owYref=gcat(squeeze(mydata(:,NC.Yref,:)),2,1);
	end
	owYref=owYref{1};
	%% put to big files
	f.ncPutYref(NC.new,'owYref',owYref,tt-1,NC.S);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calcXYref(NC,tt)
	f=funcs;
	owXY=log10OW(nc_varget(NC.currFile,'OkuboWeiss',[0 NC.yxref(1)-1 0],[inf 1 inf]),nan);
	f.ncPutXYref(NC.new,'owXYref',squeeze(owXY),tt-1,NC.S);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main(NC)
	%%
	[NC.deepestLin NC.codi]=preploop(NC);
	%%
	T=disp_progress('init','min OW''s')  ;
	for tt=1:NC.S.T
		T=disp_progress('show',T,NC.S.T);
		%% get min in z
		NC.currFile=NC.files(tt).full;
		%%
		 		calcMinZi(NC,tt)
		%%
		 		calcYref(NC,tt)
		%%
		calcXYref(NC,tt)
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OW]=log10OW(OW,dummy)
	tag=isnan(OW) | isinf(OW) | OW>=0;
	OW(tag)=dummy;
	OW=log10(-OW);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%