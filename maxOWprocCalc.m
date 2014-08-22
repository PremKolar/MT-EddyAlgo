function maxOWprocCalc
	load NC
	load logOwMean
	f=funcs;
	main(NC,logOwMean,f)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main(NC,OWmean,f)
	%% get bathymetry
	bath=getBathym(nc_varget(NC.files(1).full,'OkuboWeiss'));
	%%
	T=disp_progress('init','min OW''s')  ;
	for tt=1:NC.S.T
		T=disp_progress('show',T,NC.S.T);
		try daily=initDaily(NC,tt); catch exst; disp(exst); continue; end
		%% get min in z
		[owMin,MinZi]=spmdBlck(NC.files(tt).full,bath,f,OWmean);
		%% write daily
		f.ncPut(daily,'minOWzi',MinZi);
		f.ncPut(daily,'minOW',owMin);
		%% put to big files too
		f.ncPutBig(NC.new,'minOWzi',MinZi,tt-1,NC.S);
		f.ncPutBig(NC.new,'minOW',owMin,tt-1,NC.S);
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function daily=initDaily(NC,tt)
	daily.minOWzi.varName = NC.new.minOWzi.varName;
	daily.minOW.varName   = NC.new.minOW.varName;
	daily.minOWzi.fileName =  sprintf('%s%s_%04d.nc',NC.outdir,NC.new.minOWzi.fileName,tt);
	daily.minOW.fileName   =  sprintf('%s%s_%04d.nc',NC.outdir, NC.new.minOW.fileName ,tt);
	%%
	NC.iniNewNC(daily,'minOWzi',NC.new.dimNum(2:end),NC.new.dimName(2:end));
	NC.iniNewNC(daily,'minOW',  NC.new.dimNum(2:end),NC.new.dimName(2:end));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bath]=getBathym(OW)
	[Z,Y,X]=size(OW);
	OW2d=reshape(OW,[Z,Y*X]);
	[~,bathUpdown]=min(isnan(flipud(OW2d)),[],1);
	spmd
		bath=f.ncvOne(reshape( Z-bathUpdown + 1, [Y,X]),2);
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=funcs
	f.ncvOne = @(A,dim) getLocalPart(codistributed(A,codistributor1d(dim)));
	f.gCat = @(a,dim) gcat(squeeze(a),dim,1);
	f.ncPut=@(n,f,data)  nc_varput(n.(f).fileName ,n.(f).varName,data);
	f.ncPutBig=@(n,f,data,t,s)  nc_varput(n.(f).fileName ,n.(f).varName,data,[t,0,0],[1 s.Y s.X]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OW]=log10OW(OW,dummy)
	tag=isnan(OW) | isinf(OW) | OW>=0;
	OW(tag)=dummy;
	OW=log10(-OW);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [owMin,MinZi]=spmdBlck(currFile,mybath,f,OWmean)
	nanmaxFrom2toFloor = @(OW,bath) nanmax(OW(2:bath-1,:,:),[], 1);
	spmd
		mydata= 	f.ncvOne(log10OW(nc_varget(currFile,'OkuboWeiss'),nan),3);
		[owMin,MinZi]=nanmaxFrom2toFloor(mydata./OWmean,mybath);
		MinZi=f.gCat(MinZi-1,2); % correct for (2: ...)
		owMin=f.gCat(owMin,2);
	end
	MinZi=MinZi{1};
	owMin=owMin{1};
end
