%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 17-Jul-2014 23:52:44
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NKkk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maxOWprocMeanOW
	load NC
	f=funcs;
	logOwMean=main(NC,f);
	nc_varput(NC.new.OWmean.fileName ,NC.new.OWmean.varName,logOwMean);
	save logOwMean
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function logOwMean=main(NC,f)
	spmdmDnansumlog=@(old,new)  multiDnansum(old, log10OW(new,nan));
	spmd	
		logOwSum=f.cod(nan(NC.S.Z,NC.S.Y,NC.S.X),3);
		labBarrier;
		for tt=drange(1:NC.S.T)
			fprintf('lab%02d at tt=%02d\n',labindex,tt)	
			newOw=f.cod(nc_varget(NC.files(tt).full,'OkuboWeiss'),3);
			logOwSum=spmdmDnansumlog(logOwSum,newOw);
		end
	end
	logOwMean=gather(logOwSum/NC.S.T);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=funcs
	f.cod = @(A,dim) codistributed(A,codistributor1d(dim));
	% 	f.cod = @(A) getLocalPart(codistributed(A,codistributor1d(ndims(A))));
	f.gCat = @(a,dim) gcat(squeeze(a),dim,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OW]=log10OW(OW,dummy)
	tag=isnan(OW) | isinf(OW) | OW>=0;
	OW(tag)=dummy;
	OW=log10(-OW);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

