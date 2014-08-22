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
	nc_varput(NC.new.OWmean.fileName ,NC.new.OWmean.varName,logOwMean{1});
	save logOwMean
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function logOwMean=main(NC,f)
	T=disp_progress('init','calcing hor means of OW')  ;
	spmd
		logOwSum=f.ncvOne(nan(NC.S.Z,NC.S.Y,NC.S.X),3);
	end
	for tt=1:NC.S.T
		T=disp_progress('show',T,NC.S.T);
		%% get min in z
		spmdmDnansumlog=@(old,new)  multiDnansum(old, log10OW(new,nan));
		spmd
			newOw=f.ncvOne(nc_varget(NC.files(tt).full,'OkuboWeiss'),3);
			logOwSum=spmdmDnansumlog(logOwSum,newOw);
		end
	end
	spmd
		logOwMean=f.gCat(logOwSum/NC.S.T);
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=funcs
	f.ncvOne = @(A,dim) getLocalPart(codistributed(A,codistributor1d(dim)));
	f.gCat = @(a,dim) gcat(squeeze(a),dim,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OW]=log10OW(OW,dummy)
	tag=isnan(OW) | isinf(OW) | OW>=0;
	OW(tag)=dummy;
	OW=log10(-OW);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

