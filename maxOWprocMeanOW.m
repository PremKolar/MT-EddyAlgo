%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 17-Jul-2014 23:52:44
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NKkk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maxOWprocMeanOW
	load NC
	logOwMean=main(NC);
	nc_varput(NC.new.OWmean.fileName ,NC.new.OWmean.varName,logOwMean);
	save logOwMean
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function logOwMean=main(NC)
	codi=codistributor1d(3);
	spmd
		logOwSum=getLocalPart(nan(NC.S.Z,NC.S.Y,NC.S.X,codi));
		logOwCount=getLocalPart(ones(NC.S.Z,NC.S.Y,NC.S.X,codi));
	end
	for tt=1:NC.S.T
		tic
		fprintf('prog %02d%%\n',round(tt/NC.S.T*100))
		newOw=spmdCoDisp(nc_varget(NC.files(tt).full,'OkuboWeiss'),codi);
		toc
		tic
		spmd
			nanflag=isnan(newOw);
			newOw(nanflag)=0;
			logOwSum=logOwSum+newOw;
			logOwCount(~nanflag)=logOwCount(~nanflag)+1;
		end
		toc
	end
	spmd
		logOwSum(logOwCount==0)=nan;
		logOwSum=gcat(logOwSum./logOwCount,3,1);
	end
	logOwMean=logOwSum{1};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [out,codisp,lim]=spmdCoDisp(in,dim)
function [out]=spmdCoDisp(in,codi)
	spmd
		out=getLocalPart(codistributed(in,codi));
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%