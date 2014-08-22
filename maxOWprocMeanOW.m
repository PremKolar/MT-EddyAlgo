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
	spmd
		dispM('getting logOwSum')
		tic
		codi=codistributor1d(3);
		logOwSum=nan(NC.S.Z,NC.S.Y,NC.S.X,codi);
		toc
		dispM('getting logOwCount')
		tic
		logOwCount=ones(NC.S.Z,NC.S.Y,NC.S.X,codi);
		toc
	end
	for tt=1:NC.S.T
		tic
		fprintf('prog %02d%%\n',round(tt/NC.S.T*100))
		newOw=spmdCoDisp(nc_varget(NC.files(tt).full,'OkuboWeiss'),3);
		toc
		tic
		nanflag=isnan(newOw);
		newOw(nanflag)=0;
		logOwSum=logOwSum+newOw;
		logOwCount(~nanflag)=logOwCount(~nanflag)+1;
		toc
	end
	logOwSum(logOwCount==0)=nan;
	logOwSum=logOwSum./logOwCount;
	logOwMean=gather(logOwSum);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [out,codisp,lim]=spmdCoDisp(in,dim)
function [out]=spmdCoDisp(in,dim)
	
	spmd
		out=codistributed(in,codistributor1d(dim));
		% 		codisp=getCodistributor(out);
		% 		parti=getfield(codisp,'Partition');
		% 		LIM.a=[1 cumsum(parti(1:end-1))];
		% 		LIM.b=LIM.a + parti;
		% 		lim.a=LIM.a(labindex);
		% 		lim.b=LIM.b(labindex);
		% 		lim.len=lim.b-lim.a+1;
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%