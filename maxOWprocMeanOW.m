function maxOWprocMeanOW
	load NC
	logOwMean=main(NC);
	save logOwMean
	nc_varput(NC.new.OWmean.fileName ,NC.new.OWmean.varName,logOwMean);
	
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function logOwMean=main(NC)
	codi=codistributor1d(3);
	localNCget=@(in,codi) getLocalPart(codistributed(in,codi));
	spmd		
		logOwSum=localNCget(nc_varget(NC.files(1).full,'OkuboWeiss'),codi);
		nanflag=isnan(logOwSum);
		logOwSum(nanflag)=0;
		logOwCount=getLocalPart(ones(NC.S.Z,NC.S.Y,NC.S.X,codi));
		
		T=disp_progress('init','buildingmean')  ;
		for tt=2:NC.S.T
			T=disp_progress('show',T,NC.S.T);
			newOw=localNCget(nc_varget(NC.files(tt).full,'OkuboWeiss'),codi);
			nanflag=isnan(newOw);
			newOw(nanflag)=0;
			logOwSum=logOwSum+newOw;
			logOwCount(~nanflag)=logOwCount(~nanflag)+1;
		end
		logOwSum(logOwCount==0)=-30;
		logOwSum=gcat(logOwSum./logOwCount,3,1);
	end
	logOwMean=logOwSum{1};	
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%