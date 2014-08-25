function maxOWprocMeanOW
    load NC
    main(NC);
    %     save OwMean
    %     nc_varput(NC.new.OWmean.fileName ,NC.new.OWmean.varName,OwMean);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main(NC)
    %     codi=codistributor1d(3);
    %     localNCget=@(in,codi) getLocalPart(codistributed(in,codi));
    %     spmd
    %         OwSum= localNCget(nc_varget(NC.files(1).full,'OkuboWeiss'),codi);
    %         nanflag=isnan(OwSum);
    %         OwSum(nanflag)=0;
    %         OwCount=getLocalPart(ones(NC.S.Z,NC.S.Y,NC.S.X,codi));
    %         T=disp_progress('init','buildingmean')  ;
    %         for tt=2:NC.S.T
    %             T=disp_progress('show',T,NC.S.T);
    %             newOw=localNCget(nc_varget(NC.files(tt).full,'OkuboWeiss'),codi);
    %             nanflag=isnan(newOw);
    %             newOw(nanflag)=0;
    %             OwSum=OwSum+newOw;
    %             OwCount(~nanflag)=OwCount(~nanflag)+1;
    %         end
    %         OwSum(OwCount==1)=nan;
    %         OwSum=gcat(OwSum./OwCount,3,1);
    %     end
    %     OwMean=OwSum{1};
    %
    %
    
    
    
    tic
    owatzz=nan(NC.S.T,2400,3600);
    files=NC.files;
    T=NC.S.T;
    parfor zz=1:NC.S.Z
        flipZtoT(zz,files,T,owatzz)
        labBarrier
    end
    toc
    
    
    
    for zz=1:NC.S.Z
        zfile=sprintf('OWat%02d.mat',zz)
        load(zfile)
        log10(MEAN./MEDIAN)
        
        
    end
    
    
    
    
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function flipZtoT(zz,files,T,owatzz)
%     if exist(sprintf('OWat%02d.mat',zz),'file'),return,end
    
    for tt=1:T
        if labindex==1,            fprintf('%2.2f%%\n',((zz/41-1)+tt/T)*100),end
        owatzz(tt,:,:)=  nc_varget(files(tt).full,'OkuboWeiss',[zz 0 0],[1 2400 3600]);
    end
    
    dispM('skew')
    out.SKEW=squeeze(skewness(owatzz,0));
    
    owatzz(owatzz>=0 | isnan(owatzz)  | isinf(owatzz) ) = nan;
    dispM('median')
    out.MEDIAN=squeeze(nanmedian(owatzz,1));
    dispM('mean')
    out.MEAN=squeeze(nanmean(owatzz,1));
    dispM('std')
    out.STD=squeeze(nanstd(owatzz,1));
  
    
    
    save(sprintf('OWat%02d.mat',zz),'-struct','out')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%