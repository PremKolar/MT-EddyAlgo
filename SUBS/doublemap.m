%%%%%%%%%
% Created: 08-Apr-2013 19:50:46
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CM=doublemap(abc,cm1,cm2,centercol,alpha)
    %% get colorbardata
    if nargin<4
        alpha=1;
    end
    %% resample to fit ticks
%     dt=diff(abc([1 3]));
    dat=diff(abc([1 2]));
    dbt=diff(abc([2 3]));
    
    da=round(200*dat/dbt);
    db=200;
    
    cm1=resample(cm1,da,size(cm1,1));
%     cm1=resample(cm1,100,size(cm1,1));
    cm2=resample(cm2,200,size(cm2,1));
    CM=[cm1;flipud(cm2)];
    %% blend in the middle
    nrm=@(x) (x-min(x(:)))/max(x(:)-min(x(:)));
   
    gw1=gausswin(size(cm1,1)*2,alpha)
     gw2=gausswin(size(cm2,1)*2,alpha)
    gw=[gw1(1:length(gw1)/2); gw2(length(gw2)/2+1:end) ]
     
   gp=repmat(gw,1,3);
    centercolvec=repmat(centercol,size(CM,1),1);
    CM=nrm((1-gp).*CM + gp.*centercolvec);
    %% correct for round errors
    CM(CM<0)=0;
    CM(CM>1)=1;
    %% reset to old params
    colormap(CM);  
end
