% written by NK,1406
% 	inspired by Greg Reeves, March 2009.
% 	Division of Biology
% 	Caltech
% 	Inspired by "smooth2", written by Kelly Hilands, October 2004
% 	Applied Research Laboratory
% 	Penn State University
% 	Developed from code written by Olof Liungman, 1997
% 	Dept. of Oceanography, Earth Sciences Centre
% 	G�teborg University, Sweden
% 	E-mail: olof.liungman@oce.gu.se
function Mout = smooth2FuncOfY(Min,Nsouth,Nnorth)
    [Y,X] = size(Min);
    padwidth=max([Nsouth,Nnorth])    ;
    Mpad=padarray(Min,[padwidth padwidth],'symmetric');
    [Ypad,Xpad] = size(Mpad);
    %%
    Nm=linspace(Nsouth,Nnorth,Y);
    Nmerid=padarray((Nm),[0 padwidth],'replicate');
    %%
    
    
    
    sigmaMax=max([Nsouth,Nnorth]);    
    filterWidth=2*sigmaMax;
  
    
    
    for ii=1:Ypad
        sigma=Nmerid(ii);
        alpha(ii)=(2*filterWidth+1)./(2*sigma);
        fltrY(ii,:)=gausswin(2*filterWidth+1,alpha(ii))';        
    end
    ey = spdiags(fltrY,(-max(Nmerid):max(Nmerid)),Ypad,Ypad);
    %%
       A = isnan(Mpad);
    Mpad(A) = 0;
    Mout = ey*Mpad ;
    nrmlize = ey*(~A);
     %%
         
     for ii=1:Ypad
        fltrX = repmat(fltrY(ii,:),Xpad,1);       
        ex = spdiags(fltrX,-max(Nmerid):max(Nmerid),Xpad,Xpad);
        Mout(ii,:)=Mout(ii,:)*ex;
          nrmlize(ii,:)=nrmlize(ii,:)*ex;
    end
    
    
  
  
    nrmlize(A) = NaN;  
    Mout = Mout./nrmlize;
    Mout = Mout(padwidth+1:end-padwidth,padwidth+1:end-padwidth);
    
end











