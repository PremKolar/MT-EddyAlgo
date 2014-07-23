% http://en.wikipedia.org/wiki/Standard_deviation#Combining_standard_deviat
% ions
function [stdAB]=ComboStd(Na,Nb,meanA,meanB)
	Na(isnan(Na))=0;
	Nb(isnan(Nb))=0;
	meanA(isnan(meanA))=0;
	meanB(isnan(meanB))=0;
	%%
	Nall=Na+Nb;
	warning('off','MATLAB:divideByZero')
	stdAB=sqrt((Na.*meanA.^2 + Nb.*meanB.^2)./Nall + (Na.*Nb.*(meanA-meanB).^2)./(Nall.^2));
	warning('on','MATLAB:divideByZero')
end