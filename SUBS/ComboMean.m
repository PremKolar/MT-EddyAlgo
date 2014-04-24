function [meanAB]=ComboMean(Na,Nb,meanA,meanB)
	warning('off','MATLAB:divideByZero')
	meanA(isnan(meanA))=0;
	meanB(isnan(meanB))=0;
	meanAB=(Na.*meanA+Nb.*meanB)./(Na+Nb);
	warning('on','MATLAB:divideByZero')
end