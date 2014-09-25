function [dYdXs]=diffCentered(diffLevel,X,Y)
	Y=reshape(Y,1,[]);
	X=reshape(X,1,[]);
	d=diff(Y,diffLevel);
	xa=X( ceil(diffLevel/2+1) : end-floor(diffLevel/2));
	xb=X(floor(diffLevel/2+1) : end-ceil(diffLevel/2 ));	
	Xout=mean([xa;xb],1);	
	dYdXs=spline(Xout,d,X);
end