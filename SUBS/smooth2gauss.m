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
function matrixOut = smooth2gauss(matrixIn,Nr,Nc)
	N(1) = Nr;
	if nargin < 3, N(2) = N(1); else N(2) = Nc; end
	[row,col] = size(matrixIn);
	filL=repmat((cos(linspace(-pi/2,pi/2,2*N(1)+1))+1)/2,row,1);
	filR=repmat((cos(linspace(-pi/2,pi/2,2*N(2)+1))+1)/2,col,1);

    eL = spdiags(filL,(-N(1):N(1)),row,row);
	eR = spdiags(filR,(-N(2):N(2)),col,col);
	A = isnan(matrixIn);
	matrixIn(A) = 0;
	nrmlize = eL*(~A)*eR;
	nrmlize(A) = NaN;
	matrixOut = eL*matrixIn*eR;
	matrixOut = matrixOut./nrmlize;
end











