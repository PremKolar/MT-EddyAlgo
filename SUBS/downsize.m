function [OUT]=downsize(IN,xout,yout)

	[ydim,xdim]=size(IN);
xinc=max([round(xdim/xout),1]);
yinc=max([round(ydim/yout),1]);
OUT=IN(1:yinc:end,1:xinc:end);
end

