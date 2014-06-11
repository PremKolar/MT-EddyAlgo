function fig=ppcolor(in)	
    %% if too large, downscale; keep ratio; make double
    [Y,X]=size(in);
    YoX=Y/X;
    Y(Y>200)=200;
    X=round(Y/YoX);
    %% draw
	fig=pcolor(downsize(double(in),X,Y));
	shading flat
	colorbar	
end