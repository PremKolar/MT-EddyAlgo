%pcolor : fig=ppc(in,1) or  fig=ppc(in)
%imagesc: fig=ppc(in,2)
function fig=ppc(in,type,noshrink)
    if nargin==1,type=1; end
    if nargin<3,noshrink=0; end
    in=shapeUp(in,noshrink);
    fig=drawIt(in,type);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function in=shapeUp(in,noshrink)
    in=full(double(squeeze(in)));
    if ndims(in)>2 %#ok<ISMAT>
        in=reshape(in,size(in,1),[]);
    end
    if ~noshrink
        in=shrunky(in);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function in=shrunky(in)
    fac=floor(size(in)./[1000 1000]);
    fac(fac<1)=1;
    in=in(1:fac(1):end,1:fac(2):end);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fig=drawIt(in,type)
    switch type
        case 1
            fig=pcolor(in);
            shading flat;
        case 2
            fig=imagesc(flipud(in));
    end
    colorbar;
    axis tight;
    if numel(unique(in(:)))>2
        nstd=nanstd(in(:));
        me=nanmean(in(:));
        caxis([me-2*nstd me+2*nstd])
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
