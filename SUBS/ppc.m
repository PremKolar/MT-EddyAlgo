function fig=ppc(in)
    [a,b,c]=size(in);
    if a==1 && b~=1 && c~=1
        in=squeeze(in);
    end
    %% draw
   fig=pcolor(flipud(full(double(squeeze(in)))));
   shading flat
   colorbar
%     fig=imagesc(flipud(full(double(squeeze(in)))));
%     colorbar
end