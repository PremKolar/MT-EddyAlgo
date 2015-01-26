%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 26-Apr-2014 19:05:25
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savefig(outdir,resOut,xdim,ydim,tit,frmt,info)
    set(gcf,'renderer','painter')
    set(gcf,'Visible','off')
    if nargin < 6,	frmt='dpdf';	end
    if nargin < 7,	info=[]    ;	end
    fname=[outdir,tit];
    mkdirp(outdir);
    %% set up gcfure
    [resHere,posOld]=setupfigure(resOut,xdim,ydim);
    %% print
    fnamepdf=printStuff(frmt,fname,resOut,xdim,ydim,resHere);
    if nargin == 7,
        appendPdfMetaInfo(info,fnamepdf);
    end
    set(gcf,'Visible','on');
    set(gcf,'position',posOld);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function appendPdfMetaInfo(info,fnamepdf) %#ok<INUSL>
    structfn=sprintf('%03d_pdfinfo.mat',labindex);
    save(structfn,'info')    
    system(sprintf('pdftk %s attach_files %s output %s.tmp.pdf',fnamepdf,structfn,fnamepdf));
    system(sprintf('mv %s.tmp.pdf %s',fnamepdf,fnamepdf));    
    system(['rm ' structfn]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fnamepdf=printStuff(frmt,fname,resOut,xdim,ydim,resHere)
    if strcmp(frmt,'dpdf')
        eval(['print ',[fname,'.eps'] , ' -f -r',num2str(resOut),' -depsc ;'])
        system(['epstopdf ' fname '.eps']);
        system(['rm ' fname '.eps']);
    else
        if strcmp(fname(end-length(frmt)+1),frmt )
            fnfull=fname;
        else
            fnfull=[fname,'.',frmt(2:end)];
        end
        todo=['print ',fnfull , ' -f -r',num2str(resOut),' -',frmt,';'];
        disp(todo)
        eval(todo)
        system(['convert -density ' num2str(resHere) 'x' num2str(resHere) ' -resize ' num2str(xdim) 'x' num2str(ydim) ' -quality 100 ' fnfull ' ' fname '.pdf' ]);
    end   
    fnamepdf=[fname '.pdf'];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [resHere,posNow]=setupfigure(resOut,xdim,ydim)
    resHere=get(0,'ScreenPixelsPerInch');
    %    ratioRes=resOut/resHere;
    ratioRes=1;
    try
        posNow=get(gcf,'position');
        set(gcf,'position',[0 0 [xdim ydim]/ratioRes]);
    catch
        disp('not setting up position for docked fig')
    end
    set(gcf,'paperunits','inch','papersize',[xdim ydim]/resOut,'paperposition',[0 0 [xdim ydim]/resOut]);
    set(gca,'FontSize',12)
    set(findall(gcf,'type','text'),'FontSize',12)
end

