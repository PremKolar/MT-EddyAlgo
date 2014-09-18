%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 26-Apr-2014 19:05:25
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%savegcf(gcf,outdir,resOut,xdim,ydim,tit,frmt,closeAfter)
function savefig(outdir,resOut,xdim,ydim,tit,frmt)
    set(gcf,'renderer','painter')
    set(gcf,'Visible','off')
    if nargin < 6,	frmt='dpdf';	end
    fname=[outdir,tit];
    mkdirp(outdir);
    %% set up gcfure
    [resHere,posOld]=setupfigure(resOut,xdim,ydim);
    %% print
    set(gcf,'Visible','off');
    printStuff(frmt,fname,resOut,xdim,ydim,resHere);
    set(gcf,'Visible','on');
    set(gcf,'position',posOld);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printStuff(frmt,fname,resOut,xdim,ydim,resHere)
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

