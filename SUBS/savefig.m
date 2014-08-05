%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 26-Apr-2014 19:05:25
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savefig(outdir,resOut,xdim,ydim,tit,debug,frmt,saveFig)
    set(0,'DefaultAxesFontSize', 12)
set(gcf,'renderer','painter')     
% set(gcf,'Visible','off')  
if nargin<6,debug=false;end
    if debug,disp('yo');	return;	end
    if nargin < 7,	frmt='dpdf';	end
    if nargin < 8,	saveFig=false;	end
    fname=[outdir,tit];
    mkdirp(outdir);
    if resOut==42
        quickhack(fname)
    else
        %% set up figure
        setupfigure(resOut,xdim,ydim)
        %% print
        printStuff(frmt,fname,resOut,xdim,ydim)
    end
    if saveFig
        save(gcf,[fname,'.mat'])
    end
    close(gcf);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printStuff(frmt,fname,resOut,xdim,ydim)
    if strcmp(frmt,'dpdf')
        eval(['print ',[fname,'.eps'] , ' -f -r',num2str(resOut),' -depsc ;'])
        system(['epstopdf ' fname '.eps']);
        system(['rm ' fname '.eps']);
    else
        fnfull=[fname,'.',frmt(2:end)];
        disp(['print ',fnfull , ' -f -r',num2str(resOut),' -',frmt,';'])
        eval(['print ',fnfull , ' -f -r',num2str(resOut),' -',frmt,';'])
        system(['convert -density ' num2str(resOut) 'x' num2str(resOut) ' -resize ' num2str(xdim) 'x' num2str(ydim) ' quality 100 ' fnfull ' ' fname '.pdf' ]);
%         system(['pdfcrop ' fname '.pdf ' fname '.pdf']);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setupfigure(resOut,xdim,ydim)

    resHere=get(0,'ScreenPixelsPerInch');   
%    ratioRes=resOut/resHere;
   ratioRes=1;
    set(gcf,'position',[0 0 [xdim ydim]/ratioRes]);    
    set(gcf,'paperunits','inch','papersize',[xdim ydim]/resOut,'paperposition',[0 0 [xdim ydim]/resOut]);
   
%     
% AxesH    = findobj(gcf, 'Type', 'text');
% YLabelHC = get(AxesH, 'YLabel');
% YLabelH  = [YLabelHC{:}];
% set(YLabelH, 'String', 'Y-label')
% TitleHC  = get(AxesH, 'Title');
% TitleH   = [TitleHC{:}];
% set(TitleH, 'String', 'The title');
%     set(gca,'FontSize',12)
  
   tt= [findall(gcf,'type','text') ;findall(gca,'type','text')]  ;
  for t=tt'    
   set(t,'FontSize',12)
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function quickhack(fname)
    disp(['quick print to ' [fname,'.png']])
    print(gcf, '-dpng', [fname,'.png'])
end
