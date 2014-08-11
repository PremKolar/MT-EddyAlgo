

function ampOra
    addpath(genpath('./'))
    dbstop if error
    DD=initialise('',mfilename);    
    here=pwd;
    D={'iq2'; 'iq4'; 'iq6'; 'iq8'; 'ch400amparea'; 'ch400';'iq5-1d'};
    basedir=['/scratch/uni/ifmto/u300065/'];
    locCo = @(x) (x{labindex});
    spmd(numel(D))
        d=locCo(D);
        spmdloop(d,basedir,DD,here);
        labBarrier
    end    
    system(['pdfjam -o tmp.pdf crpd*pdf'])
    system(['pdfcrop  --margins "1 3 1 1" tmp.pdf jmd.pdf'])
    system(['okular jmd.pdf'])
end

function spmdloop(dd,bd,DD,here)
    outfile= ['tracks_' dd '.mat'];
    if ~exist(outfile,'file')
        cd([bd 'data' dd '/TRACKS/']) ;
        [~,file]=system(['ls -hlSr | tail -n 1 | tail -n 1 | grep -oE ''[^ ]+$''']);
        system(['cp ' file(1:end-1) ' ' bd 'iq2/' outfile]);
        cd(here) ;
    end
    savestuff(outfile,DD);
end

function savestuff(in,DD)
    pl= @(ss,DD) AOplots(getfield(load(ss),'trck'),DD);
    pl(in,DD)
    set(gcf,'position',[5 40 1912 437])    ;
    title(in)
    savefig('./',80,1600,500,sprintf('tmp_%d',labindex))
    system(sprintf('pdfcrop --margins "1 3 1 1" tmp_%d.pdf crpd_%d.pdf ',labindex,labindex))
end
function AOplots(track,DD)
    trck=track(2:end-1);
    try
        ar=extractfield(cell2mat(extractfield(trck,'area')),'intrp');
    catch
        return
    end
    ra=sqrt(ar/pi); %#ok<*NASGU>
    amp=extractfield(cell2mat(extractfield(cell2mat(extractfield(trck,'peak')),'amp')),'to_contour');
    age=cat(2,trck.age);
    iq=cat(2,trck.isoper);
    vol=    extractdeepfield(trck,'volume.total') ;
    nrm=@(x) (x-min(x))/max(x-min(x));
    figure(1)
    clf
    subplot(2,1,1,'align')
    hold on
    qu=((vol./ar.^(3/2)));
    quo=([1 qu(2:end)./qu(1:end-1)]);
    quo=log(quo);
    cm1=(hot);
    cm2=(hot);
    col=doublemap([-1,0,1],cm1,cm2,[.3 .3 1],3);
    col(end-2:end,:)=repmat([1 .5 0],3,1);
    col(1:3,:)=repmat([0 1 0],3,1);
    colormap(col);
    xtck=nan(size(trck));
    xtckl=num2cell(nan(size(trck)));
    set(gca,'ytick',[]);
    cb=colorbar('location','West');
    mq=max(abs(quo));
    caxis([-mq mq]);
    mqs=floor(10*mq)/10;
    ct=-mqs:2*mqs/4:mqs;
    set(cb,'ytick',ct)
    cbtckl=cellfun(@(cc) sprintf('% .1f', cc) ,num2cell((exp(ct))),'uniformoutput',false);
    cbtckl{3}=1;
    set(cb,'yticklabel',cbtckl)
    cblv=linspace(-mq, mq, size(col,1));
    [QUO,CBLV]=meshgrid(quo,cblv);
    [~,cblvPos]=min(abs(QUO-CBLV),[],1);
    for ii=1:1:numel(trck)
        y=extractdeepfield(trck(ii),'coordinates.exact.y');
        x=extractdeepfield(trck(ii),'coordinates.exact.x')-ii/numel(trck)*1000;
        plot(x,y,'color',col(cblvPos(ii),:))
        drawnow
        axis  tight
        xtck(ii) = mean(x);
        xtckl(ii) = {num2str(age(ii))};
    end
    set(gca,'xtick',[])
    a=axis;
    xAX=linspace(a(1),a(2),numel(age));
    [~,xr]=sort(xtck);
    if numel(xr)>20
        ii=round(linspace(1,numel(xr),20));
        xr=xr(ii);
    end
    set(gca,'xtick',xtck(xr))
    set(gca,'xticklabel',repmat([],1,numel(xr)))
    
    
    
    
    
    
    
    
    
    
    %%
    subplot(2,1,2,'align')
    try
        thresh=DD.thresh.amArea;
    catch
        thresh=DD.thresh.ampArea;
    end
    THR=log(repmat(thresh,length(xAX),1));
    xAXdouble=repmat(xAX',1,2);
    hold on
    
    
    
    IQ=log(iq)-mean(log(iq));
    difabs= @(a) fliplr(log([1 a(2:end)./a(1:end-1)]))';
    %      difabs= @(a) fliplr(([1 a(2:end)./a(1:end-1)]));
    A= [difabs(vol.^(2/3)), difabs(ar),difabs(amp)];
    %     outside=bar(xAX,A,'grouped');
    [AX,H1,H2] = plotyy(xAX,A,xAX,IQ,'bar','plot');
    legend('volume^{(2/3)}','area','amp','IQ','location','SouthEast');
    alld=A(:);
    my=[ceil(10*min((alld)))/10 floor(10*max((alld)))/10];
    yt=linspace(my(1),my(2),5);
    %     axis([xAX(1) xAX(end) -my my]);
    set(AX(1),'xlim',[xAX(1) xAX(end)])
    set(AX(1),'ylim',my)
    ylab = cellfun(@(cc) sprintf('% .1f', cc) , num2cell((exp(yt))), 'uniformoutput', false);
    set(AX(1),'ytick',yt)
    set(AX(1),'yticklabel',ylab)
    set(AX(1),'xtick',xtck(xr))
    set(AX(2),'xtick',[])
    set(AX(1),'xticklabel',xtckl(xr))
    set(H1(:),'Clipping','off')
    set(H2,'Clipping','off')
    set(AX(2),'ylim',[min(IQ) max(IQ)]);
    iqYt=linspace(min(IQ), max(IQ) ,5);
    iqYtL=cellfun(@(c) sprintf('%1.1f',(exp(c))),num2cell(iqYt+mean(log(iq))),'UniformOutput' ,false);
    set(AX(2),'ytick',iqYt) ;
    set(AX(2),'yticklabel',iqYtL);
    set(H2,'Clipping','off');
    xx=get(gca,'position');
    set(gca,'position',xx+[0 .1 0 0]);
    
    
    plot(xAXdouble,THR,'color','red')
    
    
    
end
function [OUT]=extractdeepfield(IN,fieldnameToAccess)
    field = textscan(fieldnameToAccess,'%s','Delimiter','.');
    fieldSize=size(field{1},1);
    switch fieldSize
        case 1
            OUT=extractfield(IN,fieldnameToAccess);
        case 2
            OUT=extractfield(cell2mat(extractfield(IN,field{1}{1})),field{1}{2} );
        case 3
            OUT=extractfield(cell2mat(extractfield(cell2mat(extractfield(IN,field{1}{1})),field{1}{2} )),field{1}{3});
        case 4
            OUT=extractfield(cell2mat(extractfield(cell2mat(extractfield( cell2mat(extractfield(IN,field{1}{1})),field{1}{2} )),field{1}{3})),field{1}{4});
    end
end

