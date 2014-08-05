

function ampOra
    addpath(genpath('./'))
   dbstop if error 
    
    savestuff('trackiq2')
    savestuff('trackiq4')
    savestuff('trackiq6')
    savestuff('trackiq8')
    savestuff('trackch400')
    
    
    system(['pdfjam -o tmp.pdf crpd*pdf'])
    system(['pdfcrop  --margins "1 3 1 1" tmp.pdf jmd.pdf'])
    system(['okular jmd.pdf'])
end



function savestuff(in)
    pl= @(ss) AOplots(load(ss));
    pl(in)
    %     system(['rm track*.pdf'])
   
        set(gcf,'position',[5 40 1912 437])
        
        
        title(in)
        savefig('./',100,1600,500,['tmp'])        
       system(['pdfcrop --margins "1 3 1 1" tmp.pdf crpd' in '.pdf'])
      
    
    
end





function AOplots(trck)
    trck=trck.trck;
    ar=extractfield(cell2mat(extractfield(trck,'area')),'total');
    ra=sqrt(ar/pi);
    amp=extractfield(cell2mat(extractfield(cell2mat(extractfield(trck,'peak')),'amp')),'to_contour')
    age=cat(2,trck.age);
    iq=cat(2,trck.isoper);
    
    %%
    
    nrm=@(x) (x-min(x))/max(x-min(x))
    
    figure(1)
    clf
    subplot(2,1,1,'align')
 
    hold on
%     qu=(amp./ra);
%    Ro=pi*1e5;    
%    qu=amp./(1-exp(-(1/Ro)*ra));
vol=    extractdeepfield(trck,'volume.total') ;
 %   sdfhdsgfhrde
  qu=((vol.^(2/3)./ar).^2)
  sleep(5) 
    quo=([1 qu(2:end)./qu(1:end-1)]);
   quo=log(quo); 
    
    col=hsv;
    
    
    
    col=col(1:end-5,:);
    col=resample(col,numel(quo),size(col,1));
    col(col<0)=0;
    col(col>1)=1;
    quon=nrm(quo)*(size(col,1)-1)+1;
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
    cbtckl=cellfun(@(cc) sprintf('% .3f', cc) ,num2cell((exp(ct))),'uniformoutput',false);
    cbtckl{3}=1;
    set(cb,'yticklabel',cbtckl)
    
    cm1=(bone);
    cm2=(winter);
    cm1=(autumn)
    cm2=(summer)
    
    
    CM=doublemap(cb,cm1,cm2,cm2(end,:));
    
    CMi=round(linspace(1,size(CM,1),numel(trck)));
    
    
    for ii=1:1:numel(trck)
        y=extractdeepfield(trck(ii),'coordinates.exact.y');
        x=extractdeepfield(trck(ii),'coordinates.exact.x')-ii/numel(trck)*1000;
        plot(x,y,'color',CM(CMi(round(quon(ii))),:))
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
    hold on
   
    difabs= @(a) fliplr(log([1 a(2:end)./a(1:end-1)]));
    
   
    a= [difabs(vol); difabs(ar)]
    bar(xAX,a','grouped')
    legend('volume','area','location','SouthEast')
    
    alld=[difabs(ra) difabs(amp)]
    my=floor(10*max(abs(alld)))/10;
    yt=-my:2*my/4:my
    
    axis([xAX(1) xAX(end) -my my])
    
    ylab = cellfun(@(cc) sprintf('% .3f', cc) , num2cell((exp(yt))), 'uniformoutput', false)
    set(gca,'ytick',yt)
    set(gca,'yticklabel',ylab)
    %    ylabel(['1e-6 radius / amplitude / iq'])
    
    
   
    set(gca,'xtick',xtck(xr))
    set(gca,'xticklabel',xtckl(xr))
    
    hold on
    
    iqs=sqrt(fliplr([(iq-min(iq))/(max(iq)-min(iq))]))*2*my-my;
    
    plot(xAX,iqs)
    
 
     AX=get(gca,'position')
    
     set(gca,'position',AX+[0 .1 0 0])
    
    
    
    
    
    
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

