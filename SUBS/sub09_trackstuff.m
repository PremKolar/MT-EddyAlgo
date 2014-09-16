function sub09_trackstuff
    load S09main II DD T    
    try
        load trackinit
    catch me
        disp(me.message)
        trackinit(DD);
    end
    %%
    senses=DD.FieldKeys.senses;
    catsen=@(f) [TR.(senses{1}).cats.(f) TR.(senses{2}).cats.(f) ];
    t2l=@(t) linspace(t(1),t(2),t(3));
    
    %%
    rad=catsen('radiusmean')/1000;
    U  =catsen('U')*100;
    age=catsen('age');
    lat=abs(catsen('lat'));
    %%
    rightyscalenum=5;
    age(end+1:end+rightyscalenum)=max(age)-0;
    lat(end+1:end+rightyscalenum)=t2l([min(lat) max(lat) rightyscalenum]);
    rad(end+1:end+rightyscalenum)=t2l([min(rad) max(rad) rightyscalenum]);
    U(end+1:end+rightyscalenum)=10;
    %%
    [~,sml2lrg] = (sort(rad))  ;
    age=age(fliplr(sml2lrg));
    lat=lat(fliplr(sml2lrg));
    rad=rad(fliplr(sml2lrg));
    U=U(fliplr(sml2lrg));
    %%
    zerage = age<=0  ;
    age(zerage)=[];
    lat(zerage)=[];
    rad(zerage)=[];
    U(zerage)=[];
    %%
    clf
    hs=scatter(age,lat,rad,U);
    axis tight
    set(gca,'XAxisLocation','bottom')
    set(gca,...
        'ytick',t2l(T.lat),...
        'xtick',t2l(T.age));
    cb=colorbar;
    cb1 = findobj(gcf,'Type','axes','Tag','Colorbar');
    cbIm = findobj(cb1,'Type','image');
    alpha(cbIm,0.5)
    set(cb,'location','north','xtick',t2l(T.vel),'xlim',T.vel([1 2]))
    doublemap([T.vel(1) 0 T.vel(2)],II.aut,II.win,[.9 1 .9],20)
    h1=gca;
    h1pos = get(h1,'Position'); % store position of first axes
    h2 = axes('Position',h1pos,...
        'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none',...
        'ytick',linspace(0,1,rightyscalenum),...
        'xtick',[],...
        'yticklabel',round(t2l([min(rad) max(rad) rightyscalenum])));
    ylabel(h2,'radius [km]')
    ylabel(h1,'lat  [{\circ}]')
    xlabel(h1,'age [d]')
    xlabel(h2,'zon. vel.  [cm/s]')
    
    set(get(gcf,'children'),'clipping','off')
    set(get(hs,'children'),'clipping','off')
    
    
    savefig(DD.path.plots,T.rez,T.width,T.height,['sct-ageLatRadU'],'dpdf');
    
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trackinit(DD)
    tsenses=fieldnames(DD.path.analyzedTracks)';
    senses=DD.FieldKeys.senses;
    lims2array=@(lims) lims(labindex,1):lims(labindex,2);
    for ss=1:2
        clear single cats
        tsen=tsenses{ss};
        sen = senses{ss};
        root=DD.path.analyzedTracks.(tsen).name;
        eds= DD.path.analyzedTracks.(tsen).files;
        
        JJ=thread_distro(DD.threads.num,numel(eds));
        spmd
            T=disp_progress('init','blubb');
            for ff= lims2array(JJ)
                T=disp_progress('calc',T,diff(JJ(labindex,:))+1);
                single(ff)=load([root eds(ff).name]);
            end
            single=gcat(single,2,1);
        end
        single=single{1};
        
        for fn=fieldnames(single)'; fn=fn{1};
            cats.(fn)=extractdeepfield(single,fn);
        end
    end
    TR.(sen).single=single;
    TR.(sen).cats=cats;
    
    save trackinit TR
end
