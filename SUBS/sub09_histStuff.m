function sub09_histStuff
    load S09main DD
    
    T(2).files=struct;
    T(1).files = dir2([DD.path.tracks.name '*' DD.FieldKeys.senses{1} '*']);
    T(2).files = dir2([DD.path.tracks.name '*' DD.FieldKeys.senses{2} '*']);
    
    for s=1:2
        N = numel(T(s).files);
        T(s).age=nan(1,N);
        for n=1:N
            T(s).age(n) = str2double(T(s).files(n).name(35:38));
        end
        
        uniO = unique(cat(2,T(:).age));
        uni = min(uniO):DD.time.delta_t*2:600;
        H{s} = histc(T(s).age,uni);
    end
    %%
    clf
    HI=reshape(cell2mat(H),numel(H{1}),[]);
    CS=cumsum(HI);
    [ax,p1,p2] = plotyy(uni,HI,uni,sum(HI,2),'bar','semilogy');
    set(p1,'BarLayout','stacked');
    legend(['anti-cyclones (' num2str(CS(end,1)) ' total)'],['cyclones (' num2str(CS(end,2)) ' total)'],'total')
    axis(ax(:),'tight')
%     xlabel('age [d]')
%     ylabel('count [1000]')
    set(ax(1),'xtick',[min(uniO) 100 200 400])
    set(ax(1),'ytick',[  3000 10000 20000])
    set(ax(1),'yticklabel',get(gca,'ytick')/1000)
    set(ax(2),'ytick',logspace(1,4,4))
    grid on
    %%
    fn=['histTrackCount'];
    savefig(DD.path.plots,200,800,500,fn,'dpdf',DD2info(DD));
    
    cpPdfTotexMT(fn);
end