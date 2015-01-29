function sub09_histStuff
    load DD
    
    T(2).files=struct;
    T(1).files = dir2([DD.path.tracks.name '*' DD.FieldKeys.senses{1} '*']);
    T(2).files = dir2([DD.path.tracks.name '*' DD.FieldKeys.senses{2} '*']);
    
    for s=1:2
        N = numel(T(s).files);
        T(s).age=nan(1,N);
        for n=1:N
            T(s).age(n) = str2double(T(s).files(n).name(35:38));
        end
        
        uni = unique(T(s).age);
        H{s} = histc(T(s).age,uni);
    end
    %%
    bar(uni,reshape(cell2mat(H),numel(H{1}),[]),'stacked')
    legend(['anti-cyclones (' num2str(sum(H{1})) ' total)'],['cyclones (' num2str(sum(H{2})) ' total)'])
    axis tight
    xlabel('age [d]')
    ylabel('count')
    set(gca,'xtick',round(linspace(min(uni),max(uni),5)))
    set(gca,'ytick',round(linspace(min(sum(cat(1,H{:}))),max(sum(cat(1,H{:}))),5)))
    %%
    fn=['histTrackCount'];
    savefig(DD.path.plots,200,800,500,fn,'dpdf',DD2info(DD));
    
    cpPdfTotexMT(fn);
end