function OverSall
    initialise
    DR='/scratch/uni/ifmto/u300065/FINAL/smAor/';
    Fnows={'iq2';'iq6';'iq6dl';'iq6dlMr';'iq6dlMrIc';'iq8';'iq4';'ch'};
    Dnows=cellfun(@(c) ['data' c],Fnows,'uniformoutput',false);
    
    Rdata=[DR 'dataB/'];
    
    for ii=1:numel(Fnows)
        dnow=Dnows{ii};
        mkdirp([DR  dnow])
        fnow=Fnows{ii};
        mkdirp([DR  fnow])
    end
    %%
    
    %%
    for ii=1:numel(Fnows)
        fnow=Fnows{ii};
        dnow=Dnows{ii};
        %         if ii==1
        %             todoPre(fnow) ;
        %         end
        todo(fnow,DR,dnow,Rdata)
        try todoPost(fnow,DR,dnow,Rdata);end %#ok<*TRYNC>
    end
end

function todoPre(fnow)
    todo=sprintf('cp INPUT%s.m INPUT.m',fnow);
    disp(todo);
    system(todo);
    S00b_prep_data
    S01_BruntVaisRossby
    S02_infer_fields
    S03_contours
end

function todo(fnow,DR,dnow,Rdata)
    todo=sprintf('cp INPUT%s.m INPUT.m',fnow);
    disp(todo);
    system(todo);
    
    S04_filter_eddies
    S05_track_eddies
    
    system(['rm -rf ' DR dnow '/' 'EDDIES'])
    system(['rm -rf ' DR dnow '/' 'TRACKS'])
    
    system(['mv ' Rdata 'EDDIES ' DR dnow '/'])
    system(['mv ' Rdata 'TRACKS ' DR dnow '/'])
    cutdir=[Rdata 'CUTS/'];
    cut=dir([cutdir 'CUT_*.mat']);
    mkdirp([DR dnow '/CUTS'])
    system(['cp ' cutdir cut(1).name ' ' DR dnow '/CUTS/'])
    
end



function todoPost(fnow,DR,dnow,Rdata)
    todo=sprintf('cp INPUT%s.m INPUT.m',fnow);
    disp(todo);
    system(todo);
    S04b_analyzeEddyThresh
    S06_init_output_maps
    S08_analyze_tracks
    S09_drawPlots
    
    %     S10_makeAnimations
    
    system(['mv ' Rdata 'ANALYZED ' DR dnow '/'])
    
    
end

