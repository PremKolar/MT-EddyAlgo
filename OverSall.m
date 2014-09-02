function OverSall
    initialise
    DR='/scratch/uni/ifmto/u300065/FINAL/smAor/';
    Fnows={'iq2';'iq6';'iq6dl';'iq6dlMr';'iq6dlMrIc';'iq8';'iq4';'ch'};
    Dnows=cellfun(@(c) ['data' c],Fnows,'uniformoutput',false);
    
    dbstop if error
    Rdata=[DR 'dataC/'];
    %%
    for ii=1:numel(Fnows)
        dnow=Dnows{ii};
        mkdirp([DR  dnow])
        fnow=Fnows{ii};
        mkdirp([DR  fnow])
    end
    %%
    for ii=1:numel(Fnows)
        fnow=Fnows{ii};
        dnow=Dnows{ii};
        catINPUT(fnow)
        if ii==1
            todoPre;
        end
        todoCore(DR,dnow,Rdata);
        %         todoPost(DR,dnow,Rdata);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function catINPUT(fnow)
    try system(['mv INPUT.m INPUTold.m']), end
    system(['more INPUTupper.m > INPUT.m'])
    system(['more INPUT' fnow '.m >> INPUT.m'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function todoPre
    S00b_prep_data
    S01_BruntVaisRossby
    S02_infer_fields
    S03_contours
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function todoCore(DR,dnow,Rdata)
    system(['rm -rf ' DR dnow '/' 'EDDIES'])
    system(['rm -rf ' DR dnow '/' 'TRACKS'])
    S04_filter_eddies
    S05_track_eddies
    system(['mv ' Rdata 'EDDIES ' DR dnow '/'])
    system(['mv ' Rdata 'TRACKS ' DR dnow '/'])
    cutdir=[Rdata 'CUTS/'];
    cut=dir([cutdir 'CUT_*.mat']);
    mkdirp([DR dnow '/CUTS'])
    system(['cp ' cutdir cut(1).name ' ' DR dnow '/CUTS/'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function todoPost(DR,dnow,Rdata)
    S04b_analyzeEddyThresh
    S06_init_output_maps
    S08_analyze_tracks
    S09_drawPlots
    system(['mv ' Rdata 'ANALYZED ' DR dnow '/'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

