%%
function threads=init_threads(threads)
    warning('off','parallel:convenience:RunningJobsExist')
    warning('off','parallel:job:DestructionOfUnavailableJob')
    [~,vers]=version;
    versN=datenum(vers);
    switch sign(versN-datenum('2014','yyyy'))
        case -1
            oldCase(threads)
        case 1
            newCase(threads)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function oldCase(threads)
    current_threads=matlabpool('size');%#ok<*DPOOL>
    max_threads = feature('numCores');
    if max_threads<threads
        disp(['you only have ' num2str(max_threads) ' cores available']);
        threads=max_threads;
    end
    if current_threads < threads
        try
            matlabpool close
        catch err
            disp(err.message)	;
        end
        if threads>1
            matlabpool(threads);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newCase(threads)
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        current_threads = 0;
    else
        current_threads = poolobj.NumWorkers;
    end
    %%
    if current_threads < threads
        try
           delete(poolobj)
        catch err
            disp(err.message)	;
        end
        if threads>1
            parpool(threads);
        end
    end
    
end