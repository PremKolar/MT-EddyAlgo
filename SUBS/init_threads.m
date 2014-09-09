function threads=init_threads(threads)
    warning('off','parallel:convenience:RunningJobsExist')
    warning('off','parallel:job:DestructionOfUnavailableJob')
    current_threads=matlabpool('size');
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