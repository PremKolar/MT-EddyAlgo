function conclude(DD,void)  %#ok<INUSD>
   main;
    %----------------------------------------------------------------------
    function main
         if ~exist('void','var'),  output;    end
    save_info(DD);
    delete ./.comm*.mat
    end
    %----------------------------------------------------------------------
    function output
        allData=dir(DD.path.root);
        allData =allData(cat(1,allData.isdir));
        allData =allData(3:end);
        [~,lastIdx]=max(cat(1,allData.datenum));
        relevantDir=allData(lastIdx);
        relevantDir.fullfile=fullfile(DD.path.root,  relevantDir.name);
        relevantDir.what=what( relevantDir.fullfile);
        %%
        inform
    end
    %----------------------------------------------------------------------
    function inform
        
        disp([' ']);
        disp([' ']);
        disp([' ']);     
        disp([mfilename ' - SUCCESS!!!']);
        disp([' ']);
        try %#ok<TRYNC>
            tobetried
        end
        disp(['time used: '])  ;
        daysUsed=toc(DD.monitor.tic)/60/60/24;
        disp(datestr(daysUsed,'dd-HH:MM:SS',0));
        disp([' ']);
        % . . . . 
        function tobetried
            disp(['Step ' DD.monitor.rootFunc.function]);
            disp(['at ' DD.monitor.rootFunc.file]);
            disp(['has produced:' ]);
            disp([' ']);           
            disp([relevantDir]);
            disp([relevantDir.what]);
            disp([' ']);           
            disp([' ']);
            disp([num2str(numel(relevantDir.what.mat)), ' files like:']);
            disp(relevantDir.what.mat(1));
            disp([' ']);           
        end
        
        
    end
end

