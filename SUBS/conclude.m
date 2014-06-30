function conclude(DD)
    allData=dir(DD.path.root);
    allData =allData(cat(1,allData.isdir));
    [~,lastIdx]=max(cat(1,allData.datenum));
    relevantDir=allData(lastIdx);
    relevantDir.fullfile=fullfile(DD.path.root,  relevantDir.name);
    relevantDir.what=what( relevantDir.fullfile);
    %%
    clc
    disp([' ']);
    disp([' ']);
    disp([' ']);
    disp(['SUCCESS!!!']);
    disp([' ']);
    disp(['Step ' DD.monitor.rootFunc.function]);
    disp(['at ' DD.monitor.rootFunc.file]);
    disp(['has produced:' ]);
    disp([' ']);
    sleep(1);
    disp([relevantDir]);
    disp([relevantDir.what]);
    disp([' ']);
    sleep(1);
    disp([' ']);
    disp([num2str(numel(relevantDir.what.mat)), ' files like:']);
    disp(relevantDir.what.mat(1));
    disp([' ']);
    sleep(1);
    disp(['time used: '])  ;
    daysUsed=toc(DD.tic)/60/60/24;
    disp(datestr(daysUsed,'dd-HH:MM:SS',0));
     disp([' ']);   
    save_info(DD);
end