function conclude(DD)
 disp([' ']);
 disp([' ']);
 disp([' ']);
 disp(['SUCCESS!!!']);
 disp(['time used: '])  ;
 daysUsed=toc(DD.tic)/60/60/24;
disp(datestr(daysUsed,'dd-HH:MM:SS',0));
end