%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 18-Jun-2014 12:00:00
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NCoverwriteornot(nc_file_name)
	%% test if file exists already
	GiveTime(nc_file_name);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GiveTime(nc_file_name)
	try
		nc_create_empty(nc_file_name,'noclobber');
	catch me
		disp(me.message)
		tt= 1;
		while tt>0
			disp('cancel if you dont want to overwrite!')
			disp(datestr(sec2day(tt),'SS'))
			tt=tt-1;
			sleep(1)
		end
		disp('overwriting')
		nc_create_empty(nc_file_name,'clobber');
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
