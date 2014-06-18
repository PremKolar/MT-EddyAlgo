%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 18-Jun-2014 12:00:00
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=NCoverwriteornot(nc_file_name,ow)
	%% test if file exists already
	out=false;
	if nargin < 2
		ask(nc_file_name)
	elseif ow
		out=true;
		notask(nc_file_name)
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ask(nc_file_name)
	try
		nc_create_empty(nc_file_name,'noclobber');
	catch me
		disp(me.message)
		reply = input('Do you want to overwrite? Y/N [Y]: ', 's');
		if isempty(reply)
			reply = 'Y';
		end
		if strcmp(reply,'Y')
			nc_create_empty(nc_file_name,'clobber');
		else
			error('exiting');
		end
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function notask(nc_file_name)
	nc_create_empty(nc_file_name,'clobber');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%