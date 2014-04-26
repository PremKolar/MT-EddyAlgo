function mkdirp(D)
	if ~exist(D,'dir')
		mkdir(D)
% 	else
% 		disp([D, ' already exists. skipping..'])
	end
end