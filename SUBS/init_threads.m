function []=init_threads(threads)
	current_threads=matlabpool('size');
	if current_threads~=threads
		try
			matlabpool close
		catch err
			disp(err.message)	;
		end
		if threads>1
			matlabpool(threads);
		end
	end