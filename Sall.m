%while true
%	try
		
%         S00b_prep_data
%         s00c_fillCorruptCuts
% 		S01_BruntVaisRossby
		S04_filter_eddies
		%%
		S05_track_eddies
		S06_init_output_maps
		S08_analyze_tracks
		S09_plotsNew
		
%	catch me
%		disp(me.message)
%		disp(' waiting 5')
%		sleep(5)
%		continue
%	end
%	break
%end
