function save_info(DD)	
%% refresh
DDtemp=get_input; 
DD.path=DDtemp.path;
%% save
save('Sall_output.mat','-struct',	'DD')