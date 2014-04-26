function save_info(DD)	
%% refresh
DDtemp=get_input; 
DD.path=DDtemp.path;
%% save
save([DD.path.root, 'DD.mat'],'-struct',	'DD')