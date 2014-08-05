function [out,file]=read_fields(dd,ff,fld,field)
	file=[dd.path.(fld).name	, dd.path.(fld).files(ff).name];
	if nargin < 4
	out=load(file);
	else
		out=load(file,field);
	end
