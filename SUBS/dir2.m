function dir2(d)
	S=dir(d);
	[base,~,~]=fileparts(S);
	for ii=1:numel(S)
		S(ii).fullname=[base S(ii).name];
	end
end