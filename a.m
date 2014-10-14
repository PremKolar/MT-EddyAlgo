function uvg = a
	repinZ = @(A,z) repmat(permute(A,[3,1,2]),[z,1,1]);
	m=matfile('thread05');
	dd.y= @(in)  diff(in,1,2);
	dd.x= @(in)  diff(in,1,3);
	uvg = getuvg(m.UV.u,m.UV.v,m.dy,m.dx,dd,repinZ,m.Z);
end
function uvg=getuvg(u,v,dy,dx,dd,repinZ,z);dF
	uvg.dUdy = inxOry(dd.y(u),'y',dy,z,repinZ);
	uvg.dUdx = inxOry(dd.x(u),'x',dx,z,repinZ);
	uvg.dVdy = inxOry(dd.y(v),'y',dy,z,repinZ);
	uvg.dVdx = inxOry(dd.x(v),'x',dx,z,repinZ);
end
function out=inxOry(in,inxy,dxy,z,repinZ);dF
	denom=repinZ(dxy,z);
	if     strcmp(inxy,'y')
		out=in( :,[1:end, end], : )./denom;
	elseif strcmp(inxy,'x')
		out= in(:, :,[1:end, end])./denom;
	end
end