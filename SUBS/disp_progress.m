function [T]=disp_progress(type,Tin,L,num_prints,uplevelratio)
	warning('off','MATLAB:divideByZero')
	if strcmp(type,'init')
		T=init(Tin);
	else
		if nargin<5, uplevelratio=nan; end
		T=later(Tin,L,num_prints,uplevelratio);
	end
	warning('on','MATLAB:divideByZero')
end

function T=init(Tin)
	T.cc=0;
	T.time=0;
	T.name=Tin;
	T.tic=tic;
end
function T=later(Tin,L,num_prints,uplevelratio)
	T=Tin;
	T.cc=T.cc+1;
	l=T.cc;
	if mod(l,ceil(L/num_prints))==0;
		%%
		T=calcu(T,l,L,uplevelratio);
		%%
		T.wb=printout(T,l,L);
	end
end
function wb=printout(T,l,L)
	disp('####')
	disp(['-',T.name,'-']);
	disp('####')
	disp(['step: ',num2str(l),'/',num2str(L)]);
	disp([ num2str(round(T.prcnt_done*100)/100),' %']);
	disp(['time so far:   ', datestr(T.time/86400,'HH:MM:SS')]);
	try
		disp(['time to go  :    ', datestr(T.time_to_go/86400,'HH:MM:SS')]);
		if ~isnan(T.full_time_to_go)
			disp(['full time tg:    ', datestr(T.full_time_to_go/86400,'dd-HH:MM:SS')]);
		end
	catch %#ok<CTCH>
		disp(['time to go:    ', 'calculating...']);
	end
	%%
	wb=spmdwaitbar(l/L,30);
end

function out=spmdwaitbar(frac,len)
	out=['[',repmat('-',1,floor(frac*len)),'>',repmat(' ',1,ceil((1-frac)*len)),']'];
	disp(out);
end

function T=calcu(T,l,L,ulr)
	T.toc=toc(T.tic);
	T.tic=tic;
	T.time=T.time+T.toc;
	T.prcnt_done=((l-1)/L)*100;
	T.time_to_go=T.time/T.prcnt_done*100-T.time;
	T.full_time_to_go=(T.time_to_go+T.time)/ulr;
end



