function [T]=disp_progress(type,Tin,L,num_prints,silent)
	warning('off','MATLAB:divideByZero')
	if strcmp(type,'init')
		T=init(Tin);
	else
		if nargin<5, silent=false;end
		T=later(Tin,L,num_prints,silent);
	end
	warning('on','MATLAB:divideByZero') %#ok<*RMWRN>
end

function T=init(Tin)
	T.cc=0;
	T.time=0;
	T.name=Tin;
	T.tic=tic;
end
function T=later(T,L,num_prints,silent)
	T.cc=T.cc+1;
	%%
	if mod(T.cc,ceil(L/num_prints))==0 && ~silent;
		T=calcu(T,L);
		printout(T,L);
	end
end
function printout(T,L)
	disp('####')
	disp(['-',T.name,'-']);
	disp('####')
	disp(['step: ',num2str(T.cc),'/',num2str(L)]);
	disp([ num2str(round(T.prcnt_done)),' %']);
	disp(['time so far:   ', datestr(T.time/86400,'dd-HH:MM:SS')]);
	if isfinite(T.time_to_go)
		disp(['time to go  :    ', datestr(T.time_to_go/86400,'dd-HH:MM:SS')]);
		spmdwaitbar(T.cc/L,30);
	else
		disp(['time to go:    ', 'calculating...']);
	end
	if isfield(T,'uplevel')
		disp(['full time togo:    ', datestr(T.uplevel.full_time_to_go/86400,'dd-HH:MM:SS')]);
		disp(T.uplevel.perDone);
		spmdwaitbar(T.uplevel.l/T.uplevel.L,30);
	end
	%%
	
end

function out=spmdwaitbar(frac,len)
	out=['[',repmat('-',1,floor(frac*len)),'>',repmat(' ',1,ceil((1-frac)*len)),']'];
	disp(out);
end

function T=calcu(T,L)
	T.time=toc(T.tic);
	T.prcnt_done=((T.cc-1)/L);
	T.time_to_go=T.time/T.prcnt_done-T.time;
	T.prcnt_done=T.prcnt_done*100;
	if isfield(T,'uplevel')
		T.uplevel.full_time_to_go=(T.time_to_go + T.time)*T.uplevel.togo;
	end
end



