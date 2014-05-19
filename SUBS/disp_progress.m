function [T]=disp_progress(type,Tin,L,num_prints)
    warning('off','MATLAB:divideByZero')
    if strcmp(type,'init')
        T=init(Tin);
    else
        T=later(Tin,L,num_prints);
    end
    warning('on','MATLAB:divideByZero') %#ok<*RMWRN>
end

function T=init(Tin)
    T.cc=0;
    T.time=0;
    T.name=Tin;
    T.tic=tic;
end
function T=later(Tin,L,num_prints)
    T=Tin;
    T.cc=T.cc+1;
    l=T.cc;
    %%
    T=calcu(T,l,L);
    if mod(l,ceil(L/num_prints))==0;
        %%
        printout(T,l,L);
    end
end
function printout(T,l,L)
    disp('####')
    disp(['-',T.name,'-']);
    disp('####')
    disp(['step: ',num2str(l),'/',num2str(L)]);
    disp([ num2str(round(T.prcnt_done)),' %']);
    disp(['time so far:   ', datestr(T.time/86400,'HH:MM:SS.FFF')]);
    if isfinite(T.time_to_go)
        disp(['time to go  :    ', datestr(T.time_to_go/86400,'HH:MM:SS.FFF')]);
        spmdwaitbar(l/L,30);
    else
        disp(['time to go:    ', 'calculating...']);
    end
    try
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

function T=calcu(T,l,L)
    T.time=toc(T.tic);
    T.prcnt_done=((l-1)/L);
    T.time_to_go=T.time/T.prcnt_done-T.time;
    T.prcnt_done=T.prcnt_done*100;
    try %#ok<*TRYNC>
        T.uplevel.full_time_to_go=(T.time_to_go + T.time)*T.uplevel.togo;
    end
end



