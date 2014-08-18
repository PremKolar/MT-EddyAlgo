%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 18-Oct-2013 11:58:41
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T]=disp_progress(type,Tin,L,num_prints)
    if nargin==3
        num_prints=L;
    end
    if strcmp(type,'conclude')
        conclude; T=[]; return
    end
    %%
    if labindex == 1
        switch type
            case 'init'
                T=init(Tin);return;
            otherwise
                T=later(Tin,L,num_prints);return;
        end
    else
        T=[];return
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function conclude
    %% wait for all workers
    labBarrier
    %% let master handle I/O
    if labindex == 1
        c=initc;
        for cc=2:c.num
            echoHis(c,cc);
        end
    end
    %% wait for all workers
    labBarrier
    %----------------------------------------------------------------------
    function c=initc
        c.files=dir('./.comm*');
        c.num=numel(c.files);
        sprintf('...going to read out all old recently collected mails from all threads(id>1)...\n' )
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function echoHis(c,cc)
    d=matfile(c.files(cc).name,'writable',true);
    tit=cell2mat(d.printstack(1,1));
    sprintf('reading disp stack for %s:\n\n', tit);
    for ii=2:numel(d.printstack)
        tot=cell2mat(d.printstack(ii,1));
        %% echo
        if numel(tot)>12
            cprintf(rainbow(1,1,1,cc-1,c.num-1),[tit ':\n' tot ' \n']); disp('')
        end
        %% recycle
        d.printstack(ii,1)={[]};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T=init(Tin)
    T.cc=0;
    T.time=0;
    T.name=Tin;
    T.tic=tic;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T=later(T,L,num_prints)
    T.cc=T.cc+1;
    %%
    if mod(T.cc,ceil(L/num_prints))==0;
        T=calcu(T,L);
        printout(T,L);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printout(T,L)
    %% init
    incol.a=[1 0 0];
    incol.b=[0 0 1];
    fracCol=T.frac*incol.a + (1-T.frac)*incol.b;
    %% build output
    strout=makeStrings;
    %% print
%     clc;
    sendoutStrings;
    %% waitbar
    spmdwaitbar(T.cc,L,30);
    %% SUBS
    %----------------------------------------------------------------------
    function sendoutStrings
        for a=1:numel( strout.a)
            disp(strout.a{a})
        end
        for b=1:numel( strout.b)
            cprintf(fracCol,strout.b{b})
        end
    end
    %----------------------------------------------------------------------
    function strout=makeStrings
        strout.a{1}='';
        strout.a{2}=['master thread says:'];
        strout.a{3}='####';
        strout.a{4}=['-',T.name,'-'];
        strout.a{5}='####';
        %%
        strout.b{1}=['step: ',num2str(T.cc),'/',num2str(L),'\n'];
        strout.b{2}=[ num2str(round(T.prcnt_done)),' percent done.\n'];
        strout.b{3}=['time so far:   ', datestr(sec2day(T.time),'dd-HH:MM:SS'),'\n'];
        if isfinite(T.time_to_go)
            strout.b{4}  = ['time to go     :    ', datestr(sec2day(T.time_to_go),'dd-HH:MM:SS'),'\n'];
        else
            strout.b{4}  = ['time to go     :    ', 'calculating...\n'];
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=spmdwaitbar(l,L,len)
    frac=l/L;
    out=['[',repmat('-',1,floor(frac*len)),'>',repmat(' ',1,ceil((1-frac)*len)),']\n'];
    cprintf(rainbow(1,1,1,l,L), out);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T=calcu(T,L)
    T.time=toc(T.tic);
    T.frac=((T.cc-1)/L);
    T.time_to_go=T.time/T.frac-T.time;
    T.prcnt_done=T.frac*100;
    if isfield(T,'uplevel')
        T.uplevel.full_time_to_go=(T.time_to_go + T.time)*T.uplevel.togo;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



