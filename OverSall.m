initialise
dataroot='/scratch/uni/ifmto/u300065/FINAL/aorStuff/dataAll/';
% dnows={'iq2';'iq4';'iq6';'iq8';'iq5';'iq5nonVoA';'ch400amparea'};

% dnows={'ch400amparea';'iq5nonVoA'};
% dnows={'iq5nonVoA'};
dnows={'iq5fortnight'};
% 
% for ii=1:numel(dnows)
%     dnow=dnows{ii};
%    system(['cp INPUT.m INPUT' dnow '.m'])
% end


for ii=1:numel(dnows)
    dnow=dnows{ii};
    sleep(5);
    todo=sprintf('cp INPUT%s.m INPUT.m',dnow);
    disp(todo);
    sleep(5);
    system(todo);
    sleep(5);
    Sall
    sleep(5);
    todo=sprintf('mkdir -p  %s%s',dataroot,dnow);
    disp(todo);
    system(todo);
    sleep(5);
    todo=sprintf('mv %sEDDIES %sTRACKS %sANALYZED %scode  %s%s/',dataroot,dataroot,dataroot,dataroot,dataroot,dnow);
    disp(todo);
    system(todo);
end
