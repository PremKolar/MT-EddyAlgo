function dispM(string)
    if labindex==1
        disp(string);
    else
        commFile=['./.comm' num2str(labindex) '.mat'];
        comm=matfile(commFile,'writable',true);
        comm.printstack(end+1,1)={string}; 
    end
end