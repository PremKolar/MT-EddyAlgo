 function cpPdfTotexMT(fn)
        [roo,dhere]=fileparts(pwd);
        ROO = fileparts(roo);
        prnt = [roo '/PLOTS/' dhere '/' fn '.pdf'];            
        c= system(['cp ' prnt ' ' ROO '/texMT/FIGS/' fn '-' dhere  '.pdf']);
        if c, error('yo'),end
    end