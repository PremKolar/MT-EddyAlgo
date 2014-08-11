%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 15-Jul-2014 13:52:44
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NKkk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maxOW
    %% init
    DD=initialise([],mfilename);
    main(DD);
    %% post process
    maxOWpostProc
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function main(DD)
    %% set up
    maxOWsetUp(DD);
    %% spmd
  maxOWmain
    %%
    maxOWprocess;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
