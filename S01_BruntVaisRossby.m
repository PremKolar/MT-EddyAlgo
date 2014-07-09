%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 20-Jun-2014 16:53:06
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S01_BruntVaisRossby
    %% init
    DD=initialise([],mfilename);
    if ~DD.switchs.RossbyStuff,return;end
    switch DD.parameters.Nknown
        case false
            S02b_fromTS
        case true
            S02b_fromRaw
    end
end