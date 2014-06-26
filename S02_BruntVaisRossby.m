%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 20-Jun-2014 16:53:06
% Computer:  GLNX86
% Matlab:  7.9
% Author:  NK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S02_BruntVaisRossby
    %% init
    DD=initialise;
    switch DD.parameters.Nknown
        case false
            S01b_fromTS
        case true
            S01b_fromRaw
    end
    conclude(DD);
end