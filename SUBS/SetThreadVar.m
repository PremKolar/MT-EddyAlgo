function [OUT]=SetThreadVar(IN)
    from=IN.threads.lims(labindex,1);
    till=IN.threads.lims(labindex,2);
    num=till-from+1;
    [OUT(1:num).daynums]=IN.checks.passed(from:till).daynums;
    [OUT(1:num).files]=IN.checks.passed(from:till).filenames;
end