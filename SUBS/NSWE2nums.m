function out=NSWE2nums(pathout,filename,geo,timestr)    
    out=strrep(filename,'SSSS',sprintf('%+03.0f',geo.south) ); %#ok<*NASGU>
    out=strrep(out, 'NNNN',sprintf('%+03.0f',geo.north) );
    out=strrep(out, 'WWWW',sprintf('%+03.0f',geo.west) );
    out=strrep(out, 'EEEE',sprintf('%+03.0f',geo.east) );
    out=[pathout, strrep(out, 'yyyymmdd',timestr)]; 
end