function MAP=map_vars
	%% user input
  MAP.geo.west=-180;
  MAP.geo.east=180;
  MAP.geo.south=-90;
  MAP.geo.north=90;
    MAP.time.delta_t = 1; % [days]
	MAP.SSH_unitFactor = 100; % eg 100 if SSH data in cm, 1/10 if in deka m etc..
	MAP.pattern.in='SsaltoDuacs__merged_msla__AVISO__ref__0.333deg__yyyymmdd.nc';
    MAP.pattern.lat='lat';
    MAP.pattern.lon='lon';
    MAP.pattern.ssh='msla';
  end

