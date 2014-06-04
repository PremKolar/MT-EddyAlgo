function MAP=map_vars
	%% user input
  MAP.geo.west=-180;
  MAP.geo.east=180;
  MAP.geo.south=-90;
  MAP.geo.north=90;
    MAP.time.delta_t = 1; % [days]
	MAP.SSH_unitFactor = 100; % eg 100 if SSH data in cm, 1/10 if in deka m etc..
MAP.pattern.in='SSH_GLB_t.t0.1_42l_CORE.yyyymmdd.nc';
    MAP.pattern.lat='U_LAT_2D';
    MAP.pattern.lon='U_LON_2D';
    MAP.pattern.ssh='SSH';
  end

