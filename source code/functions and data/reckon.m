function [latOut, lonOut] = reckon(lat1, lon1, rng, az)
    % A simplified local implementation of the reckon function if Mapping Toolbox is not available.
    % All inputs/outputs are in degrees.
    lat1 = deg2rad(lat1);
    lon1 = deg2rad(lon1);
    rng = deg2rad(rng);
    az = deg2rad(az);

    latOut = asin(sin(lat1) .* cos(rng) + cos(lat1) .* sin(rng) .* cos(az));
    lonOut = lon1 + atan2(sin(az) .* sin(rng) .* cos(lat1), cos(rng) - sin(lat1) .* sin(latOut));
    
    latOut = rad2deg(latOut);
    lonOut = rad2deg(lonOut);
end