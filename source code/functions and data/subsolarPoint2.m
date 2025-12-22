function [subsolar_lon, subsolar_lat] = subsolarPoint2(tUTC)

    if ~isnumeric(tUTC)
        JD = juliandate(tUTC);
    else
        JD = tUTC;
    end

    T = (JD - 2451545.0)/36525;

    % Sun’s mean longitude and anomaly
    L0 = mod(280.460 + 36000.770*T, 360);
    M  = mod(357.529 + 35999.050*T, 360);

    % Equation of center
    C = (1.915*sind(M)) + (0.020*sind(2*M));

    % True ecliptic longitude
    lambda = mod(L0 + C, 360);

    % Obliquity
    eps = 23.439 - 0.0000004*T;

    % Declination
    subsolar_lat = asind( sind(eps).*sind(lambda) );

    % Right ascension
    RA = atan2d( cosd(eps).*sind(lambda), cosd(lambda) );
    RA = mod(RA, 360);

    % GMST (deg)
    GMST = mod(280.46061837 + ...
               360.98564736629*(JD-2451545) + ...
               0.000387933*T.^2 - (T.^3)/38710000, 360);

    % Subsolar longitude (deg)
    subsolar_lon = RA - GMST;
    subsolar_lon = mod(subsolar_lon+180,360) ;
end
