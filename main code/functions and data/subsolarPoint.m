function [subsolar_lon, subsolar_lat] = subsolarPoint(tUTC)
    if ~isnumeric(tUTC)
        JD = juliandate(tUTC);
    else 
        JD = tUTC;
    end
    T  = (JD - 2451545.0) / 36525;

    % Mean longitude
    L0 = 280.46646 + 36000.76983*T;

    % Mean anomaly
    M = 357.52911 + 35999.05029*T;

    % Ecliptic longitude
    lambda = L0 + 1.914602*sind(M) + 0.019993*sind(2*M);

    % Obliquity
    eps = 23.439291 - 0.0130042*T;

    % Declination
    subsolar_lat = asind(sind(eps).*sind(lambda));

    % Right ascension
    RA = atan2d(cosd(eps).*sind(lambda), cosd(lambda));

    % GMST
    D = JD - 2451545;
    GMST = mod(18.697374558 + 24.06570982441908*D, 24);

    % Subsolar longitude
    subsolar_lon = mod((RA/15 - GMST)*15, 360);
    if subsolar_lon > 180
        subsolar_lon = subsolar_lon - 360;
    end
end
