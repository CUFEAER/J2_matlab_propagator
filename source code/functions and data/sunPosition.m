function [sun_ra, sun_dec] = sunPosition(datetimeUTC)

    % Convert to Julian Day
    if ~isnumeric(datetimeUTC)
    JD = juliandate(datetimeUTC);
    end
    T  = (JD - 2451545.0) / 36525;

    % Mean longitude of Sun (deg)
    L0 = 280.46646 + 36000.76983*T + 0.0003032*T.^2;

    % Mean anomaly of Sun (deg)
    M  = 357.52911 + 35999.05029*T - 0.0001537*T.^2;

    % Ecliptic longitude (deg)
    lambda = L0 + 1.914602*sind(M) + 0.019993*sind(2*M) + 0.000289*sind(3*M);

    % Obliquity of ecliptic (deg)
    eps = 23.439291 - 0.0130042*T;

    % Sun declination
    sun_dec = asind( sind(eps)*sind(lambda) );

    % Sun right ascension
    sun_ra = atan2d( cosd(eps)*sind(lambda), cosd(lambda) );
    sun_ra = mod(sun_ra, 360);
end
