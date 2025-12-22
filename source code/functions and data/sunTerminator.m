function [lonTerm, latTerm] = sunTerminator(subsolar_lon, subsolar_lat)

    lonTerm = 180:-1:-180;
    H = deg2rad(lonTerm - subsolar_lon);      % hour angle
    dec = deg2rad(subsolar_lat);

    % Terminator latitude formula:
    latTerm = rad2deg( atan( -cos(H) ./ tan(dec) ) );
end
