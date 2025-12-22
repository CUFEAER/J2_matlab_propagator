

function [LMST_hours, LAST_hours] = getLST(utc_dt, longitude_deg)
% getLocalSolarTime  Return local mean & apparent solar time
%   Inputs:
%     utc_dt        - MATLAB datetime in UTC (timezone must be 'UTC') OR numeric Julian day (see note)
%     longitude_deg - longitude in degrees (east positive)
%   Outputs:
%     LMST_hours  - Local Mean Solar Time in hours [0,24)
%     LAST_hours  - Local Apparent Solar Time (includes Equation of Time) [0,24)
%     EoT_min     - Equation of Time in minutes (approximation)
%
% Example:
%   dtUTC = datetime(2025,11,29,12,0,0,'TimeZone','UTC');
%   [LMST, LAST, EoT] = getLocalSolarTime(dtUTC, 31.2); % Cairo ~31.2E

    % If user passed numeric (Julian Day) convert to datetime UTC
    if isnumeric(utc_dt)
        % assume utc_dt is Julian day (JD)
        JD = utc_dt;
        % convert JD to datetime (UTC)
        % MATLAB: datetime( 'convert from JD' ) - implement simple conversion:
        % JD to modified julian then to datetime:
        MJD = JD - 2400000.5;
        % base date for MJD = 1858-11-17
        utc_dt = datetime(1858,11,17,0,0,0,'TimeZone','UTC') + days(MJD);
    end

    % ensure tz is UTC
    if isempty(utc_dt.TimeZone) || ~strcmp(utc_dt.TimeZone,'UTC')
        utc_dt.TimeZone = 'UTC';
    end

    % 1) UTC in hours (fractional)
    utc_hour = hours(timeofday(utc_dt));

    % 2) Local Mean Solar Time (hours)
    LMST_hours = utc_hour + longitude_deg/15;  % east positive: +longitude/15
    LMST_hours = mod(LMST_hours, 24);

    % 3) Equation of Time (approx) - needs day of year (N)
    N = day(utc_dt,'dayofyear');
    B_deg = (360/365) * (N - 81);
    B = deg2rad(B_deg);

    % EoT in minutes (approximation)
    EoT_min = 9.87 * sin(2*B) - 7.53 * cos(B) - 1.5 * sin(B);

    % 4) Local Apparent Solar Time (hours)
    LAST_hours = LMST_hours + EoT_min/60;
    LAST_hours = mod(LAST_hours, 24);
end
