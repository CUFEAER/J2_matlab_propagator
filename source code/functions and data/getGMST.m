function gmst = getGMST(epoch_val)

    % 1. Extract Year 
    year_short = floor(epoch_val / 1000);
    
    % 2. Extract Day of Year 
    day = epoch_val - (year_short * 1000);
    
    % 3. Calculate Julian Date
    year = 2000 + year_short; % Assumes 21st century (2025)
    
    % Standard Astrodynamics Formula for JD at Jan 0.0
    jd_begin = 367*year - floor(7*(year+floor((10)/12))/4) + floor(275*1/9) + 1721013.5;
    
    % Add the Day of Year to get full JD
    jd = jd_begin + day;
    
    % 4. Calculate GMST (Greenwich Mean Sidereal Time)
    T_ut1 = (jd - 2451545.0) / 36525;
    
    gmst_sec = 67310.54841 + (876600*3600 + 8640184.812866)*T_ut1 ...
               + 0.093104*T_ut1^2 - 6.2e-6*T_ut1^3;
    
    % Convert seconds to degrees and normalize to 0-360
    gmst_deg = mod(gmst_sec * (360/86400), 360);
    
    % Convert to Radians for the matrix
    gmst = deg2rad(gmst_deg);
end