function [epoch, ndot, i, W0, e, w0, M, n, bstar, epoch_date_formatted] = readTLE(filename)
    % READTLE Parses a standard 2-line element file
    fid = fopen(filename,'r');
    fgetl(fid); % Skip name
    line1 = fgetl(fid);
    line2 = fgetl(fid);
    fclose(fid);
    
    epoch = str2double(line1(19:32));
    ndot  = str2double(line1(34:43));
    % Handle exponential notation manually for bstar
    mant = line1(54:59);
    expPart = line1(60:61);
    bstar = str2double([mant(1) '.' mant(2:end) 'e' expPart]);
    
    i     = str2double(line2(9:16));
    W0    = str2double(line2(18:25));
    e     = str2double(['0.' line2(27:33)]);
    w0    = str2double(line2(35:42));
    M     = str2double(line2(44:51));
    n     = str2double(line2(53:63));
    % === Convert TLE epoch (YYDDD.DDDDD) to datetime ======================
epoch_year_short = floor(epoch / 1000);           % YY
epoch_day        = epoch - epoch_year_short*1000; % DDD.ddddd

% Convert 2-digit year to 4 digits
if epoch_year_short < 57
    epoch_year = 2000 + epoch_year_short;
else
    epoch_year = 1900 + epoch_year_short;
end

epoch_datetime = datetime(epoch_year,1,1) + days(epoch_day - 1);

% === Print epoch BEFORE showing options ===============================
fprintf("\n==========================================\n");
fprintf("   TLE Epoch Time (UTC): %s\n", datestr(epoch_datetime,'yyyy-mm-dd HH:MM:SS'));
fprintf("==========================================\n\n");
epoch_date_formatted=datestr(epoch_datetime,'yyyy-mm-dd HH:MM:SS');
end