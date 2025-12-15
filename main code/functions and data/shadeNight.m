function h = shadeNight(lonTerm, latTerm, subsolar_lon, subsolar_lat)
    % SHADENIGHT - Ultra-simple physically correct shading
    % Rule: Shade the hemisphere OPPOSITE the Sun
    
    % Sort terminator
    [lon_sorted, idx] = sort(lonTerm);
    lat_sorted = latTerm(idx);
    
    % Determine which pole is in darkness
    % If Sun is in NORTHERN hemisphere (positive lat):
    %   → South Pole is in polar night → shade SOUTHERN hemisphere
    % If Sun is in SOUTHERN hemisphere (negative lat):
    %   → North Pole is in polar night → shade NORTHERN hemisphere
    
    if subsolar_lat > 0
        % Northern Summer: Sun north → shade SOUTH
        pole_to_shade = -90;  % South Pole
        fprintf('Sun at %.1f°N → shading SOUTHERN hemisphere\n', subsolar_lat);
    else
        % Northern Winter: Sun south → shade NORTH
        pole_to_shade = 90;   % North Pole
        fprintf('Sun at %.1f°S → shading NORTHERN hemisphere\n', abs(subsolar_lat));
    end
    
    % Create polygon: terminator → correct pole
    lon_poly = [lon_sorted(:); flipud(lon_sorted(:))];
    lat_poly = [lat_sorted(:); pole_to_shade * ones(size(lon_sorted(:)))];
    
    % Close polygon
    if lon_poly(1) ~= lon_poly(end)
        lon_poly(end+1) = lon_poly(1);
        lat_poly(end+1) = lat_poly(1);
    end
    
    % Draw
    h = fill(lon_poly, lat_poly, [0 0 0], ...
             'FaceAlpha', 0.25, 'EdgeColor', 'none');
    
    % Keep terminator visible
    plot(lonTerm, latTerm, 'w-', 'LineWidth', 1.5);
    
    % Add informative markers
    plot(subsolar_lon, subsolar_lat, 'yo', ...
         'MarkerSize', 12, 'MarkerFaceColor', 'y', ...
         'DisplayName', sprintf('Sun (%.1f°)', subsolar_lat));
    
    % Mark the pole in darkness
    if pole_to_shade == 90
        plot(0, 90, 'b*', 'MarkerSize', 10, ...
             'DisplayName', 'North Pole (Polar Night)');
    else
        plot(0, -90, 'b*', 'MarkerSize', 10, ...
             'DisplayName', 'South Pole (Polar Night)');
    end
end