function h = createNightShading(subsolar_lon, subsolar_lat)
    % CREATENIGHTSHADING - Robust, physically correct night shading
    % This function generates a terminator great circle using 3D vector math,
    % which avoids singularities at the equinoxes. It then correctly
    % handles the +/-180 degree dateline by splitting the night polygon
    % into two patches if necessary, ensuring a perfect visual representation.
    %
    % Inputs:
    %   subsolar_lon: Sun's longitude [deg, -180 to 180]
    %   subsolar_lat: Sun's declination [deg, -23.44 to 23.44]
    %
    % Output:
    %   h: Handle(s) to the created patch object(s).

    % 1. Define the center of the night side (anti-solar point)
    night_center_lon = subsolar_lon + 180;
    night_center_lat = -subsolar_lat;

    % Normalize longitude to -180 to 180 range
    night_center_lon = mod(night_center_lon + 180, 360) - 180;
    
    % 2. Generate the terminator great circle (90 degrees from night center)
    % We use parametric equations for a circle on a sphere's surface.
    azimuth = linspace(0, 360, 200); % Azimuth from the night center
    [term_lat, term_lon] = reckon(night_center_lat, night_center_lon, 90, azimuth);

    % 3. Handle the Dateline Crossing
    % The raw terminator longitudes will have a large jump (e.g., from +179 to -179).
    % We must split the polygon at this jump to prevent a horizontal line artifact.
    
    % Find the wrap-around point
    lon_diff = diff(term_lon);
    wrap_idx = find(abs(lon_diff) > 180, 1);
    
    % Determine which pole is in darkness
    pole_lat = 90 * sign(-subsolar_lat);
    if pole_lat == 0, pole_lat = 90; end % Default to North if sun is on equator

    if isempty(wrap_idx)
        % No dateline crossing; the night side is fully contained in the map view.
        % This is rare but possible if the anti-solar point is at a pole.
        lon_poly = term_lon;
        lat_poly = term_lat;
        h = fill(lon_poly, lat_poly, [0 0 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    else
        % Dateline crossing found. Create two polygons.
        
        % Polygon 1: From the start of the terminator to the dateline jump
        lon1 = term_lon(1:wrap_idx);
        lat1 = term_lat(1:wrap_idx);
        edge_lon1 = lon1(end); % Longitude at the edge (+180 or -180)
        
        % Close the polygon by running along the map edge to the dark pole and back
        lon_poly1 = [lon1; edge_lon1; -edge_lon1; -edge_lon1; lon1(1)];
        lat_poly1 = [lat1; pole_lat;  pole_lat;  lat1(1);   lat1(1)];
        
        % Polygon 2: From after the dateline jump to the end
        lon2 = term_lon(wrap_idx+1:end);
        lat2 = term_lat(wrap_idx+1:end);
        edge_lon2 = lon2(1); % Longitude at the other edge
        
        % Close the second polygon similarly
        lon_poly2 = [lon2; edge_lon2; -edge_lon2; -edge_lon2; lon2(end)];
        lat_poly2 = [lat2; pole_lat;  pole_lat;  lat2(end);   lat2(end)];

        % Draw both polygons
        hold on;
        h1 = fill(lon_poly1, lat_poly1, [0 0 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        h2 = fill(lon_poly2, lat_poly2, [0 0 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        h = [h1, h2];
    end
    
    % Ensure shading is always at the bottom layer
    uistack(h, 'bottom');
    
    % Overlay terminator line for a crisp boundary (optional but nice)
    plot(term_lon, term_lat, 'w-', 'LineWidth', 0.8, 'Color', [1 1 1 0.5]);
    
    % Add sun position marker (it's better to do this here than in the main code)
    plot(subsolar_lon, subsolar_lat, 'yo', ...
         'MarkerSize', 10, 'MarkerFaceColor', 'y', ...
         'DisplayName', 'Subsolar Point');
end