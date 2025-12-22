function h_patch = patchCleanNightShading(lonTerm, latTerm)
    % Creates a clean night shading polygon without accumulation issues
    % Similar to what worked in static mode
    
    % Ensure vectors are column vectors
    lonTerm = lonTerm(:);
    latTerm = latTerm(:);
    
    % Find the westernmost point (minimum longitude) on terminator
    [~, west_idx] = min(lonTerm);
    
    % Reorder to start from westernmost point
    lonTerm = circshift(lonTerm, -west_idx+1);
    latTerm = circshift(latTerm, -west_idx+1);
    
    % Create a closed polygon for the night side
    % Strategy: Create polygon from western edge to terminator to eastern edge
    
    % Western edge to terminator start
    lon_poly = [-180; lonTerm; 180];
    lat_poly = [90; latTerm; latTerm(end)];
    
    % Complete the polygon (eastern edge back to western edge)
    lon_poly = [lon_poly; 180; -180; -180];
    lat_poly = [lat_poly; -90; -90; 90];
    
    % Create the night patch
    h_patch = patch(lon_poly, lat_poly, [0 0 0], ...  % Black color
                    'EdgeColor', 'none', ...
                    'FaceAlpha', 0.6, ...            % 60% opacity
                    'FaceColor', [0 0 0]);           % Pure black
    
    % Move to bottom layer so other graphics stay on top
    uistack(h_patch, 'bottom');
end