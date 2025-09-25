function OSI = compute_osi(responses, orientations)
    % Compute OSI using vector sum method
    R_pref = max(responses); % Preferred orientation response
    R_orth = mean(responses(abs(orientations - orientations(responses == R_pref) - pi/2) < 0.01)); % Orthogonal response
    OSI = (R_pref - R_orth) / (R_pref + R_orth);
end