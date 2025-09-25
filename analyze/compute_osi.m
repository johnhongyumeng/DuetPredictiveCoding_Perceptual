function OSI = compute_osi(responses, orientations)
    % Compute OSI using vector sum method
    R_pref = max(responses); % Preferred orientation response
    a= abs(orientations - orientations(responses == R_pref) - pi/2) < 0.01;
    b= abs(orientations - orientations(responses == R_pref) + pi/2) < 0.01;
    if sum(a)>0
        R_orth = responses(a);
    else
        R_orth = responses(b);
    end
    % OSI = (R_pref - R_orth) / (abs(R_pref) + abs(R_orth));
    OSI = R_pref - R_orth;
end