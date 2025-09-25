function CV = compute_cv(responses, orientations)
    % Compute Circular Variance using vector sum method
    vector_sum = sum(responses .* exp(1i * 2 * orientations));
    CV = 1 - abs(vector_sum) / sum(responses);
end
