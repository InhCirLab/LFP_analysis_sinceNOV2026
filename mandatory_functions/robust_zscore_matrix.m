function z = robust_zscore_matrix(data)
    % data: [frequencies x timepoints]
    z = zeros(size(data));
    scale = 1.4826; % scaling constant for MAD

    for i = 1:size(data, 1)  % loop over each frequency row
        row = data(i, :);
        med = median(row);
        mad_val = median(abs(row - med));

        if mad_val == 0
            z(i, :) = zeros(1, size(data, 2));
        else
            z(i, :) = (row - med) / (mad_val * scale);
        end
    end
end
