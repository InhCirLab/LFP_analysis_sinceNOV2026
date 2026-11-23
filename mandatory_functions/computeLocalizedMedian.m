function locS_det_median = computeLocalizedMedian(S_det, intervals)
    if isempty(intervals)
        locS_det_median = nan(1, size(S_det, 2));
    else
        locS = [];
        for i = 1:size(intervals,1)
            locS = [locS; S_det(intervals(i,1):intervals(i,2), :)];
        end
        locS_det_median = median(locS, 1, 'omitnan');
    end
end