% This is to remove the outliers in the fatigue data (BASED ON WINDOWED MAD)
function fatigueValAfterORM = ORM_windowed(fatigueVal)
    winLength = round(length(fatigueVal)/3);
    initMADVal = 1;
    numOutlierRatio = 0.25;
    while 1
        outlierIndex = isoutlier(fatigueVal, 'movmedian', winLength, 'ThresholdFactor',initMADVal);
        if sum(outlierIndex) <= length(fatigueVal) * numOutlierRatio
            break;
        else
            initMADVal = initMADVal + 0.1;
        end
    end
    fatigueValAfterORM = fatigueVal;
    fatigueValAfterORM(outlierIndex) = [];
end