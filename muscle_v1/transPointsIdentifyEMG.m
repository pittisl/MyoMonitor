% This is to identify the transition point for EMG signal
function transPoints = transPointsIdentifyEMG(rawData)
    transPoints = zeros(4,1);
    dura_stat = 3; dura_power = 0; dura_lift = 2; dura_hold = 3; dura_put = 2; dura_wait = 1; % standard test timing
    dura_total = dura_stat + dura_power + dura_lift + dura_hold + dura_put + dura_wait;
    % cut the begining steady period
    rawData = rawData(ceil(length(rawData) * 3/11):end);
    % compute differences
    diff_step = 1;
    rawData_diff = zeros(1, length(rawData) - 1);
    for i = 1:diff_step:length(rawData) - diff_step
        rawData_diff(i) = rawData(i + diff_step) - rawData(i);
    end
    % 1st transition point by difference
    p1_step = 10;
    for i = 1:length(rawData_diff)-p1_step
        if sum(abs(rawData_diff(i:i+p1_step))) >= 4 
            transPoints(1) = i + ceil(length(rawData) * 3/8) - 1; % have to add the removed part
            break
        end
    end
    % 2nd transition point
    seg_1 = rawData(ceil(length(rawData) * 2/8):floor(length(rawData) * 3/8));
    seg_2 = rawData(ceil(length(rawData) * 2.5/8):floor(length(rawData) * 3.5/8));
    seg_3 = rawData(ceil(length(rawData) * 3/8):floor(length(rawData) * 4/8));
    seg_4 = rawData(ceil(length(rawData) * 3.5/8):floor(length(rawData) * 4.5/8));
    seg_5 = rawData(ceil(length(rawData) * 4/8):floor(length(rawData) * 5/8));
    seg_6 = rawData(ceil(length(rawData) * 4.5/8):floor(length(rawData) * 5.5/8));
    maxmin_vector = [max(seg_1)-min(seg_1), max(seg_2)-min(seg_2), ...
                        max(seg_3)-min(seg_3), max(seg_4)-min(seg_4), ...
                            max(seg_5)-min(seg_5), max(seg_6)-min(seg_6)];
    p2_step = round(length(rawData) * 0.5/8);
    for i = ceil(length(rawData) * 2/8):floor(length(rawData) * 5.5/8)-p2_step
        if max(rawData(i:i+p2_step)) - min(rawData(i:i+p2_step)) <= 1.1 * median(maxmin_vector)
            transPoints(2) = i + ceil(length(rawData) * 3/8) - 1;
            break
        end
    end
    % 3rd transition point
    p3_step = p2_step;
    for i = ceil(length(rawData) * 5/8):floor(length(rawData) * 7/8)-p3_step
        if max(rawData(i:i+p3_step)) - min(rawData(i:i+p3_step)) >= 1.2 * median(maxmin_vector)
            transPoints(3) = i + ceil(length(rawData) * 3/8) - 1;
            break
        end
    end
    % 4th transition point
    p4_step = 20;
    for i = length(rawData_diff): -1: floor(length(rawData) * 6.5/8)
        if sum(abs(rawData_diff(i-p4_step:i))) >= 4 
            transPoints(4) = i + ceil(length(rawData) * 3/8) - 1; % have to add the removed part
            break
        end
        if i == floor(length(rawData) * 6.5/8)
            transPoints(4) = length(rawData) + floor(length(rawData) * 3/8) - 1; 
        end
    end
%     figure,subplot(211), plot(rawData), 
%     subplot(212),plot(rawData_diff)
end