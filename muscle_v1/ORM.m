% This is to remove the outliers in the fatigue data (BASED ON CONTINUOUS VALUE)
function fatigueValAfterORM = ORM(fatigueVal)
    factor_ratio = 3;
    while 1
        fatigueValAfterORM_1 = fatigueVal(1);
        indexAfterORM_1 = 1;
        diffVal = mad(fatigueVal);
        for i = 1 : length(fatigueVal) - 1
            if abs(fatigueVal(i+1) - fatigueValAfterORM_1(indexAfterORM_1)) <= abs(diffVal) * factor_ratio
                fatigueValAfterORM_1 = [fatigueValAfterORM_1, fatigueVal(i+1)];
                indexAfterORM_1 = indexAfterORM_1 + 1;
            end
        end
        % if more than 20% of the data has been removed, increase the ratio
        if length(fatigueVal) - length(fatigueValAfterORM_1) > length(fatigueVal) * 0.2
            factor_ratio = factor_ratio + 0.1;
            fatigueValAfterORM_1 = [];
        else
            break
        end
    end
    % figure,plot(repreValueFatigue_wo_out)
    factor_ratio = factor_ratio - 0.5;
    while 1
        fatigueValAfterORM = fatigueValAfterORM_1(end);
        diffVal = mad(fatigueValAfterORM_1);
        for i = indexAfterORM_1 : -1: 2
            if abs(fatigueValAfterORM_1(i-1) - fatigueValAfterORM(1)) <= abs(diffVal) * factor_ratio
                fatigueValAfterORM = [fatigueValAfterORM_1(i-1), fatigueValAfterORM];
            end
        end
        % if more than 20% of the data has been removed, increase the ratio
        if length(fatigueValAfterORM_1) - length(fatigueValAfterORM) > length(fatigueValAfterORM_1) * 0.2
            factor_ratio = factor_ratio + 0.1;
            fatigueValAfterORM = [];
        else
            break
        end
    end
end



%% Backup
% repreValueFatigue_wo_out = repreValueFatigue(1);
% index_wo_out = 1;
% diffVal = mad(repreValueFatigue);
% factor_ratio = 3;
% for i = 1 : length(repreValueFatigue) - 1
%     if abs(repreValueFatigue(i+1) - repreValueFatigue_wo_out(index_wo_out)) <= abs(diffVal) * factor_ratio
%         repreValueFatigue_wo_out = [repreValueFatigue_wo_out, repreValueFatigue(i+1)];
%         index_wo_out = index_wo_out + 1;
%     end
% end
% % figure,plot(repreValueFatigue_wo_out)
% repreValueFatigue_wo_out_rev = repreValueFatigue_wo_out(end);
% diffVal = mad(repreValueFatigue_wo_out);
% factor_ratio = 2.5;
% for i = index_wo_out : -1: 2
%     if abs(repreValueFatigue_wo_out(i-1) - repreValueFatigue_wo_out_rev(1)) <= abs(diffVal) * factor_ratio
%         repreValueFatigue_wo_out_rev = [repreValueFatigue_wo_out(i-1), repreValueFatigue_wo_out_rev];
%     end
% end
% % figure,plot(repreValueFatigue_wo_out_rev)
% % % use the function provided by Matlab
% % outlier_indices = isoutlier(repreValueFatigue, 'ThresholdFactor', 1.5);
% % repreValueFatigue_wo_out = repreValueFatigue;
% % repreValueFatigue_wo_out(outlier_indices) = [];