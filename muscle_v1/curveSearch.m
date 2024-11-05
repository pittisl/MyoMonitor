% This is to search the similar curve fitting results EXHAUSTIVELY
function success = curveSearch(Data_G1, Data_G2, Data_G3)
initWinLength = 5;
maxWinLength = 8;
maxPreShift_G1 = 3;
fittingOrder = 3;
% get all the curve fitting results from data group1
coefficient_G1 = zeros(fittingOrder + 1, (maxPreShift_G1 + 1)*(maxWinLength - initWinLength + 1)); % store all coefficients 1 column on fitting curve
coefficientEndWinL_G1 = zeros(2, (maxPreShift_G1 + 1)*(maxWinLength - initWinLength + 1));% store the endpoints and window lengths corresponding to coefficients in coefficient_G1 columns. 1st row, endpoint, 2nd row, window length
coefficientIndex = 1;
for endPointIndex = 0:maxPreShift_G1     % numerate from the end and shift forward
    for winLengthIndex = 1:maxWinLength - initWinLength + 1     % numerate all window lengths
        winLength = initWinLength + winLengthIndex - 1;
        [p, ~] = polyfit(1:winLength, Data_G1(end - endPointIndex - winLength + 1:end - endPointIndex), fittingOrder);
        coefficient_G1(:,coefficientIndex) = p';
        coefficientEndWinL_G1(:,coefficientIndex) = [length(Data_G1) - endPointIndex; winLength];
        coefficientIndex = coefficientIndex + 1;
    end
end

% get all the curve fitting results from data group2 and group3
coefficient_G2 = zeros(fittingOrder + 1, (length(Data_G2)-maxWinLength + 1)*(maxWinLength - initWinLength + 1)); % store all coefficients 1 column on fitting curve
coefficientEndWinL_G2 = zeros(2, (length(Data_G2)-maxWinLength + 1)*(maxWinLength - initWinLength + 1));% store the endpoints and window lengths corresponding to coefficients in coefficient_G1 columns. 1st row, endpoint, 2nd row, window length
coefficientIndex = 1;
for endPointIndex = 0:(length(Data_G2)-maxWinLength)     % numerate from the end and shift forward
    for winLengthIndex = 1:maxWinLength - initWinLength + 1     % numerate all window lengths
        winLength = initWinLength + winLengthIndex - 1;
        [p, ~] = polyfit(1:winLength, Data_G2(end - endPointIndex - winLength + 1:end - endPointIndex), fittingOrder);
        coefficient_G2(:,coefficientIndex) = p';
        coefficientEndWinL_G2(:,coefficientIndex) = [length(Data_G2) - endPointIndex; winLength];
        coefficientIndex = coefficientIndex + 1;
    end
end

coefficient_G3 = zeros(fittingOrder + 1, (length(Data_G3)-maxWinLength + 1)*(maxWinLength - initWinLength + 1)); % store all coefficients 1 column on fitting curve
coefficientEndWinL_G3 = zeros(2, (length(Data_G3)-maxWinLength + 1)*(maxWinLength - initWinLength + 1));% store the endpoints and window lengths corresponding to coefficients in coefficient_G1 columns. 1st row, endpoint, 2nd row, window length
coefficientIndex = 1;
for endPointIndex = 0:(length(Data_G3)-maxWinLength)     % numerate from the end and shift forward
    for winLengthIndex = 1:maxWinLength - initWinLength + 1     % numerate all window lengths
        winLength = initWinLength + winLengthIndex - 1;
        [p, ~] = polyfit(1:winLength, Data_G3(end - endPointIndex - winLength + 1:end - endPointIndex), fittingOrder);
        coefficient_G3(:,coefficientIndex) = p';
        coefficientEndWinL_G3(:,coefficientIndex) = [length(Data_G3) - endPointIndex; winLength];
        coefficientIndex = coefficientIndex + 1;
    end
end

% Start similarity matching with coefficients
% 1. compute similarity for all the data in group 1 and 2, total length is length(group1) * length(group2)
similarity_G1G2 = zeros(length(coefficient_G1), length(coefficient_G2));    % similarity of column x and row y
for i = 1:length(coefficient_G1)
    for j = 1:length(coefficient_G2)
        similarity_G1G2(i, j) = sum(abs(coefficient_G1(:, i) - coefficient_G2(:, j)));
    end
end
% 2. find the least 5 pairs
similaritySorted_G1G2 = sort(similarity_G1G2(:));
least5Index_G1G2 = zeros(2, 5);      % each column corresponds to a mapping to similarity_G1G2
least5_G1G2 = zeros(fittingOrder + 1, 5);
for i = 1:5
    [row_simi, col_simi] = find(similarity_G1G2 == similaritySorted_G1G2(i));
    least5Index_G1G2(:,i) = [row_simi, col_simi]';
    least5_G1G2(:,i) = mean([coefficient_G1(:,row_simi), coefficient_G2(:,col_simi)], 2); % each column is the averaged coefficients least 5 from group1 and group2
end
% 3. find the one with least difference in group 3 for each pair --- similar with step 1 and 2 combination
similarity_G12G3 = zeros(length(least5_G1G2), length(coefficient_G3));    % similarity of column x and row y
for i = 1:length(least5_G1G2)
    for j = 1:length(coefficient_G3)
        similarity_G12G3(i, j) = sum(abs(least5_G1G2(:, i) - coefficient_G3(:, j)));
    end
end
similaritySorted_G12G3 = sort(similarity_G12G3(:));
least5Index_G12G3 = zeros(2, 5);      % each column corresponds to a mapping to similarity_G1G2
for i = 1:5
    [row_simi, col_simi] = find(similarity_G12G3 == similaritySorted_G12G3(i));
    least5Index_G12G3(:,i) = [row_simi, col_simi]';
end
% 4. get the original indices in coefficient vectors of 5 most similar tuples
least5Index_G123 = zeros(3, 5);
for i = 1:5
    least5Index_G123(1:2, i) = least5Index_G1G2(:, least5Index_G12G3(1,i));
    least5Index_G123(3, i) = least5Index_G12G3(2,i);
end
% 5. derive backward getting the original window length and end point
least5IndexEndPoint = zeros(3, 5);
least5IndexWinLength = zeros(3, 5);
for i = 1:5
    least5IndexEndPoint(1, i) = coefficientEndWinL_G1(1, least5Index_G123(1, i));
    least5IndexEndPoint(2, i) = coefficientEndWinL_G2(1, least5Index_G123(2, i));
    least5IndexEndPoint(3, i) = coefficientEndWinL_G3(1, least5Index_G123(3, i));
    least5IndexWinLength(1, i) = coefficientEndWinL_G1(2, least5Index_G123(1, i));
    least5IndexWinLength(2, i) = coefficientEndWinL_G2(2, least5Index_G123(2, i));
    least5IndexWinLength(3, i) = coefficientEndWinL_G3(2, least5Index_G123(3, i));
end

% Redo the curve fitting and plot only selected window part
for i = 1:5
%     num_G1 = least5IndexEndPoint(1, i) - least5IndexWinLength(1, i) + 1:least5IndexEndPoint(1, i);
    num_G1 = 1:least5IndexWinLength(1, i);
    [p_G1, S_G1] = polyfit(num_G1, Data_G1(least5IndexEndPoint(1, i) - least5IndexWinLength(1, i) + 1:least5IndexEndPoint(1, i)), fittingOrder);
    [DataFit_G1, ~] = polyval(p_G1, num_G1, S_G1);
    [DataFit_G1_All, ~] = polyval(p_G1, -5:5, S_G1);
%     num_G2 = least5IndexEndPoint(2, i) - least5IndexWinLength(2, i) + 1:least5IndexEndPoint(2, i);
    num_G2 = 1:least5IndexWinLength(2, i);
    [p_G2, S_G2] = polyfit(num_G2, Data_G2(least5IndexEndPoint(2, i) - least5IndexWinLength(2, i) + 1:least5IndexEndPoint(2, i)), fittingOrder);
    [DataFit_G2, ~] = polyval(p_G2, num_G2, S_G2);
    [DataFit_G2_All, ~] = polyval(p_G2, -5:5, S_G2);
%     num_G3 = least5IndexEndPoint(3, i) - least5IndexWinLength(3, i) + 1:least5IndexEndPoint(3, i);
    num_G3 = 1:least5IndexWinLength(3, i);
    [p_G3, S_G3] = polyfit(num_G3, Data_G3(least5IndexEndPoint(3, i) - least5IndexWinLength(3, i) + 1:least5IndexEndPoint(3, i)), fittingOrder);
    [DataFit_G3, ~] = polyval(p_G3, num_G3, S_G3);
    [DataFit_G3_All, ~] = polyval(p_G3, -5:5, S_G3);
    figure, hold on
    plot(Data_G1, '-.b')
    plot(Data_G2, '-.r')
    plot(Data_G3, '-.m')
    plot(least5IndexEndPoint(1, i)-max(num_G1) + num_G1, DataFit_G1, 'b', 'LineWidth', 2)
    plot(least5IndexEndPoint(2, i)-max(num_G2) + num_G2, DataFit_G2, 'r', 'LineWidth', 2)
    plot(least5IndexEndPoint(3, i)-max(num_G3) + num_G3, DataFit_G3, 'm', 'LineWidth', 2)
    plot(least5IndexEndPoint(1, i)-least5IndexWinLength(1,i)-5:least5IndexEndPoint(1, i)-least5IndexWinLength(1,i) + 5, DataFit_G1_All, '--b', 'LineWidth', 2)
    plot(least5IndexEndPoint(2, i)-least5IndexWinLength(2,i)-5:least5IndexEndPoint(2, i)-least5IndexWinLength(2,i) + 5, DataFit_G2_All, '--r', 'LineWidth', 2)
    plot(least5IndexEndPoint(3, i)-least5IndexWinLength(3,i)-5:least5IndexEndPoint(3, i)-least5IndexWinLength(3,i) + 5, DataFit_G3_All, '--m', 'LineWidth', 2)
    grid on
    title('segmented 2nd-order fitting on acoustic results')
    legend('Group 1', 'Group 2', 'Group 3')
    xlabel('Workout #')
    ylabel('Fatigue Level')
    if fittingOrder == 1
        annotation('textbox', [.2 .2 .7 .7], 'String', ['a=', num2str(p_G1(1)), 10 ,'b=', num2str(p_G1(2))], 'FitBoxToText', 'on', 'Color', 'blue')
        annotation('textbox', [.2 .2 .6 .6], 'String', ['a=', num2str(p_G2(1)), 10 ,'b=', num2str(p_G2(2))], 'FitBoxToText', 'on', 'Color', 'red')
        annotation('textbox', [.2 .2 .5 .5], 'String', ['a=', num2str(p_G3(1)), 10 ,'b=', num2str(p_G3(2))], 'FitBoxToText', 'on', 'Color', 'magenta') 
    elseif fittingOrder == 2
        annotation('textbox', [.2 .2 .7 .7], 'String', ['a=', num2str(p_G1(1)), 10 ,'b=', num2str(p_G1(2)), 10 ,'c=', num2str(p_G1(3))], 'FitBoxToText', 'on', 'Color', 'blue')
        annotation('textbox', [.2 .2 .6 .6], 'String', ['a=', num2str(p_G2(1)), 10 ,'b=', num2str(p_G2(2)), 10 ,'c=', num2str(p_G2(3))], 'FitBoxToText', 'on', 'Color', 'red')
        annotation('textbox', [.2 .2 .5 .5], 'String', ['a=', num2str(p_G3(1)), 10 ,'b=', num2str(p_G3(2)), 10 ,'c=', num2str(p_G3(3))], 'FitBoxToText', 'on', 'Color', 'magenta')
    elseif fittingOrder == 3
        annotation('textbox', [.2 .2 .7 .7], 'String', ['a=', num2str(p_G1(1)), 10 ,'b=', num2str(p_G1(2)), 10 ,'c=', num2str(p_G1(3)), 10 ,'d=', num2str(p_G1(4))], 'FitBoxToText', 'on', 'Color', 'blue')
        annotation('textbox', [.2 .2 .6 .6], 'String', ['a=', num2str(p_G2(1)), 10 ,'b=', num2str(p_G2(2)), 10 ,'c=', num2str(p_G2(3)), 10 ,'d=', num2str(p_G2(4))], 'FitBoxToText', 'on', 'Color', 'red')
        annotation('textbox', [.2 .2 .5 .5], 'String', ['a=', num2str(p_G3(1)), 10 ,'b=', num2str(p_G3(2)), 10 ,'c=', num2str(p_G3(3)), 10 ,'d=', num2str(p_G3(4))], 'FitBoxToText', 'on', 'Color', 'magenta')
    end
end
end