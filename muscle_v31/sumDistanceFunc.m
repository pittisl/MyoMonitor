% This is to compute the sum distance of a complex data array to the
% centroid
function [distanceMean,distanceVar, distance] = sumDistanceFunc(channelH, realCenter, imagCenter)
    distance = zeros(1, length(channelH));
    for i = 1:length(channelH)
        distance(i) = (real(channelH(i)) - realCenter)^2 + (imag(channelH(i))- imagCenter)^2;
    end
    distanceMean = mean(distance);
%     distanceVar = std(distance);
    distanceVar = sum(distance.^3)/length(distance);
end
