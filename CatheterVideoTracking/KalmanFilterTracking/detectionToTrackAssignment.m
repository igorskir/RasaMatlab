function [assignments, unassignedTracks, unassignedDetections] = detectionToTrackAssignment(tracks,centroids)
%DETECTIONTOTRACKASSIGNMENT function uses the Munkres' version of the Hungarian algorithm to compute an assignment which minimizes the total cost

    nTracks = length(tracks);
    nDetections = size(centroids, 1);

    % Compute the cost of assigning each detection to each track.
    cost = zeros(nTracks, nDetections);
    for i = 1:nTracks
        cost(i, :) = distance(tracks(i).kalmanFilter, centroids);
    end

    % Solve the assignment problem.
    costOfNonAssignment = 20;
    [assignments, unassignedTracks, unassignedDetections] = ...
        assignDetectionsToTracks(cost, costOfNonAssignment);

end

