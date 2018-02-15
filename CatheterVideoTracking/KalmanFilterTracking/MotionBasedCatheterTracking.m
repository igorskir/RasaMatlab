function MotionBasedCatheterTracking( )
% Create System objects used for reading video, detecting moving objects,
% and displaying the results.
obj = setupSystemObjects();
    
tracks = initializeTracks(); % Create an empty array of tracks.
nextId = 1; % ID of the next track

% Detect moving objects, and track them across video frames.
while ~isDone(obj.reader)
    frame = readFrame(obj);
    [centroids, bboxes, mask] = detectObjects(frame,obj);
    predictNewLocationsOfTracks();
    [assignments, unassignedTracks, unassignedDetections] = ...
        detectionToTrackAssignment();

    updateAssignedTracks();
    updateUnassignedTracks();
    deleteLostTracks();
    createNewTracks();

    displayTrackingResults();
end

end

