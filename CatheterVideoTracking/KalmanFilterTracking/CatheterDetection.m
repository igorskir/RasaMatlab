obj = setupSystemObjects();
    
tracks = initializeTracks(); % Create an empty array of tracks.
nextId = 1; % ID of the next track

% Detect moving objects, and track them across video frames.
while ~isDone(obj.reader)
    frame = readFrame(obj);
    [centroids, bboxes, mask] = detectObjects(frame,obj);
    predictNewLocationsOfTracks(tracks);
    [assignments, unassignedTracks, unassignedDetections] = ...
        detectionToTrackAssignment(tracks,centroids);

    updateAssignedTracks(assignments);
%     updateUnassignedTracks();
%     deleteLostTracks();
%     createNewTracks();

    displayTrackingResults(obj,frame,mask,tracks);
end