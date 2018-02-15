function frame = readFrame(obj)
%READFRAME Read the next video frame from the video file.
%frame = obj.reader.step();
frame = step(obj.reader);

end

