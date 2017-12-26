function [ dir ] = TurnRight( direction )
%this functions returns  vector rotated 90 degrees left
% input [-1/0/1 -1/0/1] format
% one of the element have to be 0 (directions along x or Y axis possible)
dir=[direction(2), -direction(1)];
end

