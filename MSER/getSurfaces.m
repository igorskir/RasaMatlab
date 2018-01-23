function output = getSurfaces(input, isovalue)

if ~ndims(input) == 4
    warndlg('Waiting for 4D matrix as input data','Warning', 'modal');
    return
end

if islogical(input)
    input = double(input);
end

surface = struct('vertices', [], 'faces', []);

nTimeframe = size(input,4);
for i = 1:nTimeframe 
    surface(i) = isosurface(input(:,:,:,i),isovalue); 
end
output = surface;
end
