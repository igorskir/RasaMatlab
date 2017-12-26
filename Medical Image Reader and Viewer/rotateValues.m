function [newVal1, newVal2, newVal3] = rotateValues(varargin)

%input should be like that of 'readImages' structure, NOT dicom

[val1,val2,val3,oldIOP,newIOP] = parseImputs(varargin{:});

%old rotation (aligns to scanner)
xr = oldIOP(1);
yr = oldIOP(2);
zr = oldIOP(3);
xc = oldIOP(4);
yc = oldIOP(5);
zc = oldIOP(6);
xs = yr*zc-yc*zr;
ys = zr*xc-zc*xr;
zs = xr*yc-xc*yr;
oldRot = [xr xc xs; yr yc ys; zr zc zs];
%new rotation (transpose aligns to newIOP)
xr = newIOP(1);
yr = newIOP(2);
zr = newIOP(3);
xc = newIOP(4);
yc = newIOP(5);
zc = newIOP(6);
xs = yr*zc-yc*zr;
ys = zr*xc-zc*xr;
zs = xr*yc-xc*yr;
newRot = [xr xc xs; yr yc ys; zr zc zs];
%total rotation is product of both (basically rotational difference between old and new)
rot = oldRot'*newRot;
temp = sum(abs(rot'*[val1 0 0; 0 val2 0; 0 0 val3]),2);
newVal1 = temp(1);
newVal2 = temp(2);
newVal3 = temp(3);

end


function [val1,val2,val3,oldIOP,newIOP] = parseImputs(varargin)

if nargin==4
    val1 = double(varargin{1});
    val2 = double(varargin{2});
    val3 = double(varargin{3});
    oldIOP = double(varargin{4});
    newIOP = [1;0;0;0;1;0]; %if not input, assume to align to scanner
elseif nargin==5
    val1 = double(varargin{1});
    val2 = double(varargin{2});
    val3 = double(varargin{3});
    oldIOP = double(varargin{4});
    newIOP = double(varargin{5}(:));
% elseif isstruct(varargin{1})
%     if nargin==1
%         %if structure input (from 'readImages') use image FOVs
%         oldIOP = double(varargin{1}.IOP);
%         [a, b, c] = getImageBedRanges(varargin{1});
%         val1 = diff(a);
%         val2 = diff(b);
%         val3 = diff(c);
%         newIOP = [1;0;0;0;1;0]; %if not input, assume to align to scanner
%     elseif nargin==2
%         oldIOP = varargin{1}.IOP;
%         [a, b, c] = getImageBedRanges(varargin{1});
%         val1 = diff(a);
%         val2 = diff(b);
%         val3 = diff(c);
%         newIOP = varargin{2}(:);
%     end
else
    error('Incorrect inputs.')
end

%check new IOP
newIOP = newIOP(:);
if ~all(size(newIOP)==[6 1])
    error('Invalid patient orientation input vector size');
else
    %assure normalized vectors
    newIOP(1:3) = newIOP(1:3)/norm(newIOP(1:3));
    newIOP(4:6) = newIOP(4:6)/norm(newIOP(4:6));
    %check if valid IOP vector
    if any(isnan(newIOP)) || (round(newIOP(1:3)'*newIOP(4:6)*1e6)/1e6)~=0
        error('Invalid patient orientation input vector');
    end
end

end
