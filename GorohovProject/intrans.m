function g = intrans(f, varargin)
%INTRANS Performs intensity (gray-level) transformations.
%G = INTRANS(F,'neg') for negative
%
%G = INTRANS(F,'log', C, CLASS) computes C*log(1+F). C def to 1. Class
%offers 'uint8' or 'uint16'
%
%G = INTRANS(F,'gamma', GAM) performs a gamma transfoormation
%
%G = INTRANSE(F,'stretch', M, E) computes a contrast-stretching
%transofmation using the expressin 1./(1+(M./(F + eps)).^E). M mast be in
%the range [0,1]. Def for E is 4.
%
%For the 'neg', 'gamma', and 'stretch' transformations, double 
%input images whose maximum value is greater than 1 are scaled 
%first using MAT2GRAY. Other images are converted to double first 
% using IM2D0UBLE. For the 'log' transformation, double images are 
%transformed without being scaled; other images are converted to 
% double first using IM2D0UBLE.
%

%Verify the correct number of unputs
narginchk(2, 4)

%Store the class of the input for use later
classin = class(f);
% If the input is of class double, and it is outside the range 
% [0, 1], and the specified transformation is not 'log', convert the 
% input to the range [0, 1] 
if strcmp(class(f), 'double') && max(f(:)) > 1 && ~strcmp(varargin{1}, 'log')
    f = mat2gray(f);
else %Converte to double, regardless of class(f)
    f = im2double(f);
end
%Determine the type of transformation specified.
method = varargin{1};
%Performe the intensity transforrmation specified
switch method
    case 'neg'
        g = imcomplement(f);
    case 'log'
        if length(varargin) == 1
            c = 1;
        elseif length(varargin) == 2
            c = varargin{2};
        elseif length(varargin) == 3
            c = varargin{2};
            classin = varargin{3};
        else
            error('Incorrect number of inputs for the log option.')
        end
        g = c*(log(1 + double(f)));
    case 'gamma'
        if length(varargin) < 2
            error('Not enough inputs for the gamma option.'); 
        end
        gam = varargin{2};
        g = imadjust(f,[],[],gam);
    case 'stretch'
        if length(varargin) == 1
            % Use defaults.
            m = mean2(f);
            E = 4.0;
        elseif length(varargin) == 3
            m = varargin{2};
            E = varargin{3};
        else
            error('Incorrect number of inputs for the stretch option.')
        end
        g = 1./(1 + (m./(f + eps)).^E);
    otherwise
        error('Unknown enhancement method.');
end    
% Convert to the class of the input image
%g = changeclass(classin,g);
            
end