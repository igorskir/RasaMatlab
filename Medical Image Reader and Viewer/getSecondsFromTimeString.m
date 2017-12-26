
function varargout = getSecondsFromTimeString(varargin)

%if 2 inputs, order by assumed (maybe not actual) time, e.g. injTime, acqTime
%and it assumes that 2 times were from the same day

tString1 = varargin{1};
if ~ischar(tString1), tString1 = char(tString1); end

temp = regexp(tString1, ':', 'split');
if length(temp)==2
    temp = [temp {'00'}];
end
if numel(temp)>1 %Interfile
    hr1 = str2double(temp{1});
    min1 = str2double(temp{2});
    sec1 = str2double(temp{3});
else
    temp = regexp(tString1, '\.', 'split');
    if length(temp)==2
        temp = [temp {'00'}];
    end
    if length(temp{1})==6 %DICOM
        temp = round(str2double(tString1));
        hr1 = floor(temp / 1e4);
        min1 = floor((temp - hr1 * 1e4) / 1e2);
        sec1 = mod(temp, 1e2);
    else
        error('Account for time format.')
    end
end

if nargin==2    
    tString2 = varargin{2};
    if ~ischar(tString2), tString2 = char(tString2); end
    
    temp = regexp(tString2, ':', 'split');
    if length(temp)==2
        temp = [temp {'00'}];
    end
    if numel(temp)>1 %Interfile
        hr2 = str2double(temp{1});
        min2 = str2double(temp{2});
        sec2 = str2double(temp{3});
    else
        temp = regexp(tString2, '\.', 'split');
        if length(temp)==2
            temp = [temp {'00'}];
        end
        if length(temp{1})==6 %DICOM
            temp = round(str2double(tString2));
            hr2 = floor(temp / 1e4);
            min2 = floor((temp - hr2 * 1e4) / 1e2);
            sec2 = mod(temp, 1e2);
        else
            error('Account for time format.')
        end
    end
    
    if hr2<hr1, hr2 = hr2 + 12; end

end

varargout{1} = 3600*hr1 + 60*min1 + sec1; %in seconds
if nargin==2
    varargout{2} = 3600*hr2 + 60*min2 + sec2; %in seconds
    if varargout{2}<varargout{1}
        warning('Time 2 is before time 1');
    end
end

end


    