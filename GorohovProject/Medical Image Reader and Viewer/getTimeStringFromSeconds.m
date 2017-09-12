
function out = getTimeStringFromSeconds(t, format)

%inputs are time in seconds and format (DICOM or Interfile)

if ischar(t), t = str2double(t); end
if t<0
    error('Time in seconds cannot be negative.')
elseif t>=86400
    error('Time in seconds must be within a 24-hour period.')
end
if ~ischar(format)
    error('Must input time string format: ''DICOM'' or ''Interfile''.')
end

sec = mod(t,60);
min = mod((t - sec),3600) / 60;
hr = (t - min*60 - sec) / 3600;

switch lower(format(1))
    case 'd'
        out = [sprintf('%02d',hr) sprintf('%02d',min) sprintf('%09.6f',sec)];
    case 'i'
        out = [...
            sprintf('%02d',hr) ':' ...
            sprintf('%02d',min) ':' ...
            sprintf('%02d',round(sec))];
    otherwise
        error('String format must be ''DICOM'' or ''Interfile''.')
end

end


    