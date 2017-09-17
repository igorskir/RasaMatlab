
function SUVconv = calculateSUV(varargin)

%patWeight is in kg

[injDose, injTime, injDate, acqTime, acqDate, patWeight, actUnits, tracer] = ...
    parseInputs(varargin{:});

hl = getHalflife(tracer);

%get time difference for decay correction
if strcmpi(acqDate,injDate)
    [injT,acqT] = getSecondsFromTimeString(injTime, acqTime);
    relT = acqT - injT; %this is to acquisition start time
else
    error('Injection and acquisition were not on same day?')
end
%calculate conversion factor
if strcmpi(actUnits,'bqml')
    corrInj = injDose / exp(log(2) * relT / hl); %in Bq and seconds
    SUVconv = patWeight / corrInj; %assuming massDensity = 1kg/L
    SUVconv = SUVconv * 1e3; %because image values mL -> L
else
    error('Account for activity units.')
end

end

function [injDose, injTime, injDate, acqTime, acqDate, patWeight, actUnits, tracer] ...
    = parseInputs(varargin)

if nargin==8
    injDose = varargin{1};
    injTime = varargin{2};
    injDate = varargin{3};
    acqTime = varargin{4};
    acqDate = varargin{5};
    patWeight = varargin{6};
    actUnits = varargin{7};
    tracer = varargin{8};
elseif nargin==1 && isstruct(varargin{1})
    if isfield(varargin{1},'injectedDose')
        injDose = varargin{1}.injectedDose;
    end
    if isfield(varargin{1},'injectionTime')
        injTime = varargin{1}.injectionTime;
    end
    if isfield(varargin{1},'injectionDate')
        injDate = varargin{1}.injectionDate;
    end
    if isfield(varargin{1},'studyTime')
        acqTime = varargin{1}.studyTime;
    end
    if isfield(varargin{1},'studyDate')
        acqDate = varargin{1}.studyDate;
    end
    if isfield(varargin{1},'patientWeight')
        patWeight = varargin{1}.patientWeight;
    end
    if isfield(varargin{1},'units')
        actUnits = varargin{1}.units;
    end
    if isfield(varargin{1},'tracer')
        tracer = varargin{1}.tracer;
    end
else
    error('Wrong inputs.')
end

if ~exist('injDose','var'), error('Injected dose is not defined.'); end
if ~exist('injTime','var'), error('Injection time is not defined.'); end
if ~exist('injDate','var'), error('Injection date is not defined.'); end
if ~exist('acqTime','var'), error('Acquisition time is not defined.'); end
if ~exist('acqDate','var'), error('Acquisition date is not defined.'); end
if ~exist('patWeight','var'), error('Patient weight is not defined.'); end
if ~exist('actUnits','var'), error('Activity units are not defined.'); end
if ~exist('tracer','var'), error('Tracer is not defined.'); end

end

function hl = getHalflife(tracer)

%returns isotope halflife (in seconds)

switch lower(tracer)
    case {'fluorodeoxyglucose','fdg'}
        hl = 6582;
    case 'thymidine(flt)'
        hl = 6582;
    case 'fluorocholine'
        hl = 6582;
    case 'choline'
        hl = 6582;
    case {'sodium','sodium fluoride f^18^'}
        hl = 6582;
    otherwise
        error(['Halflife is not defined for tracer: ' tracer]);
end

end

