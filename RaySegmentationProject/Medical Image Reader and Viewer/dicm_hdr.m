function [s, info, dict] = dicm_hdr(fname, dict, iFrames)

% Return header of a dicom file in a struct.
%
% [s, err] = DICM_HDR(dicomFileName, dict, iFrames);
%
% The mandatory 1st input is the dicom file name.
%
% The optional 2nd input is a dicom dictionary returned by dicm_dict. It may
% have only part of the full dictionary, which can speed up header read
% considerably. See rename_dicm for example.
%
% The optional 3rd intput is useful for multi-frame dicom. When there are many
% frames, it may be slow to read all frames in PerFrameFunctionalGroupsSequence.
% The 3rd input can be used to specify the frames to read. By default, items for
% only 1st, 2nd and last frames are read.
%
% The optional 2nd output contains information in case of error, and will be
% empty if there is no error.
%
% DICM_HDR is like Matlab dicominfo, but is independent of Image Processing
% Toolbox. The advantage is that it decodes most private and shadow tags for
% Siemens, GE and Philips dicom, and runs faster, especially for partial header
% and multi-frame dicom.
%
% DICM_HDR can also read Philips PAR file, AFNI HEAD file and some BrainVoyager
% files, and return needed fields for dicm2nii to convert into nifti.
%
% This function was adopted from
%  'DICOM to NIfTI converter, NIfTI tool and viewer' written by Xiangrui Li

persistent dict_full;
s = []; info = '';
p.fullHdr = false; % p for parameters: only updated in main func
if nargin<2 || isempty(dict)
    if isempty(dict_full), dict_full = dicm_dict; end
    p.fullHdr = true;
    p.dict = dict_full;
    dict = [];
else
    p.dict = dict;
end

if nargin==3 && isstruct(fname) % wrapper
    s = search_MF_val(fname, dict, iFrames); % s, s1, iFrames
    return;
end

if nargin<3, iFrames = []; end
p.iFrames = iFrames;
fid = fopen(fname);
if fid<0
    info = ['File not exists: ' fname];
    return;
%     error(['File not exists: ' fname]);
end
closeFile = onCleanup(@() fclose(fid)); % auto close when done or error
fseek(fid, 0, 1); fSize = ftell(fid); fseek(fid, 0, -1);
b8 = fread(fid, 130000, 'uint8=>uint8')'; % enough for most dicom
if fSize<140 % 132 + one empty tag, ignore truncated
    info = ['Invalid file: ' fname];
    return;
%     error(['Invalid file: ' fname]);
end

iTagStart = 132; % start of first tag
isDicm = isequal(b8(129:132), 'DICM');
if ~isDicm % truncated dicom: is first group 2 or 8? not safe
    group = ch2int16(b8(1:2), 0);
    isDicm = group==2 || group==8; % truncated dicm always LE
    iTagStart = 0;
end
if ~isDicm % may be PAR/HEAD/BV file
    [~, ~, ext] = fileparts(fname);
    try
        if strcmpi(ext, '.PAR') % || strcmpi(ext, '.REC')
            [s, info] = philips_par(fname);
        elseif strcmpi(ext, '.HEAD') % || strcmpi(ext, '.BRIK')
            [s, info] = afni_head(fname);
        elseif any(strcmpi(ext, {'.vmr' '.fmr' '.dmr'})) % BrainVoyager
            [s, info] = bv_file(fname);
        else
            info = ['Unknown file type: ' fname];
            return;
%             error(['Unknown file type: ' fname]);
        end
    catch me
        info = me.message;
        return;
%         error(me.message);
    end
    if ~isempty(s), return; end
    if isempty(info)
        info = ['Not dicom file: ' fname];
%         error(['Not dicom file: ' fname]);
    end
    return;
end

p.expl = false; % default for truncated dicom
p.be = false; % little endian by default

% Get TransferSyntaxUID first, so find PixelData
i = strfind(char(b8), [char([2 0 16 0]) 'UI']); % always explicit LE
tsUID = '';
if ~isempty(i) % empty for truncated
    i = i(1) + 6;
    n = ch2int16(b8(i+(0:1)), 0);
    tsUID = dcm_str(b8(i+1+(1:n)));
    p.expl = ~strcmp(tsUID, '1.2.840.10008.1.2'); % may be wrong for some
    p.be =    strcmp(tsUID, '1.2.840.10008.1.2.2');
end

% find s.PixelData.Start so avoid search in img
% We verify iPixelData+bytes=fSize. If bytes=2^32-1, read all and use last tag
tg = char([224 127 16 0]); % PixelData, VR can be OW/OB even if expl
if p.be, tg = tg([2 1 4 3]); end
found = false;
for nb = [0 2e6 20e6 fSize] % if not enough, read more till all read
    b8 = [b8 fread(fid, nb, 'uint8=>uint8')']; %#ok
    i = strfind(char(b8), tg);
    i = i(mod(i,2)==1); % must be odd number
    for k = i(end:-1:1) % last is likely real PixelData
        p.iPixelData = k + p.expl*4 + 7; % s.PixelData.Start: 0-based
        if numel(b8)<p.iPixelData, b8 = [b8 fread(fid, 8, '*uint8')']; end %#ok
        p.bytes = ch2int32(b8(p.iPixelData+(-3:0)), p.be);
        if p.bytes==4294967295 && feof(fid), break; end % 2^32-1 compressed
        d = fSize - p.iPixelData - p.bytes; % d=0 most of time
        if d>=0 && d<16, found = true; break; end % real PixelData
    end
    if found, break; end
    if feof(fid)
        if isempty(i)
            info = ['No PixelData in ' fname];
            return;
%             error(['No PixelData in ' fname]);
        else break;
        end
    end
end

s.Filename = fopen(fid);
s.FileSize = fSize;

nTag = numel(p.dict.tag); % always search if only one tag: can find it in any SQ
toSearch = nTag<2 || (nTag<30 && ~any(strcmp(p.dict.vr, 'SQ')) && p.iPixelData<1e6);
if toSearch % search each tag if header is short and not many tags asked
    if ~isempty(tsUID), s.TransferSyntaxUID = tsUID; end % hope it is 1st tag
    ib = p.iPixelData; % will be updated each loop
    if ~isempty(p.dict.vendor)
        tg = char([8 0 112 0]); % Manufacturer
        if p.be, tg = tg([2 1 4 3]); end
        if p.expl, tg = [tg 'LO']; end
        i = strfind(char(b8(1:ib)), tg);
        i = i(mod(i,2)==1);
        if ~isempty(i)
            i = i(1) + 4 + p.expl*2; % Manufacturer should be the earlier one
            n = ch2int16(b8(i+(0:1)), p.be);
            dat = dcm_str(b8(i+1+(1:n)));
            [p, dict] = updateVendor(p, dat);
        end
    end
    
    s1 = [];
    for k = numel(p.dict.tag):-1:1 % reverse order so reduce range each loop
        group = p.dict.group(k);
        swap = p.be && group~=2;
        hasVR = p.expl || group==2;
        tg = char(typecast([group p.dict.element(k)], 'uint8'));
        if swap, tg = tg([2 1 4 3]); end
        i = strfind(char(b8(1:ib)), tg);
        i = i(mod(i,2)==1);
        if isempty(i), continue; % no this tag, next
        elseif numel(i)>1 % +1 tags found, add vr to reduce if expl
            if hasVR
                tg = [tg p.dict.vr{k}]; %#ok
                i = strfind(char(b8(1:ib)), tg);
                i = i(mod(i,2)==1);
                if numel(i)~=1, toSearch = false; break; end
            else
                toSearch = false; break; % switch to regular way
            end
        end
        ib = i-1; % make next tag search faster
        i = i + 4; % tag
        if hasVR
            vr = char(b8(i+(0:1))); i = i+2;
            if strcmp(vr, 'UN') || strcmp(vr, 'OB'), vr = p.dict.vr{k}; end
        else
            vr = p.dict.vr{k};
        end
        [n, nvr] = val_len(vr, b8(i+(0:5)), hasVR, swap); i = i+nvr;
        if n==0, continue; end % dont assign empty tag
        
        [dat, info] = read_val(b8(i+(0:n-1)), vr, swap);
        if ~isempty(info), toSearch = false; break; end % re-do in regular way
        if ~isempty(dat), s1.(p.dict.name{k}) = dat; end
    end
    if ~isempty(s1) % reverse the order: just make it look normal
        nam = fieldnames(s1);
        for k = numel(nam):-1:1, s.(nam{k}) = s1.(nam{k}); end
    end
end

i = iTagStart + 1;
while ~toSearch
    if i >= p.iPixelData
        if strcmp(name, 'PixelData') % iPixelData might be in img
            p.iPixelData = iPre + p.expl*4 + 7; % overwrite previous
            p.bytes = ch2int32(b8(p.iPixelData+(-3:0)), p.be);
        else
            info = ['End of file reached: likely error: ' s.Filename];
%             error(['End of file reached: likely error: ' s.Filename]);
        end
        break; % done or give up
    end
    iPre = i; % back it up for PixelData
    [dat, name, info, i, tg] = read_item(b8, i, p);
    if ~isempty(info), break; end
    if ~p.fullHdr && tg>p.dict.tag(end), break; end % done for partial hdr
    if isempty(dat) || isnumeric(name), continue; end
    s.(name) = dat;
    if strcmp(name, 'Manufacturer')
        [p, dict] = updateVendor(p, dat);
    elseif tg>=2621697 && ~isfield(p, 'nFrames') % BitsAllocated
        p = get_nFrames(s, p); % only make code here cleaner
    end
end

s.PixelData.Start = p.iPixelData;
s.PixelData.Bytes = p.bytes;

if isfield(s, 'CSAImageHeaderInfo') % Siemens CSA image header
    s.CSAImageHeaderInfo = read_csa(s.CSAImageHeaderInfo);
end
if isfield(s, 'CSASeriesHeaderInfo') % series header
    s.CSASeriesHeaderInfo = read_csa(s.CSASeriesHeaderInfo);
end
if isfield(s, 'ProtocolDataBlock') % GE
    s.ProtocolDataBlock = read_ProtocolDataBlock(s.ProtocolDataBlock);
end

%% nested function: update Manufacturer
    function [p, dict] = updateVendor(p, vendor)
        if ~isempty(p.dict.vendor) && strncmpi(vendor, p.dict.vendor, 2)
            dict = p.dict; % in case dicm_hdr asks 3rd output
            return;
        end
        dict_full = dicm_dict(vendor);
        if ~p.fullHdr && isfield(p.dict, 'fields')
            dict = dicm_dict(vendor, p.dict.fields);
        else
            dict = dict_full;
        end
        p.dict = dict;
    end

end % main func

%% subfunction: read dicom item. Called by dicm_hdr and read_sq
function [dat, name, info, i, tag] = read_item(b8, i, p)
dat = []; name = nan; info = ''; vr = 'CS'; % vr may be used if implicit

group = b8(i+(0:1)); i=i+2;
swap = p.be && ~isequal(group, [2 0]); % good news: no 0x0200 group
group = ch2int16(group, swap);
elmnt = ch2int16(b8(i+(0:1)), swap); i=i+2;
tag = group*65536 + elmnt;
if tag == 4294893581 % || tag == 4294893789 % FFFEE00D or FFFEE0DD
    i = i+4; % skip length
    return; % rerurn in case n is not 0
end

swap = p.be && group~=2;
hasVR = p.expl || group==2;
if hasVR, vr = char(b8(i+(0:1))); i = i+2; end % 2-byte VR

[n, nvr] = val_len(vr, b8(i+(0:5)), hasVR, swap); i = i + nvr;
if n==0, return; end % empty val

% Look up item name in dictionary
ind = find(p.dict.tag == tag, 1);
if ~isempty(ind)
    name = p.dict.name{ind};
    if strcmp(vr, 'UN') || strcmp(vr, 'OB') || ~hasVR, vr = p.dict.vr{ind}; end
elseif tag==524400 % in case not in dict
    name = 'Manufacturer';
elseif tag==131088 % need TransferSyntaxUID even if not in dict
    name = 'TransferSyntaxUID';
elseif tag==593936 % 0x0009 0010 GEIIS not dicom compliant
    i = i+n; return; % seems n is not 0xffffffff
elseif p.fullHdr
    if elmnt==0, i = i+n; return; end % skip GroupLength
    if mod(group,2), name = sprintf('Private_%04x_%04x', group, elmnt);
    else
        name = sprintf('Unknown_%04x_%04x', group, elmnt);
        % lookup unknown names
        newName = dicom_lookup(name);
        newName = '';
        if ~isempty(newName)
            name = newName;
        end
    end
    if ~hasVR, vr = 'UN'; end % not in dict, will leave as uint8
elseif n<4294967295 % no skip for SQ with length 0xffffffff
    i = i+n; return;
end

% compressed PixelData, n can be 0xffffffff
if ~hasVR && n==4294967295, vr = 'SQ'; end % best guess
if n+i>p.iPixelData && ~strcmp(vr, 'SQ'), i = i+n; return; end % PixelData or err
% fprintf('(%04x %04x) %s %6.0f %s\n', group, elmnt, vr, n, name);

if strcmp(vr, 'SQ')
    nEnd = min(i+n, p.iPixelData); % n is likely 0xffff ffff
    [dat, info, i] = read_sq(b8, i, nEnd, p, tag==1375769136); % isPerFrameSQ?
else
    [dat, info] = read_val(b8(i+(0:n-1)), vr, swap); i=i+n;
end

end

%% Subfunction: decode SQ, called by read_item (recursively)
% SQ structure:
%  while isItem (FFFE E000, Item) % Item_1, Item_2, ...
%   loop tags under the Item till FFFE E00D, ItemDelimitationItem
%   return if FFFE E0DD SequenceDelimitationItem (not checked)
function [rst, info, i] = read_sq(b8, i, nEnd, p, isPerFrameSQ)
rst = []; info = ''; tag1 = []; j = 0; % j is SQ Item index

while i<nEnd % loop through multi Item under the SQ
    tag = b8(i+([2 3 0 1])); i = i+4;
    if p.be, tag = tag([2 1 3 4]); end
    tag = ch2int32(tag, 0);
    if tag ~= 4294893568, i = i+4; return; end % only do FFFE E000, Item
    n = ch2int32(b8(i+(0:3)), p.be); i = i+4; % n may be 0xffff ffff
    n = min(i+n, nEnd);
    j = j + 1;
    
    % This 'if' block deals with partial PerFrameSQ: j and i jump to a frame.
    % The 1/2/nf frame scheme will have problem in case that tag1 in 1st frame
    % is not the first tag in other frames. Then the tags before tag1 in other
    % frames will be treated as for previous.
    if isPerFrameSQ
        if ischar(p.iFrames) % 'all' frames
            if j==1 && ~isnan(p.nFrames), rst.FrameStart = nan(1, p.nFrames); end
            rst.FrameStart(j) = i-9;
        elseif j==1, i0 = i-8; % always read 1st frame, save i0 in case of re-do
        elseif j==2 % always read 2nd frame, and find start ind for all frames
            if isnan(p.nFrames) || isempty(tag1) % 1st frame has no asked tag
                p.iFrames = 'all'; rst = []; j = 0; i = i0; % re-do the SQ
                continue; % re-do
            end
            tag1 = char(typecast(uint32(tag1), 'uint8'));
            tag1 = tag1([3 4 1 2]);
            if p.be && ~isequal(tag1(1:2),[2 0]), tag1 = tag1([2 1 4 3]); end
            ind = strfind(char(b8(i0:p.iPixelData)), tag1) + i0 - 1;
            ind = ind(mod(ind,2)==1);
            nInd = numel(ind);
            if nInd ~= p.nFrames
                tag1PerF = nInd / p.nFrames;
                if mod(tag1PerF, 1) > 0 % not integer, read all frames
                    p.iFrames = 'all'; rst = []; j = 0; i = i0; % re-do SQ
                    fprintf(2, ['Failed to determine indice for frames. ' ...
                        'Reading all frames. Maybe slow ...\n']);
                    continue;
                elseif tag1PerF>1 % more than one ind for each frame
                    ind = ind(1:tag1PerF:nInd);
                    nInd = p.nFrames;
                end
            end
            rst.FrameStart = ind-9; % 0-based
            
            iItem = 2; % initialize here. Increase when j>2
            iFrame = unique([1 2 round(p.iFrames) nInd]);
        else % overwrite j with asked iFrame, overwrite i with start ind
            iItem = iItem + 1;
            j = iFrame(iItem);
            i = ind(j); % start of tag1 for asked frame
            n = nEnd; % use very end of the sq
        end
    end
    
    Item_n = sprintf('Item_%g', j);
    while i<n % loop through multi tags under one Item
        [dat, name, info, i, tag] = read_item(b8, i, p);
        if tag == 4294893581, break; end % FFFE E00D ItemDelimitationItem
        if isempty(tag1), tag1 = tag; end % first detected tag for PerFrameSQ
        if isempty(dat) || isnumeric(name), continue; end % 0-length or skipped
        rst.(Item_n).(name) = dat;
    end
end
end

%% subfunction: cast uint8/char to double. Better performance than typecast
function d = ch2int32(u8, swap)
if swap, u8 = u8([4 3 2 1]); end
d = double(u8);
d = d(1) + d(2)*256 + d(3)*65536 + d(4)*16777216; % d = d * 256.^(0:3)';
end

function d = ch2int16(u8, swap)
if swap, u8 = u8([2 1]); end
d = double(u8);
d = d(1) + d(2)*256;
end

%% subfunction: remove trailing null and space for uint8
function ch = dcm_str(b)
while ~isempty(b) && b(end)==0, b(end) = []; end
ch = deblank(char(b));
end

%% subfunction: return value length, numel(b)=6
function [n, nvr] = val_len(vr, b, expl, swap)
len16 = 'AE AS AT CS DA DS DT FD FL IS LO LT PN SH SL SS ST TM UI UL US';
if ~expl % implicit, length irrevalent to vr (faked as CS)
    n = ch2int32(b(1:4), swap);
    nvr = 4; % bytes of VR
elseif ~isempty(strfind(len16, vr)) % length in uint16
    n = ch2int16(b(1:2), swap);
    nvr = 2;
else % length in uint32 and skip 2 bytes
    n = ch2int32(b(3:6), swap);
    nvr = 6;
end
if n==13, n = 10; end % ugly bug fix for some old dicom file
end

%% subfunction: read value, called by search method and read_item
function [dat, info] = read_val(b, vr, swap)
if strcmp(vr, 'DS') || strcmp(vr, 'IS')
    dat = sscanf(char(b), '%f\\'); % like 1\2\3
elseif ~isempty(strfind('AE AS CS DA DT LO LT PN SH ST TM UI UT', vr)) % char
    dat = dcm_str(b);
else % numeric data, UN. SQ taken care of
    fmt = vr2fmt(vr);
    if isempty(fmt)
        dat = [];
        info = sprintf('Given up: Invalid VR (%d %d)', vr);
        return;
    end
    dat = typecast(b, fmt)';
    if swap, dat = swapbytes(dat); end
end
info = '';
end

%% subfunction: numeric format str from VR
function fmt = vr2fmt(vr)
switch vr
    case 'US', fmt = 'uint16';
    case 'OB', fmt = 'uint8';
    case 'FD', fmt = 'double';
    case 'SS', fmt = 'int16';
    case 'UL', fmt = 'uint32';
    case 'SL', fmt = 'int32';
    case 'FL', fmt = 'single';
    case 'AT', fmt = 'uint16';
    case 'OW', fmt = 'uint16';
    case 'UN', fmt = 'uint8';
    otherwise, fmt = '';
end
end

%% subfunction: get nFrames into p.nFrames
function p = get_nFrames(s, p)
if isfield(s, 'NumberOfFrames')
    p.nFrames = s.NumberOfFrames; % useful for PerFrameSQ
elseif all(isfield(s, {'Columns' 'Rows' 'BitsAllocated'})) && p.bytes<4294967295
    if isfield(s, 'SamplesPerPixel'), spp = double(s.SamplesPerPixel);
    else spp = 1;
    end
    n = p.bytes * 8 / double(s.BitsAllocated);
    p.nFrames = n / (spp * double(s.Columns) * double(s.Rows));
else
    p.nFrames = nan;
end
end

%% subfunction: decode Siemens CSA image and series header
function csa = read_csa(csa)
b = csa';
if numel(b)<4 || ~strcmp(char(b(1:4)), 'SV10'), return; end % no op if not SV10
chDat = 'AE AS CS DA DT LO LT PN SH ST TM UI UN UT';
i = 8; % 'SV10' 4 3 2 1
try %#ok in case of error, we return the original csa
    nField = ch2int32(b(i+(1:4)), 0); i=i+8;
    for j = 1:nField
        i=i+68; % name(64) and vm(4)
        vr = char(b(i+(1:2))); i=i+8; % vr(4), syngodt(4)
        n = ch2int32(b(i+(1:4)), 0); i=i+8;
        if n<1, continue; end % skip name decoding, faster
        nam = char(b(i-84+(1:64)));
        ind = strfind(nam, char(0));
        nam = nam(1:ind(1)-1);
        % fprintf('%s %3g %s\n', vr, n, nam);
        
        dat = cell(n,1);
        for k = 1:n % n is often 6, but often only the first contains value
            len = ch2int32(b(i+(1:4)), 0); i=i+16;
            if len<1, i = i+(n-k)*16; break; end % rest are empty too
            foo = char(b(i+(1:len-1))); % exclude nul, need for Octave
            i = i + ceil(len/4)*4; % multiple 4-byte
            if isempty(strfind(chDat, vr))
                tmp = sscanf(foo, '%f', 1); % numeric to double
                if ~isempty(tmp), dat{k} = tmp; end
            else
                dat{k} = deblank(foo);
            end
        end
        if ischar(dat{1})
            ind = cellfun(@isempty, dat);
            dat = dat(~ind);
            if isempty(dat), continue; end
            if numel(dat)==1, dat = dat{1}; end
        else
            dat = [dat{:}]'; % cell2mat
        end
        rst.(nam) = dat;
    end
    csa = rst;
end
end

%% subfunction: decode GE ProtocolDataBlock
function ch = read_ProtocolDataBlock(ch)
n = typecast(ch(1:4), 'int32') + 4; % nBytes, zeros may be padded to make 4x
if ~all(ch(5:6) == [31 139]') || n>numel(ch), return; end % gz signature

b = gunzip_mem(ch(5:n));
b = regexp(char(b'), '(\w*)\s+"(.*)"', 'tokens', 'dotexceptnewline');
if isempty(b{1}), return; end % guzip faild or wrong format

try %#ok
    for j = 1:numel(b)
        a = sscanf(b{j}{2}, '%f', 1);
        if isempty(a), rst.(b{j}{1}) = strtrim(b{j}{2});
        else rst.(b{j}{1}) = a; % convert into num if possible
        end
    end
    ch = rst;
end
end

%% gunzip data in memory if possible.
% For a GE ProtocolDataBlock, memory / file approaches take 0.5 / 43 ms.
% When gz_bytes is large, pigz will be faster. The reversing point is about 8M.
function bytes = gunzip_mem(gz_bytes)
bytes = [];
try
    import com.mathworks.mlwidgets.io.*
    streamCopier = InterruptibleStreamCopier.getInterruptibleStreamCopier;
    baos = java.io.ByteArrayOutputStream;
    b = typecast(gz_bytes, 'int8');
    bais = java.io.ByteArrayInputStream(b);
    gzis = java.util.zip.GZIPInputStream(bais);
    streamCopier.copyStream(gzis, baos);
    bytes = typecast(baos.toByteArray, 'uint8'); % int8 to uint8
catch
    try %#ok
        tmp = tempname; % temp gz file
        fid = fopen([tmp '.gz'], 'w');
        if fid<0, return; end
        cln = onCleanup(@() delete([tmp '*'])); % delete gz and unziped files
        fwrite(fid, gz_bytes, 'uint8');
        fclose(fid);
        
        gunzipOS = nii_tool('func_handle', 'gunzipOS');
        gunzipOS([tmp '.gz']);
        
        fid = fopen(tmp);
        bytes = fread(fid, '*uint8');
        fclose(fid);
    end
end
end

%% subfunction: get fields for multiframe dicom
function s1 = search_MF_val(s, s1, iFrame)
% s1 = search_MF_val(s, s1, iFrame);
%  Arg 1: the struct returned by dicm_hdr for a multiframe dicom
%  Arg 2: a struct with fields to search, and with initial value, such as
%    zeros or nans. The number of rows indicate the number of values for the
%    tag, and columns for frames indicated by iFrame, Arg 3.
%  Arg 3: frame indice, length consistent with columns of s1 field values.
% Example:
%  s = dicm_hdr('multiFrameFile.dcm'); % read only 1,2 and last frame by default
%  s1 = struct('ImagePositionPatient', nan(3, s.NumberOfFrames)); % define size
%  s1 = search_MF_val(s, s1, 1:s.NumberOfFrames); % get values
%  This is MUCH faster than asking all frames by dicm_hdr, and avoid to get into
%  annoying SQ levels under PerFrameFuntionalGroupSequence.
% In case a tag is not found in PerFrameSQ, the code will search SharedSQ and
% common tags, and will ignore the 3th arg and fake the same value for all
% frames.

if ~isfield(s, 'PerFrameFunctionalGroupsSequence'), return; end
expl = false;
be = false;
if isfield(s, 'TransferSyntaxUID')
    expl = ~strcmp(s.TransferSyntaxUID, '1.2.840.10008.1.2');
    be =    strcmp(s.TransferSyntaxUID, '1.2.840.10008.1.2.2');
end

fStart = s.PerFrameFunctionalGroupsSequence.FrameStart; % error if no FrameStart
fid = fopen(s.Filename);
b0 = fread(fid, fStart(1), 'uint8=>char')'; % before 1st frame in PerFrameSQ
b = fread(fid, s.PixelData.Start-fStart(1), 'uint8=>char')'; % within PerFrameSQ
fclose(fid);

fStart(end+1) = s.PixelData.Start; % for later ind search
fStart = fStart - fStart(1) + 1; % 1-based index in b

flds = fieldnames(s1);
dict = dicm_dict(s.Manufacturer, flds);
len16 = 'AE AS AT CS DA DS DT FD FL IS LO LT PN SH SL SS ST TM UI UL US';
chDat = 'AE AS CS DA DT LO LT PN SH ST TM UI UT';
nf = numel(iFrame);

for i = 1:numel(flds)
    k = find(strcmp(dict.name, flds{i}), 1, 'last'); % GE has another ipp tag
    vr = dict.vr{k};
    group = dict.group(k);
    isBE = be && group~=2;
    isEX = expl || group==2;
    tg = char(typecast([group dict.element(k)], 'uint8'));
    if isBE, tg = tg([2 1 4 3]); end
    if isEX, tg = [tg vr]; end %#ok
    ind = strfind(b, tg);
    ind = ind(mod(ind,2)>0); % indice is odd
    if isempty(ind) % no tag in PerFrameSQ, try tag before PerFrameSQ
        ind = strfind(b0, tg);
        ind = ind(mod(ind,2)>0);
        if ~isempty(ind)
            k = ind(1) + numel(tg); % take 1st in case of multiple
            [n, nvr] = val_len(vr, uint8(b0(k+(0:5))), isEX, isBE); k = k + nvr;
            a = read_val(uint8(b0(k+(0:n-1))), vr, isBE);
            if ischar(a), a = {a}; end
            s1.(flds{i}) = repmat(a, 1, nf); % all frames share the same value
        end
        continue;
    end
    
    len = 4; % bytes of tag value length (uint32)
    if ~isEX % implicit, length irrevalent to VR
        ind = ind + 4; % tg(4)
    elseif ~isempty(strfind(len16, vr)) % data length in uint16
        ind = ind + 6; % tg(4), vr(2)
        len = 2;
    else % length in uint32: skip 2 bytes
        ind = ind + 8; % tg(4), vr(2), skip(2)
    end
    
    isCH = ~isempty(strfind(chDat, vr)); % char data
    isDS = strcmp(vr, 'DS') || strcmp(vr, 'IS');
    if ~isCH && ~isDS % numeric data, UN or SQ
        fmt = vr2fmt(vr);
        if isempty(fmt), continue; end % skip SQ
    end
    
    for k = 1:nf
        j = iFrame(k); % asked frame index
        j = find(ind>fStart(j) & ind<fStart(j+1), 1); % index in ind
        if isempty(j), continue; end % no tag for this frame
        if len==2, n = ch2int16(b(ind(j)+(0:1)), isBE);
        else       n = ch2int32(b(ind(j)+(0:3)), isBE);
        end
        a = b(ind(j)+len+(0:n-1));
        if isDS
            a = sscanf(a, '%f\\'); % like 1\2\3
            try s1.(flds{i})(:,k) = a; catch, end % ignore in case of error
        elseif isCH
            while ~isempty(a) && a(end)==0, a(end) = ''; end
            try s1.(flds{i}){k} = deblank(a); catch, end
        else
            a = typecast(uint8(a), fmt)';
            if isBE, a = swapbytes(a); end
            try s1.(flds{i})(:,k) = a; catch, end
        end
    end
end
end

%% subfunction: read PAR file, return struct like that from dicm_hdr.
function [s, err] = philips_par(fname)
err = '';
if numel(fname)>4 && strcmpi(fname(end+(-3:0)), '.REC')
    fname(end+(-3:0)) = '.PAR';
    if ~exist(fname, 'file'), fname(end+(-3:0)) = '.par'; end
end
fid = fopen(fname);
if fid<0, s = []; err = ['File not exist: ' fname]; return; end
str = fread(fid, inf, '*char')'; % read all as char
fname = fopen(fid); % name with full path
fclose(fid);

str = strrep(str, char(13), char(10)); % make carriage return single char(10)
while true
    ind = strfind(str, char([10 10]));
    if isempty(ind), break; end
    str(ind) = [];
end
n = numel(str);
while true
    str = strrep(str, [char(10) '.  '], [char(10) '. ']);
    n1 = numel(str);
    if n1==n, break; else n = n1; end
end

% In V4, offcentre and Angulation labeled as y z x, but actually x y z. We
% try not to use these info
key = 'image export tool';
C = regexp(str, [key '\s*(.*)'], 'tokens', 'once', 'dotexceptnewline');
if isempty(C), err = 'Not PAR file'; s = []; return; end
s.SoftwareVersion = [C{1} '\PAR'];
if strncmpi(s.SoftwareVersion, 'V3', 2)
    fprintf(2, ' V3 PAR file is not supported.\n');
    s = []; return;
end

s.PatientName = par_key('Patient name', 0);
s.StudyDescription = par_key('Examination name', 0);
[pth, nam] = fileparts(fname);
s.SeriesDescription = nam;
s.ProtocolName = par_key('Protocol name', 0);
foo = par_key('Examination date/time', 0);
foo = foo(isstrprop(foo, 'digit'));
s.AcquisitionDateTime = foo;
s.SeriesNumber = par_key('Acquisition nr');
s.SeriesInstanceUID = sprintf('%g.%s.%09.0f', s.SeriesNumber, ...
    datestr(now, 'yymmdd.HHMMSS.fff'), rand*1e9);
% s.ReconstructionNumberMR = par_key('Reconstruction nr');
% s.MRSeriesScanDuration = par_key('Scan Duration');
s.NumberOfEchoes = par_key('Max. number of echoes');
nSL = par_key('Max. number of slices/locations');
s.LocationsInAcquisition = nSL;
foo = par_key('Patient position', 0);
if isempty(foo), foo = par_key('Patient Position', 0); end
if ~isempty(foo)
    if numel(foo)>4, s.PatientPosition = foo(regexp(foo, '\<.'));
    else s.PatientPosition = foo;
    end
end
s.MRAcquisitionType = par_key('Scan mode', 0);
s.ScanningSequence = par_key('Technique', 0); % ScanningTechnique
typ = par_key('Series Type', 0); typ(isspace(typ)) = '';
s.ImageType = ['PhilipsPAR\' typ '\' s.ScanningSequence];
s.RepetitionTime = par_key('Repetition time');
s.WaterFatShift = par_key('Water Fat shift');
rotAngle = par_key('Angulation midslice'); % (ap,fh,rl) deg
rotAngle = rotAngle([3 1 2]);
posMid = par_key('Off Centre midslice'); % (ap,fh,rl) [mm]
s.Stack.Item_1.MRStackOffcentreAP = posMid(1);
s.Stack.Item_1.MRStackOffcentreFH = posMid(2);
s.Stack.Item_1.MRStackOffcentreRL = posMid(3);
posMid = posMid([3 1 2]); % better precision than those in the table
s.EPIFactor = par_key('EPI factor');
% s.DynamicSeries = par_key('Dynamic scan'); % 0 or 1
isDTI = par_key('Diffusion')>0;
if isDTI
    s.ImageType = [s.ImageType '\DIFFUSION\'];
    s.DiffusionEchoTime = par_key('Diffusion echo time'); % ms
end

foo = par_key('Preparation direction', 0); % Anterior-Posterior
if ~isempty(foo)
    foo = foo(regexp(foo, '\<.')); % 'AP'
    s.Stack.Item_1.MRStackPreparationDirection = foo;
    iPhase = strfind('LRAPFH', foo(1));
    iPhase = ceil(iPhase/2); % 1/2/3
end

% Get list of para meaning for the table, and col index of each para
i1 = strfind(str, 'IMAGE INFORMATION DEFINITION'); i1 = i1(end);
ind = strfind(str(i1:end), [char(10) '#']) + i1;
for i = 1:9 % find the empty line before column descrip
    [~, foo] = strtok(str(ind(i):ind(i+1)-2)); % remove # and char(10)
    if isempty(foo), break; end
end
j = 1;
for i = i+1:numel(ind)
    [~, foo] = strtok(str(ind(i):ind(i+1)-2));
    if isempty(foo), break; end % the end of the col label
    foo = strtrim(foo);
    i3 = strfind(foo, '<');
    i2 = strfind(foo, '(');
    if isempty(i3), i3 = i2(1); end
    colLabel{j} = strtrim(foo(1:i3(1)-1)); %#ok para name
    nCol = sscanf(foo(i2(end)+1:end), '%g');
    if isempty(nCol), nCol = 1; end
    iColumn(j) = nCol; %#ok number of columns in the table for this para
    j = j + 1;
end
iColumn = cumsum([1 iColumn]); % col start ind for corresponding colLabel
keyInLabel = @(key)strcmp(colLabel, key);
colIndex = @(key)iColumn(keyInLabel(key));

i1 = strfind(str, '= IMAGE INFORMATION ='); i1 = i1(end);
ind = strfind(str(i1:end), char(10)) + i1 + 1; % start of a line
for i = 1:9
    foo = sscanf(str(ind(i):end), '%g', 1);
    if ~isempty(foo), break; end % get the first number
end
while str(ind(i))==10, i = i+1; end % skip empty lines (only one)
str = str(ind(i):end); % now start of the table
i1 = strfind(str, char(10));
para = sscanf(str(1:i1(1)), '%g'); % 1st row
n = numel(para); % number of items each row, 41 for V4
para = sscanf(str, '%g'); % read all numbers
nImg = floor(numel(para) / n);
para = reshape(para(1:n*nImg), n, nImg)'; % whole table now
getTableVal('index in REC file', 'IndexInREC', 1:nImg);
if isfield(s, 'IndexInREC') % why Philips includes this?
    if ~all(diff(s.IndexInREC) == 1) % not needed, just avoid accident
        para = para(s.IndexInREC, :); % in the order of REC img
    end
    s = rmfield(s, 'IndexInREC');
end

s.NumberOfFrames = nImg;
nVol = nImg/nSL;
s.NumberOfTemporalPositions = nVol;

s.Dim3IsVolume = (diff(para(1:2, colIndex('slice number'))) == 0);
if s.Dim3IsVolume
    iVol = 1:nVol;
    iSL = 1:nVol:nImg;
else
    iVol = 1:nSL:nImg;
    iSL = 1:nSL;
end

% PAR/REC file may not start with SliceNumber of 1, WHY?
sl = para(iSL, colIndex('slice number'));
if any(diff(sl,2)>0), s.SliceNumber = sl; end % slice order in REC file

imgType = para(iVol, colIndex('image_type_mr')); % 0 mag; 3, phase?
if any(diff(imgType) ~= 0) % more than 1 type of image
    s.ComplexImageComponent = 'MIXED';
    s.VolumeIsPhase = (imgType==3); % one for each vol
    s.LastFile.RescaleIntercept = para(end, colIndex('rescale intercept'));
    s.LastFile.RescaleSlope = para(end, colIndex('rescale slope'));
elseif imgType(1)==0, s.ComplexImageComponent = 'MAGNITUDE';
elseif imgType(1)==3, s.ComplexImageComponent = 'PHASE';
end

% These columns should be the same for nifti-convertible images:
cols = {'image pixel size' 'recon resolution' 'image angulation' ...
    'slice thickness' 'slice gap' 'slice orientation' 'pixel spacing'};
if ~strcmp(s.ComplexImageComponent, 'MIXED')
    cols = [cols {'rescale intercept' 'rescale slope'}];
end
ind = [];
for i = 1:numel(cols)
    j = find(keyInLabel(cols{i}));
    if isempty(j), continue; end
    ind = [ind iColumn(j):iColumn(j+1)-1]; %#ok
end
foo = para(:, ind);
foo = abs(diff(foo));
if any(foo(:) > 1e-5)
    err = sprintf('Inconsistent image size, bits etc: %s', fname);
    fprintf(2, ' %s. \n', err);
    s = []; return;
end

% getTableVal('echo number', 'EchoNumber', 1:nImg);
% getTableVal('dynamic scan number', 'TemporalPositionIdentifier', 1:nImg);
getTableVal('image pixel size', 'BitsAllocated');
getTableVal('recon resolution', 'Columns');
s.Rows = s.Columns(2); s.Columns = s.Columns(1);
getTableVal('rescale intercept', 'RescaleIntercept');
getTableVal('rescale slope', 'RescaleSlope');
getTableVal('window center', 'WindowCenter', 1:nImg);
getTableVal('window width', 'WindowWidth', 1:nImg);
mx = max(s.WindowCenter + s.WindowWidth/2);
mn = min(s.WindowCenter - s.WindowWidth/2);
s.WindowCenter = round((mx+mn)/2);
s.WindowWidth = ceil(mx-mn);
getTableVal('slice thickness', 'SliceThickness');
getTableVal('echo_time', 'EchoTime');
getTableVal('image_flip_angle', 'FlipAngle');
getTableVal('number of averages', 'NumberOfAverages');
% getTableVal('trigger_time', 'TriggerTime', 1:nImg);
% getTableVal('dyn_scan_begin_time', 'TimeOfAcquisition', 1:nImg);
if isDTI
    getTableVal('diffusion_b_factor', 'B_value', iVol);
    fld = 'bvec_original';
    getTableVal('diffusion', fld, iVol);
    if isfield(s, fld), s.(fld) = s.(fld)(:, [3 1 2]); end
end
getTableVal('TURBO factor', 'TurboFactor');

% Rotation order and signs are figured out by try and err, not 100% sure
ca = cosd(rotAngle); sa = sind(rotAngle);
rx = [1 0 0; 0 ca(1) -sa(1); 0 sa(1) ca(1)]; % 3D rotation
ry = [ca(2) 0 sa(2); 0 1 0; -sa(2) 0 ca(2)];
rz = [ca(3) -sa(3) 0; sa(3) ca(3) 0; 0 0 1];
R = rx * ry * rz; % seems right for Philips

getTableVal('slice orientation', 'SliceOrientation'); % 1/2/3 for TRA/SAG/COR
iOri = mod(s.SliceOrientation+1, 3) + 1; % [1 2 3] to [3 1 2]
if iOri == 1 % Sag
    s.SliceOrientation = 'SAGITTAL';
    R(:,[1 3]) = -R(:,[1 3]);
    R = R(:, [2 3 1]);
elseif iOri == 2 % Cor
    s.SliceOrientation = 'CORONAL';
    R(:,3) = -R(:,3);
    R = R(:, [1 3 2]);
else % Tra
    s.SliceOrientation = 'TRANSVERSAL';
end

% 'pixel spacing' 'pixel spacing' and 'slice gap' have poor precision for v<=4?
% It may be wrong to use FOV, maybe due to partial Fourier?
getTableVal('pixel spacing', 'PixelSpacing');
s.PixelSpacing = s.PixelSpacing(:);
getTableVal('slice gap', 'SpacingBetweenSlices');
s.SpacingBetweenSlices = s.SpacingBetweenSlices + s.SliceThickness;

if exist('iPhase', 'var')
    foo = 'COL';
    if iPhase == (iOri==1)+1, foo = 'ROW'; end
    s.InPlanePhaseEncodingDirection = foo;
end

s.ImageOrientationPatient = R(1:6)';
R = R * diag([s.PixelSpacing; s.SpacingBetweenSlices]);
R = [R posMid; 0 0 0 1]; % 4th col is mid slice center position

getTableVal('image offcentre', 'SliceLocation');
s.SliceLocation = s.SliceLocation(iOri); % center loc for 1st slice
if sign(R(iOri,3)) ~= sign(posMid(iOri)-s.SliceLocation)
    R(:,3) = -R(:,3);
end

R(:,4) = R * [-([s.Columns s.Rows nSL]-1)/2 1]'; % vol center to corner of 1st
y = R * [0 0 nSL-1 1]'; % last slice position
s.ImagePositionPatient = R(1:3,4);
s.LastFile.ImagePositionPatient = y(1:3);
s.Manufacturer = 'Philips';
s.Filename = fullfile(pth, [nam '.REC']); % rest for dicm_img
s.PixelData.Start = 0;
s.PixelData.Bytes = s.Rows * s.Columns * nImg * s.BitsAllocated / 8;

% nested function: set field if the key is in colTable
    function getTableVal(key, fldname, iRow)
        if nargin<3, iRow = 1; end
        iCol = find(keyInLabel(key));
        if isempty(iCol), return; end
        s.(fldname) = para(iRow, iColumn(iCol):iColumn(iCol+1)-1);
    end

% nested subfunction: return value specified by key in PAR file
    function val = par_key(key, isNum)
        if nargin<2, isNum = true; end
        i0 = strfind(str, [char(10) '. ' key]); % start with new line
        if isempty(i0)
            if isNum, val = []; else val = ''; end
            return;
        end
        i0 = i0(1) + numel(key);
        ii = strfind(str(i0:end), char(10)); % next new line
        a = str(i0+(1:ii(1)-2)); % like '   : 5'
        i0 = strfind(a, ':'); % must have ':'
        val = strtrim(a(i0(1)+1:end));
        if isNum, val = sscanf(val, '%g'); end
    end
end

%% subfunction: read AFNI HEAD file, return struct like that from dicm_hdr.
function [s, err] = afni_head(fname)
persistent SN;
if isempty(SN), SN = 1; end
err = '';
if numel(fname)>5 && strcmp(fname(end+(-4:0)), '.BRIK')
    fname(end+(-4:0)) = '.HEAD';
end
fid = fopen(fname);
if fid<0, s = []; err = ['File not exist: ' fname]; return; end
str = fread(fid, inf, '*char')';
fname = fopen(fid);
fclose(fid);

i = strfind(str, 'DATASET_DIMENSIONS');
if isempty(i), s = []; err = 'Not brik header file'; return; end

% these make dicm_nii.m happy
[~, foo] = fileparts(fname);
% s.IsAFNIHEAD = true;
s.ProtocolName = foo;
s.SeriesNumber = SN; SN = SN+1; % make it unique for multilple files
s.SeriesInstanceUID = sprintf('%g.%s.%09.0f', s.SeriesNumber, ...
    datestr(now, 'yymmdd.HHMMSS.fff'), rand*1e9);
s.ImageType = ['AFNIHEAD\' afni_key('TYPESTRING')];

foo = afni_key('BYTEORDER_STRING');
if strcmp(foo(1), 'M'), err = 'BYTEORDER_STRING not supported'; s = []; return; end

foo = afni_key('BRICK_FLOAT_FACS');
if any(diff(foo)~=0), err = 'Inconsistent BRICK_FLOAT_FACS';
    s = []; return;
end
if foo(1)==0, foo = 1; end
s.RescaleSlope = foo(1);
s.RescaleIntercept = 0;

foo = afni_key('BRICK_TYPES');
if any(diff(foo)~=0), err = 'Inconsistent DataType'; s = []; return; end
foo = foo(1);
if foo == 0
    s.BitsAllocated =  8; s.PixelData.Format = '*uint8';
elseif foo == 1
    s.BitsAllocated = 16; s.PixelData.Format = '*int16';
elseif foo == 3
    s.BitsAllocated = 32; s.PixelData.Format = '*single';
else
    error('Unsupported BRICK_TYPES: %g', foo);
end

hist = afni_key('HISTORY_NOTE');
i = strfind(hist, 'Time:') + 6;
if ~isempty(i)
    dat = sscanf(hist(i:end), '%11c', 1); % Mar  1 2010
    dat = datenum(dat, 'mmm dd yyyy');
    s.AcquisitionDateTime = datestr(dat, 'yyyymmdd');
end
i = strfind(hist, 'Sequence:') + 9;
if ~isempty(i), s.ScanningSequence = strtok(hist(i:end), ' '); end
i = strfind(hist, 'Studyid:') + 8;
if ~isempty(i), s.StudyID = strtok(hist(i:end), ' '); end
% i = strfind(hist, 'Dimensions:') + 11;
% if ~isempty(i)
%     dimStr = strtok(hist(i:end), ' ') % 64x64x35x92
% end
% i = strfind(hist, 'Orientation:') + 12;
% if ~isempty(i)
%     oriStr = strtok(hist(i:end), ' ') % LAI
% end
i = strfind(hist, 'TE:') + 3;
if ~isempty(i), s.EchoTime = sscanf(hist(i:end), '%g', 1) * 1000; end

% foo = afni_key('TEMPLATE_SPACE'); % ORIG/TLRC
% INT_CMAP
foo = afni_key('SCENE_DATA');
s.TemplateSpace = foo(1)+1; %[0] 0=+orig, 1=+acpc, 2=+tlrc
if foo(2)==9, s.ImageType = [s.ImageType '\DIFFUSION\']; end
% ori = afni_key('ORIENT_SPECIFIC')+1;
% orients = [1 -1 -2 2 3 -3]; % RL LR PA AP IS SI
% ori = orients(ori) % in dicom/afni LPS,
% seems always [1 2 3], meaning AFNI re-oriented the volome

% no read/phase/slice dim info, so following 3D info are meaningless
dim = afni_key('DATASET_DIMENSIONS');
s.Columns = dim(1); s.Rows = dim(2); s.LocationsInAcquisition = dim(3);
R = afni_key('IJK_TO_DICOM_REAL'); % IJK_TO_DICOM is always straight?
if isempty(R), R = afni_key('IJK_TO_DICOM'); end
R = reshape(R, [4 3])';
s.ImagePositionPatient = R(:,4); % afni_key('ORIGIN') can be wrong
y = [R; 0 0 0 1] * [0 0 dim(3)-1 1]';
s.LastFile.ImagePositionPatient = y(1:3);
R = R(1:3, 1:3);
R = R ./ (ones(3,1) * sqrt(sum(R.^2)));
s.ImageOrientationPatient = R(1:6)';
foo = afni_key('DELTA');
s.PixelSpacing = foo(1:2);
% s.SpacingBetweenSlices = foo(3);
s.SliceThickness = foo(3);
foo = afni_key('BRICK_STATS');
foo = reshape(foo, [2 numel(foo)/2]);
mn = min(foo(1,:)); mx = max(foo(2,:));
s.WindowCenter = (mx+mn)/2;
s.WindowWidth = mx-mn;
foo = afni_key('TAXIS_FLOATS'); %[0]:0;
if ~isempty(foo), s.RepetitionTime = foo(2)*1000; end

foo = afni_key('TAXIS_NUMS'); % [0]:nvals; [1]: 0 or nSL normally
if ~isempty(foo)
    inMS = foo(3)==77001;
    foo = afni_key('TAXIS_OFFSETS');
    if inMS, foo = foo/1000; end
    if ~isempty(foo), s.MosaicRefAcqTimes = foo; end
end

foo = afni_key('DATASET_RANK'); % [3 nvals]
dim(4) = foo(2);
s.NumberOfTemporalPositions = dim(4);
% s.NumberOfFrames = dim(4)*dim(3);

s.Manufacturer = '';
s.Filename = strrep(fname, '.HEAD', '.BRIK');
s.PixelData.Start = 0; % make it work for dicm_img.m
s.PixelData.Bytes = prod(dim(1:4)) * s.BitsAllocated / 8;

% subfunction: return value specified by key in afni header str
    function val = afni_key(key)
        i1 = regexp(str, ['\nname\s*=\s*' key '\n']); % line 'name = key'
        if isempty(i1), val = []; return; end
        i1 = i1(1) + 1;
        i2 = regexp(str(1:i1), 'type\s*=\s*\w*-attribute\n');
        keyType = sscanf(str(i2(end):i1), 'type%*c=%*c%s', 1); %'string-attribute'
        i1 = find(str(i1:end)==char(10), 1, 'first') + i1;
        count = sscanf(str(i1:end), 'count%*c=%*c%g', 1);
        if strcmp(keyType, 'string-attribute')
            i1 = find(str(i1:end)=='''', 1, 'first') + i1;
            val = str(i1+(0:count-2));
        else
            i1 = find(str(i1:end)==char(10), 1, 'first') + i1;
            val = sscanf(str(i1:end), '%g', count);
        end
    end
end

%% Subfunction: read BrainVoyager vmr/fmr/dmr. Call BVQXfile
function [s, err] = bv_file(fname)
s = []; err = '';
try
    bv = BVQXfile(fname);
catch me
    err = me.message;
    if strfind(me.identifier, 'UndefinedFunction')
        fprintf(2, 'Please download BVQXtools at \n%s\n', ...
            'http://support.brainvoyager.com/available-tools/52-matlab-tools-bvxqtools.html');
    end
    return;
end

if ~isempty(bv.Trf)
    for i = 1:numel(bv.Trf)
        if ~isequal(diag(bv.Trf(i).TransformationValues), [1 1 1 1]')
            err = 'Data has been transformed: skipped.';
            return;
        end
    end
end

persistent SN subj folder % folder is used to update subj
if isempty(SN), SN = 1; subj = ''; folder = ''; end
s.Filename = bv.FilenameOnDisk;
fType = bv.filetype;
s.ImageType = ['BrainVoyagerFile\' fType];

% Find a fmr/dmr, and get subj based on dicom file name in BV format.
% Suppose BV files in the folder are for the same subj
[pth, nam] = fileparts(s.Filename);
s.SeriesDescription = nam;
if isempty(folder) || ~strcmp(folder, pth)
    folder = pth;
    subj = '';
    if strcmp(fType, 'fmr') || strcmp(fType, 'dmr')
        [~, nam] = fileparts(bv.FirstDataSourceFile);
        nam = strtok(nam, '-');
        if ~isempty(nam), subj = nam; end
    else
        fnames = dir([pth '/*.fmr']);
        if isempty(fnames), fnames = dir([pth '/*.dmr']); end
        if ~isempty(fnames)
            bv1 = BVQXfile(fullfile(pth, fnames(1).name));
            [~, nam] = fileparts(bv1.FirstDataSourceFile);
            bv1.ClearObject;
            nam = strtok(nam, '-');
            if ~isempty(nam), subj = nam; end
        end
    end
end
if ~isempty(subj), s.PatientName = subj; end

s.SoftwareVersion = sprintf('%g/BV_FileVersion', bv.FileVersion);
s.Columns = bv.NCols;
s.Rows = bv.NRows;
s.SliceThickness = bv.SliceThickness;
R = [bv.RowDirX bv.RowDirY bv.RowDirZ; bv.ColDirX bv.ColDirY bv.ColDirZ]';
s.ImageOrientationPatient = R(:);
R(:,3) = cross(R(:,1), R(:,2));
[~, ixyz] = max(abs(R)); iSL =ixyz(3);

try
    s.TemplateSpace = bv.ReferenceSpace; % 0/2/3: Scanner/ACPC/TAL
    if s.TemplateSpace==0, s.TemplateSpace = 1; end
catch
    s.TemplateSpace = 1;
end
pos = [bv.Slice1CenterX bv.Slice1CenterY bv.Slice1CenterZ
    bv.SliceNCenterX bv.SliceNCenterY bv.SliceNCenterZ]'; % for real slices

if strcmpi(fType, 'vmr')
    s.SpacingBetweenSlices = s.SliceThickness + bv.GapThickness;
    s.PixelSpacing = [bv.VoxResX bv.VoxResY]';
    if ~isempty(bv.VMRData16)
        nSL = bv.DimZ;
        s.PixelData = bv.VMRData16; % no padded zeros
    else
        v16 = [s.Filename(1:end-3) 'v16'];
        if exist(v16, 'file')
            bv16 = BVQXfile(v16);
            nSL = bv16.DimZ;
            s.PixelData = bv16.VMRData; % no padded zeros
            bv16.ClearObject;
        else % fall back the 8-bit data, and deal with padded zeros
            ix = floor((bv.DimX - s.Columns)/2);
            iy = floor((bv.DimY - s.Rows)/2);
            R3 = abs(R(iSL,3)) * s.SpacingBetweenSlices;
            nSL = round(abs(diff(pos(iSL,:))) / R3) + 1;
            iz = floor((bv.DimZ - nSL)/2);
            s.PixelData = bv.VMRData(ix+(1:s.Columns), iy+(1:s.Rows), iz+(1:nSL), :);
        end
    end
    s.LocationsInAcquisition = nSL;
    s.MRAcquisitionType = '3D'; % for dicm2nii to re-orient
elseif strcmpi(fType, 'fmr') || strcmpi(fType, 'dmr')
    s.SpacingBetweenSlices = s.SliceThickness + bv.SliceGap;
    s.PixelSpacing = [bv.InplaneResolutionX bv.InplaneResolutionY]';
    nSL = bv.NrOfSlices;
    s.LocationsInAcquisition = nSL;
    s.NumberOfTemporalPositions = bv.NrOfVolumes;
    s.RepetitionTime = bv.TR;
    s.EchoTime = bv.TE;
    if bv.TimeResolutionVerified
        switch bv.SliceAcquisitionOrder % the same as NIfTI?
            case 1, ind = 1:nSL;
            case 2, ind = nSL:-1:1;
            case 3, ind = [1:2:nSL 2:2:nSL];
            case 4, ind = [nSL:-2:1 nSL-1:-2:1];
            case 5, ind = [2:2:nSL 1:2:nSL];
            case 6, ind = [nSL-1:-2:1 nSL:-2:1];
            otherwise, ind = []; err = 'Unknown SliceAcquisitionOrder';
        end
        if ~isempty(ind)
            t = (0:s.LocationsInAcquisition-1)' * bv.InterSliceTime; % ms
            t(ind) = t;
            s.SliceTiming = t;
        end
    end
    if strcmpi(fType, 'fmr')
        bv.LoadSTC;
        s.PixelData = permute(bv.Slice(1).STCData , [1 2 4 3]);
        for i = 2:numel(bv.Slice)
            s.PixelData(:,:,i,:) = permute(bv.Slice(i).STCData , [1 2 4 3]);
        end
    else % dmr
        s.ImageType = [s.ImageType '\DIFFUSION\'];
        bv.LoadDWI;
        s.PixelData = bv.DWIData;
        if strncmpi(bv.GradientInformationAvailable, 'Y', 1)
            a = bv.GradientInformation; % nDir by 4
            s.B_value = a(:,4);
            a = a(:,1:3); % bvec
            % Following should be right in theory, but I would trust the grd
            % table which should be in dicom coodinate system, rather than the
            % confusing Gradient?DirInterpretation
            %             % 1:6 for LR RL AP PA IS SI. Default [2 3 5] by dicom LPS
            %             i1_6 = [bv.GradientXDirInterpretation ...
            %                     bv.GradientYDirInterpretation ...
            %                     bv.GradientZDirInterpretation];
            %             [xyz, ind] = sort(i1_6);
            %             if isequal(ceil(xyz/2), 1:3) % perm of 1/2/3
            %                 a = a(:,ind);
            %                 flip = xyz == [1 4 6]; % negative by dicom
            %                 a(:,flip) = -a(:,flip);
            %             else
            %                 str = sprintf(['Wrong Interpretation of gradient found: %s\n' ...
            %                        'Please check bvec and its sign.\n'], fname);
            %                 fprintf(2, str);
            %                 err = [err str];
            %             end
            s.bvec_original = a;
        end
    end
    
    % fmr/dmr are normally converted from uint16 to single
    if isfloat(s.PixelData) && isequal(floor(s.PixelData), s.PixelData) ...
            && max(s.PixelData(:))<32768 && min(s.PixelData(:))>=-32768
        s.PixelData = int16(s.PixelData);
    end
else
    err = ['Unknown BV file type: ' fType];
    s = [];
    return;
end

pos = pos - R(:,1:2) * diag(s.PixelSpacing) * [s.Columns s.Rows]'/2 * [1 1];
s.ImagePositionPatient = pos(:,1);
s.LastFile.ImagePositionPatient = pos(:,2);

% Following make dicm2nii happy
try %#ok
    [~, nam] = fileparts(bv.FirstDataSourceFile);
    serN = sscanf(nam, [subj '-%f'], 1);
    if ~isempty(serN), SN = serN; end
end
s.SeriesNumber = SN; SN = SN+1; % make it unique for multilple files
s.SeriesInstanceUID = sprintf('%g.%s.%09.0f', s.SeriesNumber, ...
    datestr(now, 'yymmdd.HHMMSS.fff'), rand*1e9);
c = class(s.PixelData);
if strcmp(c, 'double') %#ok
    s.BitsAllocated = 64;
elseif strcmp(c, 'single') %#ok
    s.BitsAllocated = 32;
else
    ind = find(isstrprop(c, 'digit'), 1);
    s.BitsAllocated = sscanf(c(ind:end), '%g');
end

end

%% Subfunction: clean up Unknown fields by LUT
function newName = dicom_lookup(oldName)

switch oldName
    case 'Unknown_0000_0000'
        newName = 'GroupLength';
    case 'Unknown_0000_0001'
        newName = 'CommandLengthToEnd';
    case 'Unknown_0000_0002'
        newName = 'AffectedSOPClassUID';
    case 'Unknown_0000_0003'
        newName = 'RequestedSOPClassUID';
    case 'Unknown_0000_0010'
        newName = 'CommandRecognitionCode';
    case 'Unknown_0000_0100'
        newName = 'CommandField';
    case 'Unknown_0000_0110'
        newName = 'MessageID';
    case 'Unknown_0000_0120'
        newName = 'MessageIDBeingRespondedTo';
    case 'Unknown_0000_0200'
        newName = 'Initiator';
    case 'Unknown_0000_0300'
        newName = 'Receiver';
    case 'Unknown_0000_0400'
        newName = 'FindLocation';
    case 'Unknown_0000_0600'
        newName = 'MoveDestination';
    case 'Unknown_0000_0700'
        newName = 'Priority';
    case 'Unknown_0000_0800'
        newName = 'DataSetType';
    case 'Unknown_0000_0850'
        newName = 'NumberOfMatches';
    case 'Unknown_0000_0860'
        newName = 'ResponseSequenceNumber';
    case 'Unknown_0000_0900'
        newName = 'Status';
    case 'Unknown_0000_0901'
        newName = 'OffendingElement';
    case 'Unknown_0000_0902'
        newName = 'ErrorComment';
    case 'Unknown_0000_0903'
        newName = 'ErrorID';
    case 'Unknown_0000_1000'
        newName = 'AffectedSOPInstanceUID';
    case 'Unknown_0000_1001'
        newName = 'RequestedSOPInstanceUID';
    case 'Unknown_0000_1002'
        newName = 'EventTypeID';
    case 'Unknown_0000_1005'
        newName = 'AttributeIdentifierList';
    case 'Unknown_0000_1008'
        newName = 'ActionTypeID';
    case 'Unknown_0000_1020'
        newName = 'NumberOfRemainingSuboperations';
    case 'Unknown_0000_1021'
        newName = 'NumberOfCompletedSuboperations';
    case 'Unknown_0000_1022'
        newName = 'NumberOfFailedSuboperations';
    case 'Unknown_0000_1023'
        newName = 'NumberOfWarningSuboperations';
    case 'Unknown_0000_1030'
        newName = 'MoveOriginatorApplicationEntityTitle';
    case 'Unknown_0000_1031'
        newName = 'MoveOriginatorMessageID';
    case 'Unknown_0000_4000'
        newName = 'DialogReceiver';
    case 'Unknown_0000_4010'
        newName = 'TerminalType';
    case 'Unknown_0000_5010'
        newName = 'MessageSetID';
    case 'Unknown_0000_5020'
        newName = 'EndMessageSet';
    case 'Unknown_0000_5110'
        newName = 'DisplayFormat';
    case 'Unknown_0000_5120'
        newName = 'PagePositionID';
    case 'Unknown_0000_5130'
        newName = 'TextFormatID';
    case 'Unknown_0000_5140'
        newName = 'NormalReverse';
    case 'Unknown_0000_5150'
        newName = 'AddGrayScale';
    case 'Unknown_0000_5160'
        newName = 'Borders';
    case 'Unknown_0000_5170'
        newName = 'Copies';
    case 'Unknown_0000_5180'
        newName = 'OldMagnificationType';
    case 'Unknown_0000_5190'
        newName = 'Erase';
    case 'Unknown_0000_51A0'
        newName = 'Print';
    case 'Unknown_0000_51B0'
        newName = 'Overlays';
    case 'Unknown_0002_0000'
        newName = 'FileMetaInformationGroupLength';
    case 'Unknown_0002_0001'
        newName = 'FileMetaInformationVersion';
    case 'Unknown_0002_0002'
        newName = 'MediaStorageSOPClassUID';
    case 'Unknown_0002_0003'
        newName = 'MediaStorageSOPInstanceUID';
    case 'Unknown_0002_0010'
        newName = 'TransferSyntaxUID';
    case 'Unknown_0002_0012'
        newName = 'ImplementationClassUID';
    case 'Unknown_0002_0013'
        newName = 'ImplementationVersionName';
    case 'Unknown_0002_0016'
        newName = 'SourceApplicationEntityTitle';
    case 'Unknown_0002_0100'
        newName = 'PrivateInformationCreatorUID';
    case 'Unknown_0002_0102'
        newName = 'PrivateInformation';
    case 'Unknown_0004_0000'
        newName = 'FileSetGroupLength';
    case 'Unknown_0004_1130'
        newName = 'FileSetID';
    case 'Unknown_0004_1141'
        newName = 'FileSetDescriptorFileID';
    case 'Unknown_0004_1142'
        newName = 'FileSetCharacterSet';
    case 'Unknown_0004_1200'
        newName = 'RootDirectoryFirstRecord';
    case 'Unknown_0004_1202'
        newName = 'RootDirectoryLastRecord';
    case 'Unknown_0004_1212'
        newName = 'FileSetConsistencyFlag';
    case 'Unknown_0004_1220'
        newName = 'DirectoryRecordSequence';
    case 'Unknown_0004_1400'
        newName = 'NextDirectoryRecordOffset';
    case 'Unknown_0004_1410'
        newName = 'RecordInUseFlag';
    case 'Unknown_0004_1420'
        newName = 'LowerLevelDirectoryOffset';
    case 'Unknown_0004_1430'
        newName = 'DirectoryRecordType';
    case 'Unknown_0004_1432'
        newName = 'PrivateRecordUID';
    case 'Unknown_0004_1500'
        newName = 'ReferencedFileID';
    case 'Unknown_0004_1504'
        newName = 'MRDRDirectoryRecordOffset';
    case 'Unknown_0004_1510'
        newName = 'ReferencedSOPClassUIDInFile';
    case 'Unknown_0004_1511'
        newName = 'ReferencedSOPInstanceUIDInFile';
    case 'Unknown_0004_1512'
        newName = 'ReferencedTransferSyntaxUIDInFile';
    case 'Unknown_0004_151A'
        newName = 'ReferencedRelatedGeneralSOPClassUIDInFile';
    case 'Unknown_0004_1600'
        newName = 'NumberOfReferences';
    case 'Unknown_0008_0000'
        newName = 'IdentifyingGroupLength';
    case 'Unknown_0008_0001'
        newName = 'LengthToEnd';
    case 'Unknown_0008_0005'
        newName = 'SpecificCharacterSet';
    case 'Unknown_0008_0008'
        newName = 'ImageType';
    case 'Unknown_0008_0010'
        newName = 'RecognitionCode';
    case 'Unknown_0008_0012'
        newName = 'InstanceCreationDate';
    case 'Unknown_0008_0013'
        newName = 'InstanceCreationTime';
    case 'Unknown_0008_0014'
        newName = 'InstanceCreatorUID';
    case 'Unknown_0008_0016'
        newName = 'SOPClassUID';
    case 'Unknown_0008_0018'
        newName = 'SOPInstanceUID';
    case 'Unknown_0008_001A'
        newName = 'RelatedGeneralSOPClassUID';
    case 'Unknown_0008_001B'
        newName = 'OriginalSpecializedSOPClassUID';
    case 'Unknown_0008_0020'
        newName = 'StudyDate';
    case 'Unknown_0008_0021'
        newName = 'SeriesDate';
    case 'Unknown_0008_0022'
        newName = 'AcquisitionDate';
    case 'Unknown_0008_0023'
        newName = 'ContentDate';
    case 'Unknown_0008_0024'
        newName = 'OverlayDate';
    case 'Unknown_0008_0025'
        newName = 'CurveDate';
    case 'Unknown_0008_002A'
        newName = 'AcquisitionDateTime';
    case 'Unknown_0008_0030'
        newName = 'StudyTime';
    case 'Unknown_0008_0031'
        newName = 'SeriesTime';
    case 'Unknown_0008_0032'
        newName = 'AcquisitionTime';
    case 'Unknown_0008_0033'
        newName = 'ContentTime';
    case 'Unknown_0008_0034'
        newName = 'OverlayTime';
    case 'Unknown_0008_0035'
        newName = 'CurveTime';
    case 'Unknown_0008_0040'
        newName = 'OldDataSetType';
    case 'Unknown_0008_0041'
        newName = 'OldDataSetSubtype';
    case 'Unknown_0008_0042'
        newName = 'NuclearMedicineSeriesTypeRetired';
    case 'Unknown_0008_0050'
        newName = 'AccessionNumber';
    case 'Unknown_0008_0052'
        newName = 'QueryRetrieveLevel';
    case 'Unknown_0008_0054'
        newName = 'RetrieveAETitle';
    case 'Unknown_0008_0056'
        newName = 'InstanceAvailability';
    case 'Unknown_0008_0058'
        newName = 'FailedSOPInstanceUIDList';
    case 'Unknown_0008_0060'
        newName = 'Modality';
    case 'Unknown_0008_0061'
        newName = 'ModalitiesInStudy';
    case 'Unknown_0008_0062'
        newName = 'SOPClassesInStudy';
    case 'Unknown_0008_0064'
        newName = 'ConversionType';
    case 'Unknown_0008_0068'
        newName = 'PresentationIntentType';
    case 'Unknown_0008_0070'
        newName = 'Manufacturer';
    case 'Unknown_0008_0080'
        newName = 'InstitutionName';
    case 'Unknown_0008_0081'
        newName = 'InstitutionAddress';
    case 'Unknown_0008_0082'
        newName = 'InstitutionCodeSequence';
    case 'Unknown_0008_0090'
        newName = 'ReferringPhysicianName';
    case 'Unknown_0008_0092'
        newName = 'ReferringPhysicianAddress';
    case 'Unknown_0008_0094'
        newName = 'ReferringPhysicianTelephoneNumber';
    case 'Unknown_0008_0096'
        newName = 'ReferringPhysicianIdentificationSequence';
    case 'Unknown_0008_0100'
        newName = 'CodeValue';
    case 'Unknown_0008_0102'
        newName = 'CodingSchemeDesignator';
    case 'Unknown_0008_0103'
        newName = 'CodingSchemeVersion';
    case 'Unknown_0008_0104'
        newName = 'CodeMeaning';
    case 'Unknown_0008_0105'
        newName = 'MappingResource';
    case 'Unknown_0008_0106'
        newName = 'ContextGroupVersion';
    case 'Unknown_0008_0107'
        newName = 'ContextGroupLocalVersion';
    case 'Unknown_0008_010B'
        newName = 'ContextGroupExtensionFlag';
    case 'Unknown_0008_010C'
        newName = 'CodingSchemeUID';
    case 'Unknown_0008_010D'
        newName = 'ContextGroupExtensionCreatorUID';
    case 'Unknown_0008_010F'
        newName = 'ContextIdentifier';
    case 'Unknown_0008_0110'
        newName = 'CodingSchemeIdentificationSequence';
    case 'Unknown_0008_0112'
        newName = 'CodingSchemeRegistry';
    case 'Unknown_0008_0114'
        newName = 'CodingSchemeExternalID';
    case 'Unknown_0008_0115'
        newName = 'CodingSchemeName';
    case 'Unknown_0008_0116'
        newName = 'CodingSchemeResponsibleOrganization';
    case 'Unknown_0008_0201'
        newName = 'TimezoneOffsetFromUTC';
    case 'Unknown_0008_1000'
        newName = 'NetworkID';
    case 'Unknown_0008_1010'
        newName = 'StationName';
    case 'Unknown_0008_1030'
        newName = 'StudyDescription';
    case 'Unknown_0008_1032'
        newName = 'ProcedureCodeSequence';
    case 'Unknown_0008_103E'
        newName = 'SeriesDescription';
    case 'Unknown_0008_1040'
        newName = 'InstitutionalDepartmentName';
    case 'Unknown_0008_1048'
        newName = 'PhysicianOfRecord';
    case 'Unknown_0008_1049'
        newName = 'PhysicianOfRecordIdentificationSequence';
    case 'Unknown_0008_1050'
        newName = 'PerformingPhysicianName';
    case 'Unknown_0008_1052'
        newName = 'PerformingPhysicianIdentificationSequence';
    case 'Unknown_0008_1060'
        newName = 'PhysicianReadingStudy';
    case 'Unknown_0008_1062'
        newName = 'PhysicianReadingStudyIdentificationSequence';
    case 'Unknown_0008_1070'
        newName = 'OperatorName';
    case 'Unknown_0008_1072'
        newName = 'OperatorIdentificationSequence';
    case 'Unknown_0008_1080'
        newName = 'AdmittingDiagnosesDescription';
    case 'Unknown_0008_1084'
        newName = 'AdmittingDiagnosesCodeSequence';
    case 'Unknown_0008_1090'
        newName = 'ManufacturerModelName';
    case 'Unknown_0008_1100'
        newName = 'ReferencedResultsSequence';
    case 'Unknown_0008_1110'
        newName = 'ReferencedStudySequence';
    case 'Unknown_0008_1111'
        newName = 'ReferencedPerformedProcedureStepSequence';
    case 'Unknown_0008_1115'
        newName = 'ReferencedSeriesSequence';
    case 'Unknown_0008_1120'
        newName = 'ReferencedPatientSequence';
    case 'Unknown_0008_1125'
        newName = 'ReferencedVisitSequence';
    case 'Unknown_0008_1130'
        newName = 'ReferencedOverlaySequence';
    case 'Unknown_0008_113A'
        newName = 'ReferencedWaveformSequence';
    case 'Unknown_0008_1140'
        newName = 'ReferencedImageSequence';
    case 'Unknown_0008_1145'
        newName = 'ReferencedCurveSequence';
    case 'Unknown_0008_114A'
        newName = 'ReferencedInstanceSequence';
    case 'Unknown_0008_114B'
        newName = 'ReferencedRealWorldValueMappingInstanceSequence';
    case 'Unknown_0008_1150'
        newName = 'ReferencedSOPClassUID';
    case 'Unknown_0008_1155'
        newName = 'ReferencedSOPInstanceUID';
    case 'Unknown_0008_115A'
        newName = 'SOPClassesSupported';
    case 'Unknown_0008_1160'
        newName = 'ReferencedFrameNumber';
    case 'Unknown_0008_1195'
        newName = 'TransactionUID';
    case 'Unknown_0008_1197'
        newName = 'FailureReason';
    case 'Unknown_0008_1198'
        newName = 'FailedSOPSequence';
    case 'Unknown_0008_1199'
        newName = 'ReferencedSOPSequence';
    case 'Unknown_0008_1200'
        newName = 'StudiesContainingOtherReferencedInstancesSequence';
    case 'Unknown_0008_1250'
        newName = 'RelatedSeriesSequence';
    case 'Unknown_0008_2110'
        newName = 'OldLossyImageCompression';
    case 'Unknown_0008_2111'
        newName = 'DerivationDescription';
    case 'Unknown_0008_2112'
        newName = 'SourceImageSequence';
    case 'Unknown_0008_2120'
        newName = 'StageName';
    case 'Unknown_0008_2122'
        newName = 'StageNumber';
    case 'Unknown_0008_2124'
        newName = 'NumberOfStages';
    case 'Unknown_0008_2127'
        newName = 'ViewName';
    case 'Unknown_0008_2128'
        newName = 'ViewNumber';
    case 'Unknown_0008_2129'
        newName = 'NumberOfEventTimers';
    case 'Unknown_0008_212A'
        newName = 'NumberOfViewsInStage';
    case 'Unknown_0008_2130'
        newName = 'EventElapsedTime';
    case 'Unknown_0008_2132'
        newName = 'EventTimerName';
    case 'Unknown_0008_2142'
        newName = 'StartTrim';
    case 'Unknown_0008_2143'
        newName = 'StopTrim';
    case 'Unknown_0008_2144'
        newName = 'RecommendedDisplayFrameRate';
    case 'Unknown_0008_2200'
        newName = 'TransducerPosition';
    case 'Unknown_0008_2204'
        newName = 'TransducerOrientation';
    case 'Unknown_0008_2208'
        newName = 'AnatomicStructure';
    case 'Unknown_0008_2218'
        newName = 'AnatomicRegionSequence';
    case 'Unknown_0008_2220'
        newName = 'AnatomicRegionModifierSequence';
    case 'Unknown_0008_2228'
        newName = 'PrimaryAnatomicStructureSequence';
    case 'Unknown_0008_2229'
        newName = 'AnatomicStructureSpaceOrRegionSequence';
    case 'Unknown_0008_2230'
        newName = 'PrimaryAnatomicStructureModifierSequence';
    case 'Unknown_0008_2240'
        newName = 'TransducerPositionSequence';
    case 'Unknown_0008_2242'
        newName = 'TransducerPositionModifierSequence';
    case 'Unknown_0008_2244'
        newName = 'TransducerOrientationSequence';
    case 'Unknown_0008_2246'
        newName = 'TransducerOrientationModifierSequence';
    case 'Unknown_0008_2251'
        newName = 'AnatomicStructureSpaceOrRegionCodeSequenceTrial';
    case 'Unknown_0008_2253'
        newName = 'AnatomicPortalOfEntranceCodeSequenceTrial';
    case 'Unknown_0008_2255'
        newName = 'AnatomicApproachDirectionCodeSequenceTrial';
    case 'Unknown_0008_2256'
        newName = 'AnatomicPerspectiveDescriptionTrial';
    case 'Unknown_0008_2257'
        newName = 'AnatomicPerspectiveCodeSequenceTrial';
    case 'Unknown_0008_2258'
        newName = 'AnatomicLocationOfExaminingInstrumentDescriptionTrial';
    case 'Unknown_0008_2259'
        newName = 'AnatomicLocationOfExaminingInstrumentCodeSequenceTrial';
    case 'Unknown_0008_225A'
        newName = 'AnatomicStructureSpaceOrRegionModifierCodeSequenceTrial';
    case 'Unknown_0008_225C'
        newName = 'OnAxisBackgroundAnatomicStructureCodeSequenceTrial';
    case 'Unknown_0008_3001'
        newName = 'AlternateRepresentationSequence';
    case 'Unknown_0008_3010'
        newName = 'IrradiationEventUID';
    case 'Unknown_0008_4000'
        newName = 'IdentifyingComments';
    case 'Unknown_0008_9007'
        newName = 'FrameType';
    case 'Unknown_0008_9092'
        newName = 'ReferencedImageEvidenceSequence';
    case 'Unknown_0008_9121'
        newName = 'ReferencedRawDataSequence';
    case 'Unknown_0008_9123'
        newName = 'CreatorVersionUID';
    case 'Unknown_0008_9124'
        newName = 'DerivationImageSequence';
    case 'Unknown_0008_9154'
        newName = 'SourceImageEvidenceSequence';
    case 'Unknown_0008_9205'
        newName = 'PixelPresentation';
    case 'Unknown_0008_9206'
        newName = 'VolumetricProperties';
    case 'Unknown_0008_9207'
        newName = 'VolumeBasedCalculationTechnique';
    case 'Unknown_0008_9208'
        newName = 'ComplexImageComponent';
    case 'Unknown_0008_9209'
        newName = 'AcquisitionContrast';
    case 'Unknown_0008_9215'
        newName = 'DerivationCodeSequence';
    case 'Unknown_0008_9237'
        newName = 'ReferencedGrayscalePresentationStateSequence';
    case 'Unknown_0008_9410'
        newName = 'ReferencedOtherPlaneSequence';
    case 'Unknown_0008_9458'
        newName = 'FrameDisplaySequence';
    case 'Unknown_0008_9459'
        newName = 'RecommendedDisplayFrameRateInFloat';
    case 'Unknown_0008_9460'
        newName = 'SkipFrameRangeFlag';
    case 'Unknown_0010_0000'
        newName = 'PatientGroupLength';
    case 'Unknown_0010_0010'
        newName = 'PatientName';
    case 'Unknown_0010_0020'
        newName = 'PatientID';
    case 'Unknown_0010_0021'
        newName = 'IssuerOfPatientID';
    case 'Unknown_0010_0022'
        newName = 'TypeOfPatientID';
    case 'Unknown_0010_0030'
        newName = 'PatientBirthDate';
    case 'Unknown_0010_0032'
        newName = 'PatientBirthTime';
    case 'Unknown_0010_0040'
        newName = 'PatientSex';
    case 'Unknown_0010_0050'
        newName = 'PatientInsurancePlanCodeSequence';
    case 'Unknown_0010_0101'
        newName = 'PatientPrimaryLanguageCodeSequence';
    case 'Unknown_0010_0102'
        newName = 'PatientPrimaryLanguageModifierCodeSequence';
    case 'Unknown_0010_1000'
        newName = 'OtherPatientID';
    case 'Unknown_0010_1001'
        newName = 'OtherPatientName';
    case 'Unknown_0010_1002'
        newName = 'OtherPatientIDSequence';
    case 'Unknown_0010_1005'
        newName = 'PatientBirthName';
    case 'Unknown_0010_1010'
        newName = 'PatientAge';
    case 'Unknown_0010_1020'
        newName = 'PatientSize';
    case 'Unknown_0010_1030'
        newName = 'PatientWeight';
    case 'Unknown_0010_1040'
        newName = 'PatientAddress';
    case 'Unknown_0010_1050'
        newName = 'InsurancePlanIdentification';
    case 'Unknown_0010_1060'
        newName = 'PatientMotherBirthName';
    case 'Unknown_0010_1080'
        newName = 'MilitaryRank';
    case 'Unknown_0010_1081'
        newName = 'BranchOfService';
    case 'Unknown_0010_1090'
        newName = 'MedicalRecordLocator';
    case 'Unknown_0010_2000'
        newName = 'MedicalAlerts';
    case 'Unknown_0010_2110'
        newName = 'ContrastAllergies';
    case 'Unknown_0010_2150'
        newName = 'CountryOfResidence';
    case 'Unknown_0010_2152'
        newName = 'RegionOfResidence';
    case 'Unknown_0010_2154'
        newName = 'PatientTelephoneNumber';
    case 'Unknown_0010_2160'
        newName = 'EthnicGroup';
    case 'Unknown_0010_2180'
        newName = 'Occupation';
    case 'Unknown_0010_21A0'
        newName = 'SmokingStatus';
    case 'Unknown_0010_21B0'
        newName = 'AdditionalPatientHistory';
    case 'Unknown_0010_21C0'
        newName = 'PregnancyStatus';
    case 'Unknown_0010_21D0'
        newName = 'LastMenstrualDate';
    case 'Unknown_0010_21F0'
        newName = 'PatientReligiousPreference';
    case 'Unknown_0010_2201'
        newName = 'PatientSpeciesDescription';
    case 'Unknown_0010_2202'
        newName = 'PatientSpeciesCodeSequence';
    case 'Unknown_0010_2203'
        newName = 'PatientSexNeutered';
    case 'Unknown_0010_2292'
        newName = 'PatientBreedDescription';
    case 'Unknown_0010_2293'
        newName = 'PatientBreedCodeSequence';
    case 'Unknown_0010_2294'
        newName = 'BreedRegistrationSequence';
    case 'Unknown_0010_2295'
        newName = 'BreedRegistrationNumber';
    case 'Unknown_0010_2296'
        newName = 'BreedRegistryCodeSequence';
    case 'Unknown_0010_2297'
        newName = 'ResponsiblePerson';
    case 'Unknown_0010_2298'
        newName = 'ResponsiblePersonRole';
    case 'Unknown_0010_2299'
        newName = 'ResponsibleOrganization';
    case 'Unknown_0010_4000'
        newName = 'PatientComments';
    case 'Unknown_0010_9431'
        newName = 'ExaminedBodyThickness';
    case 'Unknown_0012_0000'
        newName = 'ClinicalTrialGroupLength';
    case 'Unknown_0012_0010'
        newName = 'ClinicalTrialSponsorName';
    case 'Unknown_0012_0020'
        newName = 'ClinicalTrialProtocolID';
    case 'Unknown_0012_0021'
        newName = 'ClinicalTrialProtocolName';
    case 'Unknown_0012_0030'
        newName = 'ClinicalTrialSiteID';
    case 'Unknown_0012_0031'
        newName = 'ClinicalTrialSiteName';
    case 'Unknown_0012_0040'
        newName = 'ClinicalTrialSubjectID';
    case 'Unknown_0012_0042'
        newName = 'ClinicalTrialSubjectReadingID';
    case 'Unknown_0012_0050'
        newName = 'ClinicalTrialTimePointID';
    case 'Unknown_0012_0051'
        newName = 'ClinicalTrialTimePointDescription';
    case 'Unknown_0012_0060'
        newName = 'ClinicalTrialCoordinatingCenterName';
    case 'Unknown_0012_0062'
        newName = 'PatientIdentityRemoved';
    case 'Unknown_0012_0063'
        newName = 'DeidentificationMethod';
    case 'Unknown_0012_0064'
        newName = 'DeidentificationMethodCodeSequence';
    case 'Unknown_0018_0000'
        newName = 'AcquisitionGroupLength';
    case 'Unknown_0018_0010'
        newName = 'ContrastBolusAgent';
    case 'Unknown_0018_0012'
        newName = 'ContrastBolusAgentSequence';
    case 'Unknown_0018_0014'
        newName = 'ContrastBolusAdministrationRouteSequence';
    case 'Unknown_0018_0015'
        newName = 'BodyPartExamined';
    case 'Unknown_0018_0020'
        newName = 'ScanningSequence';
    case 'Unknown_0018_0021'
        newName = 'SequenceVariant';
    case 'Unknown_0018_0022'
        newName = 'ScanOptions';
    case 'Unknown_0018_0023'
        newName = 'MRAcquisitionType';
    case 'Unknown_0018_0024'
        newName = 'SequenceName';
    case 'Unknown_0018_0025'
        newName = 'AngioFlag';
    case 'Unknown_0018_0026'
        newName = 'InterventionDrugInformationSequence';
    case 'Unknown_0018_0027'
        newName = 'InterventionDrugStopTime';
    case 'Unknown_0018_0028'
        newName = 'InterventionDrugDose';
    case 'Unknown_0018_0029'
        newName = 'InterventionDrugCodeSequence';
    case 'Unknown_0018_002A'
        newName = 'AdditionalDrugSequence';
    case 'Unknown_0018_0030'
        newName = 'Radionuclide';
    case 'Unknown_0018_0031'
        newName = 'Radiopharmaceutical';
    case 'Unknown_0018_0032'
        newName = 'EnergyWindowCenterline';
    case 'Unknown_0018_0033'
        newName = 'EnergyWindowTotalWidth';
    case 'Unknown_0018_0034'
        newName = 'InterventionDrugName';
    case 'Unknown_0018_0035'
        newName = 'InterventionDrugStartTime';
    case 'Unknown_0018_0036'
        newName = 'InterventionSequence';
    case 'Unknown_0018_0037'
        newName = 'TherapyType';
    case 'Unknown_0018_0038'
        newName = 'InterventionStatus';
    case 'Unknown_0018_0039'
        newName = 'TherapyDescription';
    case 'Unknown_0018_003A'
        newName = 'InterventionDescription';
    case 'Unknown_0018_0040'
        newName = 'CineRate';
    case 'Unknown_0018_0050'
        newName = 'SliceThickness';
    case 'Unknown_0018_0060'
        newName = 'KVP';
    case 'Unknown_0018_0070'
        newName = 'CountsAccumulated';
    case 'Unknown_0018_0071'
        newName = 'AcquisitionTerminationCondition';
    case 'Unknown_0018_0072'
        newName = 'EffectiveDuration';
    case 'Unknown_0018_0073'
        newName = 'AcquisitionStartCondition';
    case 'Unknown_0018_0074'
        newName = 'AcquisitionStartConditionData';
    case 'Unknown_0018_0075'
        newName = 'AcquisitionTerminationConditionData';
    case 'Unknown_0018_0080'
        newName = 'RepetitionTime';
    case 'Unknown_0018_0081'
        newName = 'EchoTime';
    case 'Unknown_0018_0082'
        newName = 'InversionTime';
    case 'Unknown_0018_0083'
        newName = 'NumberOfAverages';
    case 'Unknown_0018_0084'
        newName = 'ImagingFrequency';
    case 'Unknown_0018_0085'
        newName = 'ImagedNucleus';
    case 'Unknown_0018_0086'
        newName = 'EchoNumber';
    case 'Unknown_0018_0087'
        newName = 'MagneticFieldStrength';
    case 'Unknown_0018_0088'
        newName = 'SpacingBetweenSlices';
    case 'Unknown_0018_0089'
        newName = 'NumberOfPhaseEncodingSteps';
    case 'Unknown_0018_0090'
        newName = 'DataCollectionDiameter';
    case 'Unknown_0018_0091'
        newName = 'EchoTrainLength';
    case 'Unknown_0018_0093'
        newName = 'PercentSampling';
    case 'Unknown_0018_0094'
        newName = 'PercentPhaseFieldOfView';
    case 'Unknown_0018_0095'
        newName = 'PixelBandwidth';
    case 'Unknown_0018_1000'
        newName = 'DeviceSerialNumber';
    case 'Unknown_0018_1002'
        newName = 'DeviceUID';
    case 'Unknown_0018_1004'
        newName = 'PlateID';
    case 'Unknown_0018_1005'
        newName = 'GeneratorID';
    case 'Unknown_0018_1006'
        newName = 'GridID';
    case 'Unknown_0018_1007'
        newName = 'CassetteID';
    case 'Unknown_0018_1008'
        newName = 'GantryID';
    case 'Unknown_0018_1010'
        newName = 'SecondaryCaptureDeviceID';
    case 'Unknown_0018_1011'
        newName = 'HardcopyCreationDeviceID';
    case 'Unknown_0018_1012'
        newName = 'DateOfSecondaryCapture';
    case 'Unknown_0018_1014'
        newName = 'TimeOfSecondaryCapture';
    case 'Unknown_0018_1016'
        newName = 'SecondaryCaptureDeviceManufacturer';
    case 'Unknown_0018_1017'
        newName = 'HardcopyDeviceManufacturer';
    case 'Unknown_0018_1018'
        newName = 'SecondaryCaptureDeviceManufacturerModelName';
    case 'Unknown_0018_1019'
        newName = 'SecondaryCaptureDeviceSoftwareVersion';
    case 'Unknown_0018_101A'
        newName = 'HardcopyDeviceSoftwareVersion';
    case 'Unknown_0018_101B'
        newName = 'HardcopyDeviceManufacturerModelName';
    case 'Unknown_0018_1020'
        newName = 'SoftwareVersion';
    case 'Unknown_0018_1022'
        newName = 'VideoImageFormatAcquired';
    case 'Unknown_0018_1023'
        newName = 'DigitalImageFormatAcquired';
    case 'Unknown_0018_1030'
        newName = 'ProtocolName';
    case 'Unknown_0018_1040'
        newName = 'ContrastBolusRoute';
    case 'Unknown_0018_1041'
        newName = 'ContrastBolusVolume';
    case 'Unknown_0018_1042'
        newName = 'ContrastBolusStartTime';
    case 'Unknown_0018_1043'
        newName = 'ContrastBolusStopTime';
    case 'Unknown_0018_1044'
        newName = 'ContrastBolusTotalDose';
    case 'Unknown_0018_1045'
        newName = 'SyringeCounts';
    case 'Unknown_0018_1046'
        newName = 'ContrastFlowRate';
    case 'Unknown_0018_1047'
        newName = 'ContrastFlowDuration';
    case 'Unknown_0018_1048'
        newName = 'ContrastBolusIngredient';
    case 'Unknown_0018_1049'
        newName = 'ContrastBolusIngredientConcentration';
    case 'Unknown_0018_1050'
        newName = 'SpatialResolution';
    case 'Unknown_0018_1060'
        newName = 'TriggerTime';
    case 'Unknown_0018_1061'
        newName = 'TriggerSourceOrType';
    case 'Unknown_0018_1062'
        newName = 'NominalInterval';
    case 'Unknown_0018_1063'
        newName = 'FrameTime';
    case 'Unknown_0018_1064'
        newName = 'FramingType';
    case 'Unknown_0018_1065'
        newName = 'FrameTimeVector';
    case 'Unknown_0018_1066'
        newName = 'FrameDelay';
    case 'Unknown_0018_1067'
        newName = 'ImageTriggerDelay';
    case 'Unknown_0018_1068'
        newName = 'MultiplexGroupTimeOffset';
    case 'Unknown_0018_1069'
        newName = 'TriggerTimeOffset';
    case 'Unknown_0018_106A'
        newName = 'SynchronizationTrigger';
    case 'Unknown_0018_106C'
        newName = 'SynchronizationChannel';
    case 'Unknown_0018_106E'
        newName = 'TriggerSamplePosition';
    case 'Unknown_0018_1070'
        newName = 'RadiopharmaceuticalRoute';
    case 'Unknown_0018_1071'
        newName = 'RadiopharmaceuticalVolume';
    case 'Unknown_0018_1072'
        newName = 'RadiopharmaceuticalStartTime';
    case 'Unknown_0018_1073'
        newName = 'RadiopharmaceuticalStopTime';
    case 'Unknown_0018_1074'
        newName = 'RadionuclideTotalDose';
    case 'Unknown_0018_1075'
        newName = 'RadionuclideHalfLife';
    case 'Unknown_0018_1076'
        newName = 'RadionuclidePositronFraction';
    case 'Unknown_0018_1077'
        newName = 'RadiopharmaceuticalSpecificActivity';
    case 'Unknown_0018_1078'
        newName = 'RadiopharmaceuticalStartDatetime';
    case 'Unknown_0018_1079'
        newName = 'RadiopharmaceuticalStopDatetime';
    case 'Unknown_0018_1080'
        newName = 'BeatRejectionFlag';
    case 'Unknown_0018_1081'
        newName = 'LowRRValue';
    case 'Unknown_0018_1082'
        newName = 'HighRRValue';
    case 'Unknown_0018_1083'
        newName = 'IntervalsAcquired';
    case 'Unknown_0018_1084'
        newName = 'IntervalsRejected';
    case 'Unknown_0018_1085'
        newName = 'PVCRejection';
    case 'Unknown_0018_1086'
        newName = 'SkipBeats';
    case 'Unknown_0018_1088'
        newName = 'HeartRate';
    case 'Unknown_0018_1090'
        newName = 'CardiacNumberOfImages';
    case 'Unknown_0018_1094'
        newName = 'TriggerWindow';
    case 'Unknown_0018_1100'
        newName = 'ReconstructionDiameter';
    case 'Unknown_0018_1110'
        newName = 'DistanceSourceToDetector';
    case 'Unknown_0018_1111'
        newName = 'DistanceSourceToPatient';
    case 'Unknown_0018_1114'
        newName = 'EstimatedRadiographicMagnificationFactor';
    case 'Unknown_0018_1120'
        newName = 'GantryDetectorTilt';
    case 'Unknown_0018_1121'
        newName = 'GantryDetectorSlew';
    case 'Unknown_0018_1130'
        newName = 'TableHeight';
    case 'Unknown_0018_1131'
        newName = 'TableTraverse';
    case 'Unknown_0018_1134'
        newName = 'TableMotion';
    case 'Unknown_0018_1135'
        newName = 'TableVerticalIncrement';
    case 'Unknown_0018_1136'
        newName = 'TableLateralIncrement';
    case 'Unknown_0018_1137'
        newName = 'TableLongitudinalIncrement';
    case 'Unknown_0018_1138'
        newName = 'TableAngle';
    case 'Unknown_0018_113A'
        newName = 'TableType';
    case 'Unknown_0018_1140'
        newName = 'RotationDirection';
    case 'Unknown_0018_1141'
        newName = 'AngularPosition';
    case 'Unknown_0018_1142'
        newName = 'RadialPosition';
    case 'Unknown_0018_1143'
        newName = 'ScanArc';
    case 'Unknown_0018_1144'
        newName = 'AngularStep';
    case 'Unknown_0018_1145'
        newName = 'CenterOfRotationOffset';
    case 'Unknown_0018_1146'
        newName = 'RotationOffset';
    case 'Unknown_0018_1147'
        newName = 'FieldOfViewShape';
    case 'Unknown_0018_1149'
        newName = 'FieldOfViewDimensions';
    case 'Unknown_0018_1150'
        newName = 'ExposureTime';
    case 'Unknown_0018_1151'
        newName = 'XrayTubeCurrent';
    case 'Unknown_0018_1152'
        newName = 'Exposure';
    case 'Unknown_0018_1153'
        newName = 'ExposureInuAs';
    case 'Unknown_0018_1154'
        newName = 'AveragePulseWidth';
    case 'Unknown_0018_1155'
        newName = 'RadiationSetting';
    case 'Unknown_0018_1156'
        newName = 'RectificationType';
    case 'Unknown_0018_115A'
        newName = 'RadiationMode';
    case 'Unknown_0018_115E'
        newName = 'ImageAndFluoroscopyAreaDoseProduct';
    case 'Unknown_0018_1160'
        newName = 'FilterType';
    case 'Unknown_0018_1161'
        newName = 'TypeOfFilters';
    case 'Unknown_0018_1162'
        newName = 'IntensifierSize';
    case 'Unknown_0018_1164'
        newName = 'ImagerPixelSpacing';
    case 'Unknown_0018_1166'
        newName = 'Grid';
    case 'Unknown_0018_1170'
        newName = 'GeneratorPower';
    case 'Unknown_0018_1180'
        newName = 'CollimatorGridName';
    case 'Unknown_0018_1181'
        newName = 'CollimatorType';
    case 'Unknown_0018_1182'
        newName = 'FocalDistance';
    case 'Unknown_0018_1183'
        newName = 'XFocusCenter';
    case 'Unknown_0018_1184'
        newName = 'YFocusCenter';
    case 'Unknown_0018_1190'
        newName = 'FocalSpot';
    case 'Unknown_0018_1191'
        newName = 'AnodeTargetMaterial';
    case 'Unknown_0018_11A0'
        newName = 'BodyPartThickness';
    case 'Unknown_0018_11A2'
        newName = 'CompressionForce';
    case 'Unknown_0018_1200'
        newName = 'DateOfLastCalibration';
    case 'Unknown_0018_1201'
        newName = 'TimeOfLastCalibration';
    case 'Unknown_0018_1210'
        newName = 'ConvolutionKernel';
    case 'Unknown_0018_1240'
        newName = 'UpperLowerPixelValues';
    case 'Unknown_0018_1242'
        newName = 'ActualFrameDuration';
    case 'Unknown_0018_1243'
        newName = 'CountRate';
    case 'Unknown_0018_1244'
        newName = 'PreferredPlaybackSequencing';
    case 'Unknown_0018_1250'
        newName = 'ReceiveCoilName';
    case 'Unknown_0018_1251'
        newName = 'TransmitCoilName';
    case 'Unknown_0018_1260'
        newName = 'PlateType';
    case 'Unknown_0018_1261'
        newName = 'PhosphorType';
    case 'Unknown_0018_1300'
        newName = 'ScanVelocity';
    case 'Unknown_0018_1301'
        newName = 'WholeBodyTechnique';
    case 'Unknown_0018_1302'
        newName = 'ScanLength';
    case 'Unknown_0018_1310'
        newName = 'AcquisitionMatrix';
    case 'Unknown_0018_1312'
        newName = 'InPlanePhaseEncodingDirection';
    case 'Unknown_0018_1314'
        newName = 'FlipAngle';
    case 'Unknown_0018_1315'
        newName = 'VariableFlipAngleFlag';
    case 'Unknown_0018_1316'
        newName = 'SAR';
    case 'Unknown_0018_1318'
        newName = 'dBdt';
    case 'Unknown_0018_1400'
        newName = 'AcquisitionDeviceProcessingDescription';
    case 'Unknown_0018_1401'
        newName = 'AcquisitionDeviceProcessingCode';
    case 'Unknown_0018_1402'
        newName = 'CassetteOrientation';
    case 'Unknown_0018_1403'
        newName = 'CassetteSize';
    case 'Unknown_0018_1404'
        newName = 'ExposuresOnPlate';
    case 'Unknown_0018_1405'
        newName = 'RelativeXrayExposure';
    case 'Unknown_0018_1450'
        newName = 'ColumnAngulation';
    case 'Unknown_0018_1460'
        newName = 'TomoLayerHeight';
    case 'Unknown_0018_1470'
        newName = 'TomoAngle';
    case 'Unknown_0018_1480'
        newName = 'TomoTime';
    case 'Unknown_0018_1490'
        newName = 'TomoType';
    case 'Unknown_0018_1491'
        newName = 'TomoClass';
    case 'Unknown_0018_1495'
        newName = 'NumberOfTomosynthesisSourceImages';
    case 'Unknown_0018_1500'
        newName = 'PositionerMotion';
    case 'Unknown_0018_1508'
        newName = 'PositionerType';
    case 'Unknown_0018_1510'
        newName = 'PositionerPrimaryAngle';
    case 'Unknown_0018_1511'
        newName = 'PositionerSecondaryAngle';
    case 'Unknown_0018_1520'
        newName = 'PositionerPrimaryAngleIncrement';
    case 'Unknown_0018_1521'
        newName = 'PositionerSecondaryAngleIncrement';
    case 'Unknown_0018_1530'
        newName = 'DetectorPrimaryAngle';
    case 'Unknown_0018_1531'
        newName = 'DetectorSecondaryAngle';
    case 'Unknown_0018_1600'
        newName = 'ShutterShape';
    case 'Unknown_0018_1602'
        newName = 'ShutterLeftVerticalEdge';
    case 'Unknown_0018_1604'
        newName = 'ShutterRightVerticalEdge';
    case 'Unknown_0018_1606'
        newName = 'ShutterUpperHorizontalEdge';
    case 'Unknown_0018_1608'
        newName = 'ShutterLowerHorizontalEdge';
    case 'Unknown_0018_1610'
        newName = 'CenterOfCircularShutter';
    case 'Unknown_0018_1612'
        newName = 'RadiusOfCircularShutter';
    case 'Unknown_0018_1620'
        newName = 'VerticesOfPolygonalShutter';
    case 'Unknown_0018_1622'
        newName = 'ShutterPresentationValue';
    case 'Unknown_0018_1623'
        newName = 'ShutterOverlayGroup';
    case 'Unknown_0018_1624'
        newName = 'ShutterPresentationColorCIELabValue';
    case 'Unknown_0018_1700'
        newName = 'CollimatorShape';
    case 'Unknown_0018_1702'
        newName = 'CollimatorLeftVerticalEdge';
    case 'Unknown_0018_1704'
        newName = 'CollimatorRightVerticalEdge';
    case 'Unknown_0018_1706'
        newName = 'CollimatorUpperHorizontalEdge';
    case 'Unknown_0018_1708'
        newName = 'CollimatorLowerHorizontalEdge';
    case 'Unknown_0018_1710'
        newName = 'CenterOfCircularCollimator';
    case 'Unknown_0018_1712'
        newName = 'RadiusOfCircularCollimator';
    case 'Unknown_0018_1720'
        newName = 'VerticesOfPolygonalCollimator';
    case 'Unknown_0018_1800'
        newName = 'AcquisitionTimeSynchronized';
    case 'Unknown_0018_1801'
        newName = 'TimeSource';
    case 'Unknown_0018_1802'
        newName = 'TimeDistributionProtocol';
    case 'Unknown_0018_1803'
        newName = 'NTPSourceAddress';
    case 'Unknown_0018_2001'
        newName = 'PageNumberVector';
    case 'Unknown_0018_2002'
        newName = 'FrameLabelVector';
    case 'Unknown_0018_2003'
        newName = 'FramePrimaryAngleVector';
    case 'Unknown_0018_2004'
        newName = 'FrameSecondaryAngleVector';
    case 'Unknown_0018_2005'
        newName = 'SliceLocationVector';
    case 'Unknown_0018_2006'
        newName = 'DisplayWindowLabelVector';
    case 'Unknown_0018_2010'
        newName = 'NominalScannedPixelSpacing';
    case 'Unknown_0018_2020'
        newName = 'DigitizingDeviceTransportDirection';
    case 'Unknown_0018_2030'
        newName = 'RotationOfScannedFilm';
    case 'Unknown_0018_3100'
        newName = 'IVUSAcquisition';
    case 'Unknown_0018_3101'
        newName = 'IVUSPullbackRate';
    case 'Unknown_0018_3102'
        newName = 'IVUSGatedRate';
    case 'Unknown_0018_3103'
        newName = 'IVUSPullbackStartFrameNumber';
    case 'Unknown_0018_3104'
        newName = 'IVUSPullbackStopFrameNumber';
    case 'Unknown_0018_3105'
        newName = 'LesionNumber';
    case 'Unknown_0018_4000'
        newName = 'AcquisitionComments';
    case 'Unknown_0018_5000'
        newName = 'OutputPower';
    case 'Unknown_0018_5010'
        newName = 'TransducerData';
    case 'Unknown_0018_5012'
        newName = 'FocusDepth';
    case 'Unknown_0018_5020'
        newName = 'ProcessingFunction';
    case 'Unknown_0018_5021'
        newName = 'PostprocessingFunction';
    case 'Unknown_0018_5022'
        newName = 'MechanicalIndex';
    case 'Unknown_0018_5024'
        newName = 'BoneThermalIndex';
    case 'Unknown_0018_5026'
        newName = 'CranialThermalIndex';
    case 'Unknown_0018_5027'
        newName = 'SoftTissueThermalIndex';
    case 'Unknown_0018_5028'
        newName = 'SoftTissueFocusThermalIndex';
    case 'Unknown_0018_5029'
        newName = 'SoftTissueSurfaceThermalIndex';
    case 'Unknown_0018_5030'
        newName = 'DynamicRange';
    case 'Unknown_0018_5040'
        newName = 'TotalGain';
    case 'Unknown_0018_5050'
        newName = 'DepthOfScanField';
    case 'Unknown_0018_5100'
        newName = 'PatientPosition';
    case 'Unknown_0018_5101'
        newName = 'ViewPosition';
    case 'Unknown_0018_5104'
        newName = 'ProjectionEponymousNameCodeSequence';
    case 'Unknown_0018_5210'
        newName = 'ImageTransformationMatrix';
    case 'Unknown_0018_5212'
        newName = 'ImageTranslationVector';
    case 'Unknown_0018_6000'
        newName = 'Sensitivity';
    case 'Unknown_0018_6011'
        newName = 'SequenceOfUltrasoundRegions';
    case 'Unknown_0018_6012'
        newName = 'RegionSpatialFormat';
    case 'Unknown_0018_6014'
        newName = 'RegionDataType';
    case 'Unknown_0018_6016'
        newName = 'RegionFlags';
    case 'Unknown_0018_6018'
        newName = 'RegionLocationMinX0';
    case 'Unknown_0018_601A'
        newName = 'RegionLocationMinY0';
    case 'Unknown_0018_601C'
        newName = 'RegionLocationMaxX1';
    case 'Unknown_0018_601E'
        newName = 'RegionLocationMaxY1';
    case 'Unknown_0018_6020'
        newName = 'ReferencePixelX0';
    case 'Unknown_0018_6022'
        newName = 'ReferencePixelY0';
    case 'Unknown_0018_6024'
        newName = 'PhysicalUnitsXDirection';
    case 'Unknown_0018_6026'
        newName = 'PhysicalUnitsYDirection';
    case 'Unknown_0018_6028'
        newName = 'ReferencePixelPhysicalValueX';
    case 'Unknown_0018_602A'
        newName = 'ReferencePixelPhysicalValueY';
    case 'Unknown_0018_602C'
        newName = 'PhysicalDeltaX';
    case 'Unknown_0018_602E'
        newName = 'PhysicalDeltaY';
    case 'Unknown_0018_6030'
        newName = 'TransducerFrequency';
    case 'Unknown_0018_6031'
        newName = 'TransducerType';
    case 'Unknown_0018_6032'
        newName = 'PulseRepetitionFrequency';
    case 'Unknown_0018_6034'
        newName = 'DopplerCorrectionAngle';
    case 'Unknown_0018_6036'
        newName = 'SteeringAngle';
    case 'Unknown_0018_6038'
        newName = 'DopplerSampleVolumeXPositionRetired';
    case 'Unknown_0018_6039'
        newName = 'DopplerSampleVolumeXPosition';
    case 'Unknown_0018_603A'
        newName = 'DopplerSampleVolumeYPositionRetired';
    case 'Unknown_0018_603B'
        newName = 'DopplerSampleVolumeYPosition';
    case 'Unknown_0018_603C'
        newName = 'TMLinePositionX0Retired';
    case 'Unknown_0018_603D'
        newName = 'TMLinePositionX0';
    case 'Unknown_0018_603E'
        newName = 'TMLinePositionY0Retired';
    case 'Unknown_0018_603F'
        newName = 'TMLinePositionY0';
    case 'Unknown_0018_6040'
        newName = 'TMLinePositionX1Retired';
    case 'Unknown_0018_6041'
        newName = 'TMLinePositionX1';
    case 'Unknown_0018_6042'
        newName = 'TMLinePositionY1Retired';
    case 'Unknown_0018_6043'
        newName = 'TMLinePositionY1';
    case 'Unknown_0018_6044'
        newName = 'PixelComponentOrganization';
    case 'Unknown_0018_6046'
        newName = 'PixelComponentMask';
    case 'Unknown_0018_6048'
        newName = 'PixelComponentRangeStart';
    case 'Unknown_0018_604A'
        newName = 'PixelComponentRangeStop';
    case 'Unknown_0018_604C'
        newName = 'PixelComponentPhysicalUnits';
    case 'Unknown_0018_604E'
        newName = 'PixelComponentDataType';
    case 'Unknown_0018_6050'
        newName = 'NumberOfTableBreakPoints';
    case 'Unknown_0018_6052'
        newName = 'TableOfXBreakPoints';
    case 'Unknown_0018_6054'
        newName = 'TableOfYBreakPoints';
    case 'Unknown_0018_6056'
        newName = 'NumberOfTableEntries';
    case 'Unknown_0018_6058'
        newName = 'TableOfPixelValues';
    case 'Unknown_0018_605A'
        newName = 'TableOfParameterValues';
    case 'Unknown_0018_6060'
        newName = 'RWaveTimeVector';
    case 'Unknown_0018_7000'
        newName = 'DetectorConditionsNominalFlag';
    case 'Unknown_0018_7001'
        newName = 'DetectorTemperature';
    case 'Unknown_0018_7004'
        newName = 'DetectorType';
    case 'Unknown_0018_7005'
        newName = 'DetectorConfiguration';
    case 'Unknown_0018_7006'
        newName = 'DetectorDescription';
    case 'Unknown_0018_7008'
        newName = 'DetectorMode';
    case 'Unknown_0018_700A'
        newName = 'DetectorID';
    case 'Unknown_0018_700C'
        newName = 'DateOfLastDetectorCalibration';
    case 'Unknown_0018_700E'
        newName = 'TimeOfLastDetectorCalibration';
    case 'Unknown_0018_7010'
        newName = 'ExposuresOnDetectorSinceLastCalibration';
    case 'Unknown_0018_7011'
        newName = 'ExposuresOnDetectorSinceManufactured';
    case 'Unknown_0018_7012'
        newName = 'DetectorTimeSinceLastExposure';
    case 'Unknown_0018_7014'
        newName = 'DetectorActiveTime';
    case 'Unknown_0018_7016'
        newName = 'DetectorActivationOffsetFromExposure';
    case 'Unknown_0018_701A'
        newName = 'DetectorBinning';
    case 'Unknown_0018_7020'
        newName = 'DetectorElementPhysicalSize';
    case 'Unknown_0018_7022'
        newName = 'DetectorElementSpacing';
    case 'Unknown_0018_7024'
        newName = 'DetectorActiveShape';
    case 'Unknown_0018_7026'
        newName = 'DetectorActiveDimensions';
    case 'Unknown_0018_7028'
        newName = 'DetectorActiveOrigin';
    case 'Unknown_0018_702A'
        newName = 'DetectorManufacturerName';
    case 'Unknown_0018_702B'
        newName = 'DetectorManufacturerModelName';
    case 'Unknown_0018_7030'
        newName = 'FieldOfViewOrigin';
    case 'Unknown_0018_7032'
        newName = 'FieldOfViewRotation';
    case 'Unknown_0018_7034'
        newName = 'FieldOfViewHorizontalFlip';
    case 'Unknown_0018_7040'
        newName = 'GridAbsorbingMaterial';
    case 'Unknown_0018_7041'
        newName = 'GridSpacingMaterial';
    case 'Unknown_0018_7042'
        newName = 'GridThickness';
    case 'Unknown_0018_7044'
        newName = 'GridPitch';
    case 'Unknown_0018_7046'
        newName = 'GridAspectRatio';
    case 'Unknown_0018_7048'
        newName = 'GridPeriod';
    case 'Unknown_0018_704C'
        newName = 'GridFocalDistance';
    case 'Unknown_0018_7050'
        newName = 'FilterMaterial';
    case 'Unknown_0018_7052'
        newName = 'FilterThicknessMinimum';
    case 'Unknown_0018_7054'
        newName = 'FilterThicknessMaximum';
    case 'Unknown_0018_7060'
        newName = 'ExposureControlMode';
    case 'Unknown_0018_7062'
        newName = 'ExposureControlModeDescription';
    case 'Unknown_0018_7064'
        newName = 'ExposureStatus';
    case 'Unknown_0018_7065'
        newName = 'PhototimerSetting';
    case 'Unknown_0018_8150'
        newName = 'ExposureTimeInuS';
    case 'Unknown_0018_8151'
        newName = 'XrayTubeCurrentInuA';
    case 'Unknown_0018_9004'
        newName = 'ContentQualification';
    case 'Unknown_0018_9005'
        newName = 'PulseSequenceName';
    case 'Unknown_0018_9006'
        newName = 'MRImagingModifierSequence';
    case 'Unknown_0018_9008'
        newName = 'EchoPulseSequence';
    case 'Unknown_0018_9009'
        newName = 'InversionRecovery';
    case 'Unknown_0018_9010'
        newName = 'FlowCompensation';
    case 'Unknown_0018_9011'
        newName = 'MultipleSpinEcho';
    case 'Unknown_0018_9012'
        newName = 'MultiplanarExcitation';
    case 'Unknown_0018_9014'
        newName = 'PhaseContrast';
    case 'Unknown_0018_9015'
        newName = 'TimeOfFlightContrast';
    case 'Unknown_0018_9016'
        newName = 'Spoiling';
    case 'Unknown_0018_9017'
        newName = 'SteadyStatePulseSequence';
    case 'Unknown_0018_9018'
        newName = 'EchoPlanarPulseSequence';
    case 'Unknown_0018_9019'
        newName = 'TagAngleFirstAxis';
    case 'Unknown_0018_9020'
        newName = 'MagnetizationTransfer';
    case 'Unknown_0018_9021'
        newName = 'T2Preparation';
    case 'Unknown_0018_9022'
        newName = 'BloodSignalNulling';
    case 'Unknown_0018_9024'
        newName = 'SaturationRecovery';
    case 'Unknown_0018_9025'
        newName = 'SpectrallySelectedSuppression';
    case 'Unknown_0018_9026'
        newName = 'SpectrallySelectedExcitation';
    case 'Unknown_0018_9027'
        newName = 'SpatialPresaturation';
    case 'Unknown_0018_9028'
        newName = 'Tagging';
    case 'Unknown_0018_9029'
        newName = 'OversamplingPhase';
    case 'Unknown_0018_9030'
        newName = 'TagSpacingFirstDimension';
    case 'Unknown_0018_9032'
        newName = 'GeometryOfKSpaceTraversal';
    case 'Unknown_0018_9033'
        newName = 'SegmentedKSpaceTraversal';
    case 'Unknown_0018_9034'
        newName = 'RectilinearPhaseEncodeReordering';
    case 'Unknown_0018_9035'
        newName = 'TagThickness';
    case 'Unknown_0018_9036'
        newName = 'PartialFourierDirection';
    case 'Unknown_0018_9037'
        newName = 'CardiacSynchronizationTechnique';
    case 'Unknown_0018_9041'
        newName = 'ReceiveCoilManufacturerName';
    case 'Unknown_0018_9042'
        newName = 'MRReceiveCoilSequence';
    case 'Unknown_0018_9043'
        newName = 'ReceiveCoilType';
    case 'Unknown_0018_9044'
        newName = 'QuadratureReceiveCoil';
    case 'Unknown_0018_9045'
        newName = 'MultiCoilDefinitionSequence';
    case 'Unknown_0018_9046'
        newName = 'MultiCoilConfiguration';
    case 'Unknown_0018_9047'
        newName = 'MultiCoilElementName';
    case 'Unknown_0018_9048'
        newName = 'MultiCoilElementUsed';
    case 'Unknown_0018_9049'
        newName = 'MRTransmitCoilSequence';
    case 'Unknown_0018_9050'
        newName = 'TransmitCoilManufacturerName';
    case 'Unknown_0018_9051'
        newName = 'TransmitCoilType';
    case 'Unknown_0018_9052'
        newName = 'SpectralWidth';
    case 'Unknown_0018_9053'
        newName = 'ChemicalShiftReference';
    case 'Unknown_0018_9054'
        newName = 'VolumeLocalizationTechnique';
    case 'Unknown_0018_9058'
        newName = 'MRAcquisitionFrequencyEncodingSteps';
    case 'Unknown_0018_9059'
        newName = 'Decoupling';
    case 'Unknown_0018_9060'
        newName = 'DecoupledNucleus';
    case 'Unknown_0018_9061'
        newName = 'DecouplingFrequency';
    case 'Unknown_0018_9062'
        newName = 'DecouplingMethod';
    case 'Unknown_0018_9063'
        newName = 'DecouplingChemicalShiftReference';
    case 'Unknown_0018_9064'
        newName = 'KSpaceFiltering';
    case 'Unknown_0018_9065'
        newName = 'TimeDomainFiltering';
    case 'Unknown_0018_9066'
        newName = 'NumberOfZeroFills';
    case 'Unknown_0018_9067'
        newName = 'BaselineCorrection';
    case 'Unknown_0018_9069'
        newName = 'ParallelReductionFactorInPlane';
    case 'Unknown_0018_9070'
        newName = 'CardiacRRIntervalSpecified';
    case 'Unknown_0018_9073'
        newName = 'AcquisitionDuration';
    case 'Unknown_0018_9074'
        newName = 'FrameAcquisitionDatetime';
    case 'Unknown_0018_9075'
        newName = 'DiffusionDirectionality';
    case 'Unknown_0018_9076'
        newName = 'DiffusionGradientDirectionSequence';
    case 'Unknown_0018_9077'
        newName = 'ParallelAcquisition';
    case 'Unknown_0018_9078'
        newName = 'ParallelAcquisitionTechnique';
    case 'Unknown_0018_9079'
        newName = 'InversionTimes';
    case 'Unknown_0018_9080'
        newName = 'MetaboliteMapDescription';
    case 'Unknown_0018_9081'
        newName = 'PartialFourier';
    case 'Unknown_0018_9082'
        newName = 'EffectiveEchoTime';
    case 'Unknown_0018_9083'
        newName = 'MetaboliteCodeSequence';
    case 'Unknown_0018_9084'
        newName = 'ChemicalShiftSequence';
    case 'Unknown_0018_9085'
        newName = 'CardiacSignalSource';
    case 'Unknown_0018_9087'
        newName = 'DiffusionBValue';
    case 'Unknown_0018_9089'
        newName = 'DiffusionGradientOrientation';
    case 'Unknown_0018_9090'
        newName = 'VelocityEncodingDirection';
    case 'Unknown_0018_9091'
        newName = 'VelocityEncodingMinimumValue';
    case 'Unknown_0018_9093'
        newName = 'NumberOfKSpaceTrajectories';
    case 'Unknown_0018_9094'
        newName = 'CoverageOfKSpace';
    case 'Unknown_0018_9095'
        newName = 'SpectroscopyAcquisitionPhaseRows';
    case 'Unknown_0018_9098'
        newName = 'TransmitterFrequency';
    case 'Unknown_0018_9100'
        newName = 'ResonantNucleus';
    case 'Unknown_0018_9101'
        newName = 'FrequencyCorrection';
    case 'Unknown_0018_9103'
        newName = 'MRSpectroscopyFOVGeometrySequence';
    case 'Unknown_0018_9104'
        newName = 'SlabThickness';
    case 'Unknown_0018_9105'
        newName = 'SlabOrientation';
    case 'Unknown_0018_9106'
        newName = 'MidSlabPosition';
    case 'Unknown_0018_9107'
        newName = 'MRSpatialSaturationSequence';
    case 'Unknown_0018_9112'
        newName = 'MRTimingAndRelatedParametersSequence';
    case 'Unknown_0018_9114'
        newName = 'MREchoSequence';
    case 'Unknown_0018_9115'
        newName = 'MRModifierSequence';
    case 'Unknown_0018_9117'
        newName = 'MRDiffusionSequence';
    case 'Unknown_0018_9118'
        newName = 'CardiacTriggerSequence';
    case 'Unknown_0018_9119'
        newName = 'MRAveragesSequence';
    case 'Unknown_0018_9125'
        newName = 'MRFOVGeometrySequence';
    case 'Unknown_0018_9126'
        newName = 'VolumeLocalizationSequence';
    case 'Unknown_0018_9127'
        newName = 'SpectroscopyAcquisitionDataColumns';
    case 'Unknown_0018_9147'
        newName = 'DiffusionAnisotropyType';
    case 'Unknown_0018_9151'
        newName = 'FrameReferenceDatetime';
    case 'Unknown_0018_9152'
        newName = 'MRMetaboliteMapSequence';
    case 'Unknown_0018_9155'
        newName = 'ParallelReductionFactorOutOfPlane';
    case 'Unknown_0018_9159'
        newName = 'SpectroscopyAcquisitionOutOfPlanePhaseSteps';
    case 'Unknown_0018_9166'
        newName = 'BulkMotionStatus';
    case 'Unknown_0018_9168'
        newName = 'ParallelReductionFactorSecondInPlane';
    case 'Unknown_0018_9169'
        newName = 'CardiacBeatRejectionTechnique';
    case 'Unknown_0018_9170'
        newName = 'RespiratoryMotionCompensationTechnique';
    case 'Unknown_0018_9171'
        newName = 'RespiratorySignalSource';
    case 'Unknown_0018_9172'
        newName = 'BulkMotionCompensationTechnique';
    case 'Unknown_0018_9173'
        newName = 'BulkMotionSignalSource';
    case 'Unknown_0018_9174'
        newName = 'ApplicableSafetyStandardAgency';
    case 'Unknown_0018_9175'
        newName = 'ApplicableSafetyStandardDescription';
    case 'Unknown_0018_9176'
        newName = 'OperatingModeSequence';
    case 'Unknown_0018_9177'
        newName = 'OperatingModeType';
    case 'Unknown_0018_9178'
        newName = 'OperatingMode';
    case 'Unknown_0018_9179'
        newName = 'SpecificAbsorptionRateDefinition';
    case 'Unknown_0018_9180'
        newName = 'GradientOutputType';
    case 'Unknown_0018_9181'
        newName = 'SpecificAbsorptionRateValue';
    case 'Unknown_0018_9182'
        newName = 'GradientOutput';
    case 'Unknown_0018_9183'
        newName = 'FlowCompensationDirection';
    case 'Unknown_0018_9184'
        newName = 'TaggingDelay';
    case 'Unknown_0018_9185'
        newName = 'RespiratoryMotionCompensationTechniqueDescription';
    case 'Unknown_0018_9186'
        newName = 'RespiratorySignalSourceID';
    case 'Unknown_0018_9195'
        newName = 'ChemicalShiftMinimumIntegrationLimitInHz';
    case 'Unknown_0018_9196'
        newName = 'ChemicalShiftMaximumIntegrationLimitInHz';
    case 'Unknown_0018_9197'
        newName = 'MRVelocityEncodingSequence';
    case 'Unknown_0018_9198'
        newName = 'FirstOrderPhaseCorrection';
    case 'Unknown_0018_9199'
        newName = 'WaterReferencedPhaseCorrection';
    case 'Unknown_0018_9200'
        newName = 'MRSpectroscopyAcquisitionType';
    case 'Unknown_0018_9214'
        newName = 'RespiratoryCyclePosition';
    case 'Unknown_0018_9217'
        newName = 'VelocityEncodingMaximumValue';
    case 'Unknown_0018_9218'
        newName = 'TagSpacingSecondDimension';
    case 'Unknown_0018_9219'
        newName = 'TagAngleSecondAxis';
    case 'Unknown_0018_9220'
        newName = 'FrameAcquisitionDuration';
    case 'Unknown_0018_9226'
        newName = 'MRImageFrameTypeSequence';
    case 'Unknown_0018_9227'
        newName = 'MRSpectroscopyFrameTypeSequence';
    case 'Unknown_0018_9231'
        newName = 'MRAcquisitionPhaseEncodingStepsInPlane';
    case 'Unknown_0018_9232'
        newName = 'MRAcquisitionPhaseEncodingStepsOutOfPlane';
    case 'Unknown_0018_9234'
        newName = 'SpectroscopyAcquisitionPhaseColumns';
    case 'Unknown_0018_9236'
        newName = 'CardiacCyclePosition';
    case 'Unknown_0018_9239'
        newName = 'SpecificAbsorptionRateSequence';
    case 'Unknown_0018_9240'
        newName = 'RFEchoTrainLength';
    case 'Unknown_0018_9241'
        newName = 'GradientEchoTrainLength';
    case 'Unknown_0018_9295'
        newName = 'ChemicalShiftMinimumIntegrationLimitInPPM';
    case 'Unknown_0018_9296'
        newName = 'ChemicalShiftMaximumIntegrationLimitInPPM';
    case 'Unknown_0018_9301'
        newName = 'CTAcquisitionTypeSequence';
    case 'Unknown_0018_9302'
        newName = 'AcquisitionType';
    case 'Unknown_0018_9303'
        newName = 'TubeAngle';
    case 'Unknown_0018_9304'
        newName = 'CTAcquisitionDetailsSequence';
    case 'Unknown_0018_9305'
        newName = 'RevolutionTime';
    case 'Unknown_0018_9306'
        newName = 'SingleCollimationWidth';
    case 'Unknown_0018_9307'
        newName = 'TotalCollimationWidth';
    case 'Unknown_0018_9308'
        newName = 'CTTableDynamicsSequence';
    case 'Unknown_0018_9309'
        newName = 'TableSpeed';
    case 'Unknown_0018_9310'
        newName = 'TableFeedPerRotation';
    case 'Unknown_0018_9311'
        newName = 'SpiralPitchFactor';
    case 'Unknown_0018_9312'
        newName = 'CTGeometrySequence';
    case 'Unknown_0018_9313'
        newName = 'DataCollectionCenterPatient';
    case 'Unknown_0018_9314'
        newName = 'CTReconstructionSequence';
    case 'Unknown_0018_9315'
        newName = 'ReconstructionAlgorithm';
    case 'Unknown_0018_9316'
        newName = 'ConvolutionKernelGroup';
    case 'Unknown_0018_9317'
        newName = 'ReconstructionFieldOfView';
    case 'Unknown_0018_9318'
        newName = 'ReconstructionTargetCenterPatient';
    case 'Unknown_0018_9319'
        newName = 'ReconstructionAngle';
    case 'Unknown_0018_9320'
        newName = 'ImageFilter';
    case 'Unknown_0018_9321'
        newName = 'CTExposureSequence';
    case 'Unknown_0018_9322'
        newName = 'ReconstructionPixelSpacing';
    case 'Unknown_0018_9323'
        newName = 'ExposureModulationType';
    case 'Unknown_0018_9324'
        newName = 'EstimatedDoseSaving';
    case 'Unknown_0018_9325'
        newName = 'CTXrayDetailsSequence';
    case 'Unknown_0018_9326'
        newName = 'CTPositionSequence';
    case 'Unknown_0018_9327'
        newName = 'TablePosition';
    case 'Unknown_0018_9328'
        newName = 'ExposureTimeInms';
    case 'Unknown_0018_9329'
        newName = 'CTImageFrameTypeSequence';
    case 'Unknown_0018_9330'
        newName = 'XrayTubeCurrentInmA';
    case 'Unknown_0018_9332'
        newName = 'ExposureInmAs';
    case 'Unknown_0018_9333'
        newName = 'ConstantVolumeFlag';
    case 'Unknown_0018_9334'
        newName = 'FluoroscopyFlag';
    case 'Unknown_0018_9335'
        newName = 'DistanceSourceToDataCollectionCenter';
    case 'Unknown_0018_9337'
        newName = 'ContrastBolusAgentNumber';
    case 'Unknown_0018_9338'
        newName = 'ContrastBolusIngredientCodeSequence';
    case 'Unknown_0018_9340'
        newName = 'ContrastAdministrationProfileSequence';
    case 'Unknown_0018_9341'
        newName = 'ContrastBolusUsageSequence';
    case 'Unknown_0018_9342'
        newName = 'ContrastBolusAgentAdministered';
    case 'Unknown_0018_9343'
        newName = 'ContrastBolusAgentDetected';
    case 'Unknown_0018_9344'
        newName = 'ContrastBolusAgentPhase';
    case 'Unknown_0018_9345'
        newName = 'CTDIvol';
    case 'Unknown_0018_9401'
        newName = 'ProjectionPixelCalibrationSequence';
    case 'Unknown_0018_9402'
        newName = 'DistanceSourceToIsocenter';
    case 'Unknown_0018_9403'
        newName = 'DistanceObjectToTableTop';
    case 'Unknown_0018_9404'
        newName = 'ObjectPixelSpacingInCenterOfBeam';
    case 'Unknown_0018_9405'
        newName = 'PositionerPositionSequence';
    case 'Unknown_0018_9406'
        newName = 'TablePositionSequence';
    case 'Unknown_0018_9407'
        newName = 'CollimatorShapeSequence';
    case 'Unknown_0018_9412'
        newName = 'XAXRFFrameCharacteristicsSequence';
    case 'Unknown_0018_9417'
        newName = 'FrameAcquisitionSequence';
    case 'Unknown_0018_9420'
        newName = 'XRayReceptorType';
    case 'Unknown_0018_9423'
        newName = 'AcquisitionProtocolName';
    case 'Unknown_0018_9424'
        newName = 'AcquisitionProtocolDescription';
    case 'Unknown_0018_9425'
        newName = 'ContrastBolusIngredientOpaque';
    case 'Unknown_0018_9426'
        newName = 'DistanceReceptorPlaneToDetectorHousing';
    case 'Unknown_0018_9427'
        newName = 'IntensifierActiveShape';
    case 'Unknown_0018_9428'
        newName = 'IntensifierActiveDimensions';
    case 'Unknown_0018_9429'
        newName = 'PhysicalDetectorSize';
    case 'Unknown_0018_9430'
        newName = 'PositionOfIsocenterProjection';
    case 'Unknown_0018_9432'
        newName = 'FieldOfViewSequence';
    case 'Unknown_0018_9433'
        newName = 'FieldOfViewDescription';
    case 'Unknown_0018_9434'
        newName = 'ExposureControlSensingRegionsSequence';
    case 'Unknown_0018_9435'
        newName = 'ExposureControlSensingRegionShape';
    case 'Unknown_0018_9436'
        newName = 'ExposureControlSensingRegionLeftVerticalEdge';
    case 'Unknown_0018_9437'
        newName = 'ExposureControlSensingRegionRightVerticalEdge';
    case 'Unknown_0018_9438'
        newName = 'ExposureControlSensingRegionUpperHorizontalEdge';
    case 'Unknown_0018_9439'
        newName = 'ExposureControlSensingRegionLowerHorizontalEdge';
    case 'Unknown_0018_9440'
        newName = 'CenterOfCircularExposureControlSensingRegion';
    case 'Unknown_0018_9441'
        newName = 'RadiusOfCircularExposureControlSensingRegion';
    case 'Unknown_0018_9442'
        newName = 'VerticesOfPolygonalExposureControlSensingRegion';
    case 'Unknown_0018_9447'
        newName = 'ColumnAngulationPatient';
    case 'Unknown_0018_9449'
        newName = 'BeamAngle';
    case 'Unknown_0018_9451'
        newName = 'FrameDetectorParametersSequence';
    case 'Unknown_0018_9452'
        newName = 'CalculatedAnatomyThickness';
    case 'Unknown_0018_9455'
        newName = 'CalibrationSequence';
    case 'Unknown_0018_9456'
        newName = 'ObjectThicknessSequence';
    case 'Unknown_0018_9457'
        newName = 'PlaneIdentification';
    case 'Unknown_0018_9461'
        newName = 'FieldOfViewDimensionsInFloat';
    case 'Unknown_0018_9462'
        newName = 'IsocenterReferenceSystemSequence';
    case 'Unknown_0018_9463'
        newName = 'PositionerIsocenterPrimaryAngle';
    case 'Unknown_0018_9464'
        newName = 'PositionerIsocenterSecondaryAngle';
    case 'Unknown_0018_9465'
        newName = 'PositionerIsocenterDetectorRotationAngle';
    case 'Unknown_0018_9466'
        newName = 'TableXPositionToIsocenter';
    case 'Unknown_0018_9467'
        newName = 'TableYPositionToIsocenter';
    case 'Unknown_0018_9468'
        newName = 'TableZPositionToIsocenter';
    case 'Unknown_0018_9469'
        newName = 'TableHorizontalRotationAngle';
    case 'Unknown_0018_9470'
        newName = 'TableHeadTiltAngle';
    case 'Unknown_0018_9471'
        newName = 'TableCradleTiltAngle';
    case 'Unknown_0018_9472'
        newName = 'FrameDisplayShutterSequence';
    case 'Unknown_0018_9473'
        newName = 'AcquiredImageAreaDoseProduct';
    case 'Unknown_0018_9474'
        newName = 'CArmPositionerTabletopRelationship';
    case 'Unknown_0018_9476'
        newName = 'XRayGeometrySequence';
    case 'Unknown_0018_9477'
        newName = 'IrradiationEventIdentificationSequence';
    case 'Unknown_0018_A001'
        newName = 'ContributingEquipmentSequence';
    case 'Unknown_0018_A002'
        newName = 'ContributionDateTime';
    case 'Unknown_0018_A003'
        newName = 'ContributionDescription';
    case 'Unknown_0020_0000'
        newName = 'RelationshipGroupLength';
    case 'Unknown_0020_000D'
        newName = 'StudyInstanceUID';
    case 'Unknown_0020_000E'
        newName = 'SeriesInstanceUID';
    case 'Unknown_0020_0010'
        newName = 'StudyID';
    case 'Unknown_0020_0011'
        newName = 'SeriesNumber';
    case 'Unknown_0020_0012'
        newName = 'AcquisitionNumber';
    case 'Unknown_0020_0013'
        newName = 'InstanceNumber';
    case 'Unknown_0020_0014'
        newName = 'IsotopeNumber';
    case 'Unknown_0020_0015'
        newName = 'PhaseNumber';
    case 'Unknown_0020_0016'
        newName = 'IntervalNumber';
    case 'Unknown_0020_0017'
        newName = 'TimeSlotNumber';
    case 'Unknown_0020_0018'
        newName = 'AngleNumber';
    case 'Unknown_0020_0019'
        newName = 'ItemNumber';
    case 'Unknown_0020_0020'
        newName = 'PatientOrientation';
    case 'Unknown_0020_0022'
        newName = 'OverlayNumber';
    case 'Unknown_0020_0024'
        newName = 'CurveNumber';
    case 'Unknown_0020_0026'
        newName = 'LUTNumber';
    case 'Unknown_0020_0030'
        newName = 'ImagePosition';
    case 'Unknown_0020_0032'
        newName = 'ImagePositionPatient';
    case 'Unknown_0020_0035'
        newName = 'ImageOrientation';
    case 'Unknown_0020_0037'
        newName = 'ImageOrientationPatient';
    case 'Unknown_0020_0050'
        newName = 'Location';
    case 'Unknown_0020_0052'
        newName = 'FrameOfReferenceUID';
    case 'Unknown_0020_0060'
        newName = 'Laterality';
    case 'Unknown_0020_0062'
        newName = 'ImageLaterality';
    case 'Unknown_0020_0070'
        newName = 'ImageGeometryType';
    case 'Unknown_0020_0080'
        newName = 'MaskingImage';
    case 'Unknown_0020_0100'
        newName = 'TemporalPositionIdentifier';
    case 'Unknown_0020_0105'
        newName = 'NumberOfTemporalPositions';
    case 'Unknown_0020_0110'
        newName = 'TemporalResolution';
    case 'Unknown_0020_0200'
        newName = 'SynchronizationFrameOfReferenceUID';
    case 'Unknown_0020_1000'
        newName = 'SeriesInStudy';
    case 'Unknown_0020_1001'
        newName = 'AcquisitionsInSeries';
    case 'Unknown_0020_1002'
        newName = 'ImagesInAcquisition';
    case 'Unknown_0020_1003'
        newName = 'ImagesInSeries';
    case 'Unknown_0020_1004'
        newName = 'AcquisitionsInStudy';
    case 'Unknown_0020_1005'
        newName = 'ImagesInStudy';
    case 'Unknown_0020_1020'
        newName = 'Reference';
    case 'Unknown_0020_1040'
        newName = 'PositionReferenceIndicator';
    case 'Unknown_0020_1041'
        newName = 'SliceLocation';
    case 'Unknown_0020_1070'
        newName = 'OtherStudyNumbers';
    case 'Unknown_0020_1200'
        newName = 'NumberOfPatientRelatedStudies';
    case 'Unknown_0020_1202'
        newName = 'NumberOfPatientRelatedSeries';
    case 'Unknown_0020_1204'
        newName = 'NumberOfPatientRelatedInstances';
    case 'Unknown_0020_1206'
        newName = 'NumberOfStudyRelatedSeries';
    case 'Unknown_0020_1208'
        newName = 'NumberOfStudyRelatedInstances';
    case 'Unknown_0020_1209'
        newName = 'NumberOfSeriesRelatedInstances';
    case 'Unknown_0020_31xx'
        newName = 'SourceImageID';
    case 'Unknown_0020_3401'
        newName = 'ModifyingDeviceID';
    case 'Unknown_0020_3402'
        newName = 'ModifiedImageID';
    case 'Unknown_0020_3403'
        newName = 'ModifiedImageDate';
    case 'Unknown_0020_3404'
        newName = 'ModifyingDeviceManufacturer';
    case 'Unknown_0020_3405'
        newName = 'ModifiedImageTime';
    case 'Unknown_0020_3406'
        newName = 'ModifiedImageDescription';
    case 'Unknown_0020_4000'
        newName = 'ImageComments';
    case 'Unknown_0020_5000'
        newName = 'OriginalImageIdentification';
    case 'Unknown_0020_5002'
        newName = 'OriginalImageIdentificationNomenclature';
    case 'Unknown_0020_9056'
        newName = 'StackID';
    case 'Unknown_0020_9057'
        newName = 'InStackPositionNumber';
    case 'Unknown_0020_9071'
        newName = 'FrameAnatomySequence';
    case 'Unknown_0020_9072'
        newName = 'FrameLaterality';
    case 'Unknown_0020_9111'
        newName = 'FrameContentSequence';
    case 'Unknown_0020_9113'
        newName = 'PlanePositionSequence';
    case 'Unknown_0020_9116'
        newName = 'PlaneOrientationSequence';
    case 'Unknown_0020_9128'
        newName = 'TemporalPositionIndex';
    case 'Unknown_0020_9153'
        newName = 'CardiacTriggerDelayTime';
    case 'Unknown_0020_9156'
        newName = 'FrameAcquisitionNumber';
    case 'Unknown_0020_9157'
        newName = 'DimensionIndexValues';
    case 'Unknown_0020_9158'
        newName = 'FrameComments';
    case 'Unknown_0020_9161'
        newName = 'ConcatenationUID';
    case 'Unknown_0020_9162'
        newName = 'InConcatenationNumber';
    case 'Unknown_0020_9163'
        newName = 'InConcatenationTotalNumber';
    case 'Unknown_0020_9164'
        newName = 'DimensionOrganizationUID';
    case 'Unknown_0020_9165'
        newName = 'DimensionIndexPointer';
    case 'Unknown_0020_9167'
        newName = 'FunctionalGroupPointer';
    case 'Unknown_0020_9213'
        newName = 'DimensionIndexPrivateCreator';
    case 'Unknown_0020_9221'
        newName = 'DimensionOrganizationSequence';
    case 'Unknown_0020_9222'
        newName = 'DimensionIndexSequence';
    case 'Unknown_0020_9228'
        newName = 'ConcatenationFrameOffsetNumber';
    case 'Unknown_0020_9238'
        newName = 'FunctionalGroupPrivateCreator';
    case 'Unknown_0020_9251'
        newName = 'RRIntervalTimeMeasured';
    case 'Unknown_0020_9253'
        newName = 'RespiratoryTriggerSequence';
    case 'Unknown_0020_9254'
        newName = 'RespiratoryIntervalTime';
    case 'Unknown_0020_9255'
        newName = 'RespiratoryTriggerDelayTime';
    case 'Unknown_0020_9256'
        newName = 'RespiratoryTriggerDelayThreshold';
    case 'Unknown_0020_9421'
        newName = 'DimensionDescriptionLabel';
    case 'Unknown_0020_9450'
        newName = 'PatientOrientationInFrameSequence';
    case 'Unknown_0020_9453'
        newName = 'FrameLabel';
    case 'Unknown_0022_0001'
        newName = 'LightPathFilterPassThroughWavelength';
    case 'Unknown_0022_0002'
        newName = 'LightPathFilterPassBand';
    case 'Unknown_0022_0003'
        newName = 'ImagePathFilterPassThroughWavelength';
    case 'Unknown_0022_0004'
        newName = 'ImagePathFilterPassBand';
    case 'Unknown_0022_0005'
        newName = 'PatientEyeMovementCommanded';
    case 'Unknown_0022_0006'
        newName = 'PatientEyeMovementCommandedCodeSequence';
    case 'Unknown_0022_0007'
        newName = 'SphericalLensPower';
    case 'Unknown_0022_0008'
        newName = 'CylinderLensPower';
    case 'Unknown_0022_0009'
        newName = 'CylinderAxis';
    case 'Unknown_0022_000A'
        newName = 'EmmetropicMagnification';
    case 'Unknown_0022_000B'
        newName = 'IntraOcularPressure';
    case 'Unknown_0022_000C'
        newName = 'HorizontalFieldOfView';
    case 'Unknown_0022_000D'
        newName = 'PupilDilated';
    case 'Unknown_0022_000E'
        newName = 'DegreeOfDilation';
    case 'Unknown_0022_0010'
        newName = 'StereoBaselineAngle';
    case 'Unknown_0022_0011'
        newName = 'StereoBaselineDisplacement';
    case 'Unknown_0022_0012'
        newName = 'StereoHorizontalPixelOffset';
    case 'Unknown_0022_0013'
        newName = 'StereoVerticalPixelOffset';
    case 'Unknown_0022_0014'
        newName = 'StereoRotation';
    case 'Unknown_0022_0015'
        newName = 'AcquisitionDeviceTypeCodeSequence';
    case 'Unknown_0022_0016'
        newName = 'IlluminationTypeCodeSequence';
    case 'Unknown_0022_0017'
        newName = 'LightPathFilterTypeStackCodeSequence';
    case 'Unknown_0022_0018'
        newName = 'ImagePathFilterTypeStackCodeSequence';
    case 'Unknown_0022_0019'
        newName = 'LensesCodeSequence';
    case 'Unknown_0022_001A'
        newName = 'ChannelDescriptionCodeSequence';
    case 'Unknown_0022_001B'
        newName = 'RefractiveStateSequence';
    case 'Unknown_0022_001C'
        newName = 'MydriaticAgentCodeSequence';
    case 'Unknown_0022_001D'
        newName = 'RelativeImagePositionCodeSequence';
    case 'Unknown_0022_0020'
        newName = 'StereoPairsSequence';
    case 'Unknown_0022_0021'
        newName = 'LeftImageSequence';
    case 'Unknown_0022_0022'
        newName = 'RightImageSequence';
    case 'Unknown_0028_0000'
        newName = 'ImagePresentationGroupLength';
    case 'Unknown_0028_0002'
        newName = 'SamplesPerPixel';
    case 'Unknown_0028_0003'
        newName = 'SamplesPerPixelUsed';
    case 'Unknown_0028_0004'
        newName = 'PhotometricInterpretation';
    case 'Unknown_0028_0005'
        newName = 'ImageDimensions';
    case 'Unknown_0028_0006'
        newName = 'PlanarConfiguration';
    case 'Unknown_0028_0008'
        newName = 'NumberOfFrames';
    case 'Unknown_0028_0009'
        newName = 'FrameIncrementPointer';
    case 'Unknown_0028_000A'
        newName = 'FrameDimensionPointer';
    case 'Unknown_0028_0010'
        newName = 'Rows';
    case 'Unknown_0028_0011'
        newName = 'Columns';
    case 'Unknown_0028_0012'
        newName = 'Planes';
    case 'Unknown_0028_0014'
        newName = 'UltrasoundColorDataPresent';
    case 'Unknown_0028_0030'
        newName = 'PixelSpacing';
    case 'Unknown_0028_0031'
        newName = 'ZoomFactor';
    case 'Unknown_0028_0032'
        newName = 'ZoomCenter';
    case 'Unknown_0028_0034'
        newName = 'PixelAspectRatio';
    case 'Unknown_0028_0040'
        newName = 'ImageFormat';
    case 'Unknown_0028_0050'
        newName = 'ManipulatedImage';
    case 'Unknown_0028_0051'
        newName = 'CorrectedImage';
    case 'Unknown_0028_005F'
        newName = 'CompressionRecognitionCode';
    case 'Unknown_0028_0060'
        newName = 'CompressionCode';
    case 'Unknown_0028_0061'
        newName = 'CompressionOriginator';
    case 'Unknown_0028_0062'
        newName = 'CompressionLabel';
    case 'Unknown_0028_0063'
        newName = 'CompressionDescription';
    case 'Unknown_0028_0065'
        newName = 'CompressionSequence';
    case 'Unknown_0028_0066'
        newName = 'CompressionStepPointers';
    case 'Unknown_0028_0068'
        newName = 'RepeatInterval';
    case 'Unknown_0028_0069'
        newName = 'BitsGrouped';
    case 'Unknown_0028_0070'
        newName = 'PerimeterTable';
    case 'Unknown_0028_0071'
        newName = 'PerimeterValue';
    case 'Unknown_0028_0080'
        newName = 'PredictorRows';
    case 'Unknown_0028_0081'
        newName = 'PredictorColumns';
    case 'Unknown_0028_0082'
        newName = 'PredictorConstants';
    case 'Unknown_0028_0090'
        newName = 'BlockedPixels';
    case 'Unknown_0028_0091'
        newName = 'BlockRows';
    case 'Unknown_0028_0092'
        newName = 'BlockColumns';
    case 'Unknown_0028_0093'
        newName = 'RowOverlap';
    case 'Unknown_0028_0094'
        newName = 'ColumnOverlap';
    case 'Unknown_0028_0100'
        newName = 'BitsAllocated';
    case 'Unknown_0028_0101'
        newName = 'BitsStored';
    case 'Unknown_0028_0102'
        newName = 'HighBit';
    case 'Unknown_0028_0103'
        newName = 'PixelRepresentation';
    case 'Unknown_0028_0104'
        newName = 'SmallestValidPixelValue';
    case 'Unknown_0028_0105'
        newName = 'LargestValidPixelValue';
    case 'Unknown_0028_0106'
        newName = 'SmallestImagePixelValue';
    case 'Unknown_0028_0107'
        newName = 'LargestImagePixelValue';
    case 'Unknown_0028_0108'
        newName = 'SmallestPixelValueInSeries';
    case 'Unknown_0028_0109'
        newName = 'LargestPixelValueInSeries';
    case 'Unknown_0028_0110'
        newName = 'SmallestPixelValueInPlane';
    case 'Unknown_0028_0111'
        newName = 'LargestPixelValueInPlane';
    case 'Unknown_0028_0120'
        newName = 'PixelPaddingValue';
    case 'Unknown_0028_0121'
        newName = 'PixelPaddingRangeLimit';
    case 'Unknown_0028_0200'
        newName = 'ImageLocation';
    case 'Unknown_0028_0300'
        newName = 'QualityControlImage';
    case 'Unknown_0028_0301'
        newName = 'BurnedInAnnotation';
    case 'Unknown_0028_0400'
        newName = 'TransformLabel';
    case 'Unknown_0028_0401'
        newName = 'TransformVersionNumber';
    case 'Unknown_0028_0402'
        newName = 'PixelSpacingCalibrationType';
    case 'Unknown_0028_0403'
        newName = 'SequenceOfCompressedData';
    case 'Unknown_0028_0404'
        newName = 'PixelSpacingCalibrationDescription';
    case 'Unknown_0028_0700'
        newName = 'DCTLabel';
    case 'Unknown_0028_0701'
        newName = 'DataBlockDescription';
    case 'Unknown_0028_0702'
        newName = 'DataBlock';
    case 'Unknown_0028_0710'
        newName = 'NormalizationFactorFormat';
    case 'Unknown_0028_0720'
        newName = 'ZonalMapNumberFormat';
    case 'Unknown_0028_0721'
        newName = 'ZonalMapLocation';
    case 'Unknown_0028_0722'
        newName = 'ZonalMapFormat';
    case 'Unknown_0028_0730'
        newName = 'AdaptiveMapFormat';
    case 'Unknown_0028_0740'
        newName = 'CodeNumberFormat';
    case 'Unknown_0028_1040'
        newName = 'PixelIntensityRelationship';
    case 'Unknown_0028_1041'
        newName = 'PixelIntensityRelationshipSign';
    case 'Unknown_0028_1050'
        newName = 'WindowCenter';
    case 'Unknown_0028_1051'
        newName = 'WindowWidth';
    case 'Unknown_0028_1052'
        newName = 'RescaleIntercept';
    case 'Unknown_0028_1053'
        newName = 'RescaleSlope';
    case 'Unknown_0028_1054'
        newName = 'RescaleType';
    case 'Unknown_0028_1055'
        newName = 'WindowCenterWidthExplanation';
    case 'Unknown_0028_1056'
        newName = 'VOILUTFunction';
    case 'Unknown_0028_1080'
        newName = 'GrayScale';
    case 'Unknown_0028_1090'
        newName = 'RecommendedViewingMode';
    case 'Unknown_0028_1100'
        newName = 'GrayLookupTableDescriptor';
    case 'Unknown_0028_1101'
        newName = 'RedPaletteColorLookupTableDescriptor';
    case 'Unknown_0028_1102'
        newName = 'GreenPaletteColorLookupTableDescriptor';
    case 'Unknown_0028_1103'
        newName = 'BluePaletteColorLookupTableDescriptor';
    case 'Unknown_0028_1111'
        newName = 'LargeRedPaletteColorLookupTableDescriptor';
    case 'Unknown_0028_1112'
        newName = 'LargeGreenPaletteColorLookupTableDescriptor';
    case 'Unknown_0028_1113'
        newName = 'LargeBluePaletteColorLookupTableDescriptor';
    case 'Unknown_0028_1199'
        newName = 'PaletteColorLookupTableUID';
    case 'Unknown_0028_1200'
        newName = 'GrayLookupTableData';
    case 'Unknown_0028_1201'
        newName = 'RedPaletteColorLookupTableData';
    case 'Unknown_0028_1202'
        newName = 'GreenPaletteColorLookupTableData';
    case 'Unknown_0028_1203'
        newName = 'BluePaletteColorLookupTableData';
    case 'Unknown_0028_1211'
        newName = 'LargeRedPaletteColorLookupTableData';
    case 'Unknown_0028_1212'
        newName = 'LargeGreenPaletteColorLookupTableData';
    case 'Unknown_0028_1213'
        newName = 'LargeBluePaletteColorLookupTableData';
    case 'Unknown_0028_1214'
        newName = 'LargePaletteColorLookupTableUID';
    case 'Unknown_0028_1221'
        newName = 'SegmentedRedPaletteColorLookupTableData';
    case 'Unknown_0028_1222'
        newName = 'SegmentedGreenPaletteColorLookupTableData';
    case 'Unknown_0028_1223'
        newName = 'SegmentedBluePaletteColorLookupTableData';
    case 'Unknown_0028_1300'
        newName = 'ImplantPresent';
    case 'Unknown_0028_1350'
        newName = 'PartialView';
    case 'Unknown_0028_1351'
        newName = 'PartialViewDescription';
    case 'Unknown_0028_1352'
        newName = 'PartialViewCodeSequence';
    case 'Unknown_0028_135A'
        newName = 'SpatialLocationsPreserved';
    case 'Unknown_0028_2000'
        newName = 'ICCProfile';
    case 'Unknown_0028_2110'
        newName = 'LossyImageCompression';
    case 'Unknown_0028_2112'
        newName = 'LossyImageCompressionRatio';
    case 'Unknown_0028_2114'
        newName = 'LossyImageCompressionMethod';
    case 'Unknown_0028_3000'
        newName = 'ModalityLUTSequence';
    case 'Unknown_0028_3002'
        newName = 'LUTDescriptor';
    case 'Unknown_0028_3003'
        newName = 'LUTExplanation';
    case 'Unknown_0028_3004'
        newName = 'ModalityLUTType';
    case 'Unknown_0028_3006'
        newName = 'LUTData';
    case 'Unknown_0028_3010'
        newName = 'VOILUTSequence';
    case 'Unknown_0028_3110'
        newName = 'SoftcopyVOILUTSequence';
    case 'Unknown_0028_4000'
        newName = 'ImagePresentationComments';
    case 'Unknown_0028_5000'
        newName = 'BiplaneAcquisitionSequence';
    case 'Unknown_0028_6010'
        newName = 'RepresentativeFrameNumber';
    case 'Unknown_0028_6020'
        newName = 'FrameNumbersOfInterest';
    case 'Unknown_0028_6022'
        newName = 'FrameOfInterestDescription';
    case 'Unknown_0028_6023'
        newName = 'FrameOfInterestType';
    case 'Unknown_0028_6030'
        newName = 'MaskPointer';
    case 'Unknown_0028_6040'
        newName = 'RWavePointer';
    case 'Unknown_0028_6100'
        newName = 'MaskSubtractionSequence';
    case 'Unknown_0028_6101'
        newName = 'MaskOperation';
    case 'Unknown_0028_6102'
        newName = 'ApplicableFrameRange';
    case 'Unknown_0028_6110'
        newName = 'MaskFrameNumbers';
    case 'Unknown_0028_6112'
        newName = 'ContrastFrameAveraging';
    case 'Unknown_0028_6114'
        newName = 'MaskSubPixelShift';
    case 'Unknown_0028_6120'
        newName = 'TIDOffset';
    case 'Unknown_0028_6190'
        newName = 'MaskOperationExplanation';
    case 'Unknown_0028_7FE0'
        newName = 'PixelDataProviderURL';
    case 'Unknown_0028_9001'
        newName = 'DataPointRows';
    case 'Unknown_0028_9002'
        newName = 'DataPointColumns';
    case 'Unknown_0028_9003'
        newName = 'SignalDomainColumns';
    case 'Unknown_0028_9099'
        newName = 'LargestMonochromePixelValue';
    case 'Unknown_0028_9108'
        newName = 'DataRepresentation';
    case 'Unknown_0028_9110'
        newName = 'PixelMeasuresSequence';
    case 'Unknown_0028_9132'
        newName = 'FrameVOILUTSequence';
    case 'Unknown_0028_9145'
        newName = 'PixelValueTransformationSequence';
    case 'Unknown_0028_9235'
        newName = 'SignalDomainRows';
    case 'Unknown_0028_9411'
        newName = 'DisplayFilterPercentage';
    case 'Unknown_0028_9415'
        newName = 'FramePixelShiftSequence';
    case 'Unknown_0028_9416'
        newName = 'SubtractionItemID';
    case 'Unknown_0028_9422'
        newName = 'PixelIntensityRelationshipLUTSequence';
    case 'Unknown_0028_9443'
        newName = 'FramePixelDataPropertiesSequence';
    case 'Unknown_0028_9444'
        newName = 'GeometricalProperties';
    case 'Unknown_0028_9445'
        newName = 'GeometricMaximumDistortion';
    case 'Unknown_0028_9446'
        newName = 'ImageProcessingApplied';
    case 'Unknown_0028_9454'
        newName = 'MaskSelectionMode';
    case 'Unknown_0028_9474'
        newName = 'LUTFunction';
    case 'Unknown_0032_0000'
        newName = 'StudyGroupLength';
    case 'Unknown_0032_000A'
        newName = 'StudyStatusID';
    case 'Unknown_0032_000C'
        newName = 'StudyPriorityID';
    case 'Unknown_0032_0012'
        newName = 'StudyIDIssuer';
    case 'Unknown_0032_0032'
        newName = 'StudyVerifiedDate';
    case 'Unknown_0032_0033'
        newName = 'StudyVerifiedTime';
    case 'Unknown_0032_0034'
        newName = 'StudyReadDate';
    case 'Unknown_0032_0035'
        newName = 'StudyReadTime';
    case 'Unknown_0032_1000'
        newName = 'ScheduledStudyStartDate';
    case 'Unknown_0032_1001'
        newName = 'ScheduledStudyStartTime';
    case 'Unknown_0032_1010'
        newName = 'ScheduledStudyStopDate';
    case 'Unknown_0032_1011'
        newName = 'ScheduledStudyStopTime';
    case 'Unknown_0032_1020'
        newName = 'ScheduledStudyLocation';
    case 'Unknown_0032_1021'
        newName = 'ScheduledStudyLocationAETitle';
    case 'Unknown_0032_1030'
        newName = 'ReasonForStudy';
    case 'Unknown_0032_1031'
        newName = 'RequestingPhysicianIdentificationSequence';
    case 'Unknown_0032_1032'
        newName = 'RequestingPhysician';
    case 'Unknown_0032_1033'
        newName = 'RequestingService';
    case 'Unknown_0032_1040'
        newName = 'StudyArrivalDate';
    case 'Unknown_0032_1041'
        newName = 'StudyArrivalTime';
    case 'Unknown_0032_1050'
        newName = 'StudyCompletionDate';
    case 'Unknown_0032_1051'
        newName = 'StudyCompletionTime';
    case 'Unknown_0032_1055'
        newName = 'StudyComponentStatusID';
    case 'Unknown_0032_1060'
        newName = 'RequestedProcedureDescription';
    case 'Unknown_0032_1064'
        newName = 'RequestedProcedureCodeSequence';
    case 'Unknown_0032_1070'
        newName = 'RequestedContrastAgent';
    case 'Unknown_0032_4000'
        newName = 'StudyComments';
    case 'Unknown_0038_0000'
        newName = 'VisitGroupLength';
    case 'Unknown_0038_0004'
        newName = 'ReferencedPatientAliasSequence';
    case 'Unknown_0038_0008'
        newName = 'VisitStatusID';
    case 'Unknown_0038_0010'
        newName = 'AdmissionID';
    case 'Unknown_0038_0011'
        newName = 'IssuerOfAdmissionID';
    case 'Unknown_0038_0016'
        newName = 'RouteOfAdmissions';
    case 'Unknown_0038_001A'
        newName = 'ScheduledAdmissionDate';
    case 'Unknown_0038_001B'
        newName = 'ScheduledAdmissionTime';
    case 'Unknown_0038_001C'
        newName = 'ScheduledDischargeDate';
    case 'Unknown_0038_001D'
        newName = 'ScheduledDischargeTime';
    case 'Unknown_0038_001E'
        newName = 'ScheduledPatientInstitutionResidence';
    case 'Unknown_0038_0020'
        newName = 'AdmittingDate';
    case 'Unknown_0038_0021'
        newName = 'AdmittingTime';
    case 'Unknown_0038_0030'
        newName = 'DischargeDate';
    case 'Unknown_0038_0032'
        newName = 'DischargeTime';
    case 'Unknown_0038_0040'
        newName = 'DischargeDiagnosisDescription';
    case 'Unknown_0038_0044'
        newName = 'DischargeDiagnosisCodeSequence';
    case 'Unknown_0038_0050'
        newName = 'SpecialNeeds';
    case 'Unknown_0038_0100'
        newName = 'PertinentDocumentsSequence';
    case 'Unknown_0038_0300'
        newName = 'CurrentPatientLocation';
    case 'Unknown_0038_0400'
        newName = 'PatientInstitutionResidence';
    case 'Unknown_0038_0500'
        newName = 'PatientState';
    case 'Unknown_0038_0502'
        newName = 'PatientClinicalTrialParticipationSequence';
    case 'Unknown_0038_4000'
        newName = 'VisitComments';
    case 'Unknown_003A_0004'
        newName = 'WaveformOriginality';
    case 'Unknown_003A_0005'
        newName = 'NumberOfWaveformChannels';
    case 'Unknown_003A_0010'
        newName = 'NumberOfWaveformSamples';
    case 'Unknown_003A_001A'
        newName = 'SamplingFrequency';
    case 'Unknown_003A_0020'
        newName = 'MultiplexGroupLabel';
    case 'Unknown_003A_0200'
        newName = 'ChannelDefinitionSequence';
    case 'Unknown_003A_0202'
        newName = 'WaveformChannelNumber';
    case 'Unknown_003A_0203'
        newName = 'ChannelLabel';
    case 'Unknown_003A_0205'
        newName = 'ChannelStatus';
    case 'Unknown_003A_0208'
        newName = 'ChannelSourceSequence';
    case 'Unknown_003A_0209'
        newName = 'ChannelSourceModifiersSequence';
    case 'Unknown_003A_020A'
        newName = 'SourceWaveformSequence';
    case 'Unknown_003A_020C'
        newName = 'ChannelDerivationDescription';
    case 'Unknown_003A_0210'
        newName = 'ChannelSensitivity';
    case 'Unknown_003A_0211'
        newName = 'ChannelSensitivityUnitsSequence';
    case 'Unknown_003A_0212'
        newName = 'ChannelSensitivityCorrectionFactor';
    case 'Unknown_003A_0213'
        newName = 'ChannelBaseline';
    case 'Unknown_003A_0214'
        newName = 'ChannelTimeSkew';
    case 'Unknown_003A_0215'
        newName = 'ChannelSampleSkew';
    case 'Unknown_003A_0218'
        newName = 'ChannelOffset';
    case 'Unknown_003A_021A'
        newName = 'WaveformBitsStored';
    case 'Unknown_003A_0220'
        newName = 'FilterLowFrequency';
    case 'Unknown_003A_0221'
        newName = 'FilterHighFrequency';
    case 'Unknown_003A_0222'
        newName = 'NotchFilterFrequency';
    case 'Unknown_003A_0223'
        newName = 'NotchFilterBandwidth';
    case 'Unknown_003A_0300'
        newName = 'MultiplexedAudioChannelsDescriptionCodeSequence';
    case 'Unknown_003A_0301'
        newName = 'ChannelIdentificationCode';
    case 'Unknown_003A_0302'
        newName = 'ChannelMode';
    case 'Unknown_0040_0001'
        newName = 'ScheduledStationAETitle';
    case 'Unknown_0040_0002'
        newName = 'ScheduledProcedureStepStartDate';
    case 'Unknown_0040_0003'
        newName = 'ScheduledProcedureStepStartTime';
    case 'Unknown_0040_0004'
        newName = 'ScheduledProcedureStepEndDate';
    case 'Unknown_0040_0005'
        newName = 'ScheduledProcedureStepEndTime';
    case 'Unknown_0040_0006'
        newName = 'ScheduledPerformingPhysicianName';
    case 'Unknown_0040_0007'
        newName = 'ScheduledProcedureStepDescription';
    case 'Unknown_0040_0008'
        newName = 'ScheduledProtocolCodeSequence';
    case 'Unknown_0040_0009'
        newName = 'ScheduledProcedureStepID';
    case 'Unknown_0040_000A'
        newName = 'StageCodeSequence';
    case 'Unknown_0040_000B'
        newName = 'ScheduledPerformingPhysicianIdentificationSequence';
    case 'Unknown_0040_0010'
        newName = 'ScheduledStationName';
    case 'Unknown_0040_0011'
        newName = 'ScheduledProcedureStepLocation';
    case 'Unknown_0040_0012'
        newName = 'PreMedication';
    case 'Unknown_0040_0020'
        newName = 'ScheduledProcedureStepStatus';
    case 'Unknown_0040_0100'
        newName = 'ScheduledProcedureStepSequence';
    case 'Unknown_0040_0220'
        newName = 'ReferencedNonImageCompositeSOPInstanceSequence';
    case 'Unknown_0040_0241'
        newName = 'PerformedStationAETitle';
    case 'Unknown_0040_0242'
        newName = 'PerformedStationName';
    case 'Unknown_0040_0243'
        newName = 'PerformedLocation';
    case 'Unknown_0040_0244'
        newName = 'PerformedProcedureStepStartDate';
    case 'Unknown_0040_0245'
        newName = 'PerformedProcedureStepStartTime';
    case 'Unknown_0040_0250'
        newName = 'PerformedProcedureStepEndDate';
    case 'Unknown_0040_0251'
        newName = 'PerformedProcedureStepEndTime';
    case 'Unknown_0040_0252'
        newName = 'PerformedProcedureStepStatus';
    case 'Unknown_0040_0253'
        newName = 'PerformedProcedureStepID';
    case 'Unknown_0040_0254'
        newName = 'PerformedProcedureStepDescription';
    case 'Unknown_0040_0255'
        newName = 'PerformedProcedureTypeDescription';
    case 'Unknown_0040_0260'
        newName = 'PerformedProtocolCodeSequence';
    case 'Unknown_0040_0270'
        newName = 'ScheduledStepAttributesSequence';
    case 'Unknown_0040_0275'
        newName = 'RequestAttributesSequence';
    case 'Unknown_0040_0280'
        newName = 'CommentsOnPerformedProcedureStep';
    case 'Unknown_0040_0281'
        newName = 'PerformedProcedureStepDiscontinuationReasonCodeSequence';
    case 'Unknown_0040_0293'
        newName = 'QuantitySequence';
    case 'Unknown_0040_0294'
        newName = 'Quantity';
    case 'Unknown_0040_0295'
        newName = 'MeasuringUnitsSequence';
    case 'Unknown_0040_0296'
        newName = 'BillingItemSequence';
    case 'Unknown_0040_0300'
        newName = 'TotalTimeOfFluoroscopy';
    case 'Unknown_0040_0301'
        newName = 'TotalNumberOfExposures';
    case 'Unknown_0040_0302'
        newName = 'EntranceDose';
    case 'Unknown_0040_0303'
        newName = 'ExposedArea';
    case 'Unknown_0040_0306'
        newName = 'DistanceSourceToEntrance';
    case 'Unknown_0040_0307'
        newName = 'DistanceSourceToSupport';
    case 'Unknown_0040_030E'
        newName = 'ExposureDoseSequence';
    case 'Unknown_0040_0310'
        newName = 'CommentsOnRadiationDose';
    case 'Unknown_0040_0312'
        newName = 'XRayOutput';
    case 'Unknown_0040_0314'
        newName = 'HalfValueLayer';
    case 'Unknown_0040_0316'
        newName = 'OrganDose';
    case 'Unknown_0040_0318'
        newName = 'OrganExposed';
    case 'Unknown_0040_0320'
        newName = 'BillingProcedureStepSequence';
    case 'Unknown_0040_0321'
        newName = 'FilmConsumptionSequence';
    case 'Unknown_0040_0324'
        newName = 'BillingSuppliesAndDevicesSequence';
    case 'Unknown_0040_0330'
        newName = 'ReferencedProcedureStepSequence';
    case 'Unknown_0040_0340'
        newName = 'PerformedSeriesSequence';
    case 'Unknown_0040_0400'
        newName = 'CommentsOnScheduledProcedureStep';
    case 'Unknown_0040_0440'
        newName = 'ProtocolContextSequence';
    case 'Unknown_0040_0441'
        newName = 'ContentItemModifierSequence';
    case 'Unknown_0040_050A'
        newName = 'SpecimenAccessionNumber';
    case 'Unknown_0040_0550'
        newName = 'SpecimenSequence';
    case 'Unknown_0040_0551'
        newName = 'SpecimenIdentifier';
    case 'Unknown_0040_0552'
        newName = 'SpecimenDescriptionSequenceTrial';
    case 'Unknown_0040_0553'
        newName = 'SpecimenDescriptionTrial';
    case 'Unknown_0040_0555'
        newName = 'AcquisitionContextSequence';
    case 'Unknown_0040_0556'
        newName = 'AcquisitionContextDescription';
    case 'Unknown_0040_059A'
        newName = 'SpecimenTypeCodeSequence';
    case 'Unknown_0040_06FA'
        newName = 'SlideIdentifier';
    case 'Unknown_0040_071A'
        newName = 'ImageCenterPointCoordinatesSequence';
    case 'Unknown_0040_072A'
        newName = 'XOffsetInSlideCoordinateSystem';
    case 'Unknown_0040_073A'
        newName = 'YOffsetInSlideCoordinateSystem';
    case 'Unknown_0040_074A'
        newName = 'ZOffsetInSlideCoordinateSystem';
    case 'Unknown_0040_08D8'
        newName = 'PixelSpacingSequence';
    case 'Unknown_0040_08DA'
        newName = 'CoordinateSystemAxisCodeSequence';
    case 'Unknown_0040_08EA'
        newName = 'MeasurementUnitsCodeSequence';
    case 'Unknown_0040_09F8'
        newName = 'VitalStainCodeSequenceTrial';
    case 'Unknown_0040_1001'
        newName = 'RequestedProcedureID';
    case 'Unknown_0040_1002'
        newName = 'ReasonForRequestedProcedure';
    case 'Unknown_0040_1003'
        newName = 'RequestedProcedurePriority';
    case 'Unknown_0040_1004'
        newName = 'PatientTransportArrangements';
    case 'Unknown_0040_1005'
        newName = 'RequestedProcedureLocation';
    case 'Unknown_0040_1006'
        newName = 'PlacerOrderNumberOfProcedure';
    case 'Unknown_0040_1007'
        newName = 'FillerOrderNumberOfProcedure';
    case 'Unknown_0040_1008'
        newName = 'ConfidentialityCode';
    case 'Unknown_0040_1009'
        newName = 'ReportingPriority';
    case 'Unknown_0040_100A'
        newName = 'ReasonForRequestedProcedureCodeSequence';
    case 'Unknown_0040_1010'
        newName = 'NamesOfIntendedRecipientsOfResults';
    case 'Unknown_0040_1011'
        newName = 'IntendedRecipientsOfResultsIdentificationSequence';
    case 'Unknown_0040_1101'
        newName = 'PersonIdentificationCodeSequence';
    case 'Unknown_0040_1102'
        newName = 'PersonAddress';
    case 'Unknown_0040_1103'
        newName = 'PersonTelephoneNumbers';
    case 'Unknown_0040_1400'
        newName = 'RequestedProcedureComments';
    case 'Unknown_0040_2001'
        newName = 'ReasonForImagingServiceRequest';
    case 'Unknown_0040_2004'
        newName = 'IssueDateOfImagingServiceRequest';
    case 'Unknown_0040_2005'
        newName = 'IssueTimeOfImagingServiceRequest';
    case 'Unknown_0040_2006'
        newName = 'PlacerOrderNumberOfImagingServiceRequestRetired';
    case 'Unknown_0040_2007'
        newName = 'FillerOrderNumberOfImagingServiceRequestRetired';
    case 'Unknown_0040_2008'
        newName = 'OrderEnteredBy';
    case 'Unknown_0040_2009'
        newName = 'OrderEntererLocation';
    case 'Unknown_0040_2010'
        newName = 'OrderCallbackPhoneNumber';
    case 'Unknown_0040_2016'
        newName = 'PlacerOrderNumberOfImagingServiceRequest';
    case 'Unknown_0040_2017'
        newName = 'FillerOrderNumberOfImagingServiceRequest';
    case 'Unknown_0040_2400'
        newName = 'ImagingServiceRequestComments';
    case 'Unknown_0040_3001'
        newName = 'ConfidentialityConstraintOnPatientDataDescription';
    case 'Unknown_0040_4001'
        newName = 'GeneralPurposeScheduledProcedureStepStatus';
    case 'Unknown_0040_4002'
        newName = 'GeneralPurposePerformedProcedureStepStatus';
    case 'Unknown_0040_4003'
        newName = 'GeneralPurposeScheduledProcedureStepPriority';
    case 'Unknown_0040_4004'
        newName = 'ScheduledProcessingApplicationsCodeSequence';
    case 'Unknown_0040_4005'
        newName = 'ScheduledProcedureStepStartDateAndTime';
    case 'Unknown_0040_4006'
        newName = 'MultipleCopiesFlag';
    case 'Unknown_0040_4007'
        newName = 'PerformedProcessingApplicationsCodeSequence';
    case 'Unknown_0040_4009'
        newName = 'HumanPerformerCodeSequence';
    case 'Unknown_0040_4010'
        newName = 'ScheduledProcedureStepModificationDateAndTime';
    case 'Unknown_0040_4011'
        newName = 'ExpectedCompletionDateAndTime';
    case 'Unknown_0040_4015'
        newName = 'ResultingGeneralPurposePerformedProcedureStepsSequence';
    case 'Unknown_0040_4016'
        newName = 'ReferencedGeneralPurposeScheduledProcedureStepSequence';
    case 'Unknown_0040_4018'
        newName = 'ScheduledWorkitemCodeSequence';
    case 'Unknown_0040_4019'
        newName = 'PerformedWorkitemCodeSequence';
    case 'Unknown_0040_4020'
        newName = 'InputAvailabilityFlag';
    case 'Unknown_0040_4021'
        newName = 'InputInformationSequence';
    case 'Unknown_0040_4022'
        newName = 'RelevantInformationSequence';
    case 'Unknown_0040_4023'
        newName = 'ReferencedGeneralPurposeScheduledProcedureStepTransactionUID';
    case 'Unknown_0040_4025'
        newName = 'ScheduledStationNameCodeSequence';
    case 'Unknown_0040_4026'
        newName = 'ScheduledStationClassCodeSequence';
    case 'Unknown_0040_4027'
        newName = 'ScheduledStationGeographicLocationCodeSequence';
    case 'Unknown_0040_4028'
        newName = 'PerformedStationNameCodeSequence';
    case 'Unknown_0040_4029'
        newName = 'PerformedStationClassCodeSequence';
    case 'Unknown_0040_4030'
        newName = 'PerformedStationGeographicLocationCodeSequence';
    case 'Unknown_0040_4031'
        newName = 'RequestedSubsequentWorkitemCodeSequence';
    case 'Unknown_0040_4032'
        newName = 'NonDICOMOutputCodeSequence';
    case 'Unknown_0040_4033'
        newName = 'OutputInformationSequence';
    case 'Unknown_0040_4034'
        newName = 'ScheduledHumanPerformersSequence';
    case 'Unknown_0040_4035'
        newName = 'ActualHumanPerformersSequence';
    case 'Unknown_0040_4036'
        newName = 'HumanPerformersOrganization';
    case 'Unknown_0040_4037'
        newName = 'HumanPerformersName';
    case 'Unknown_0040_8302'
        newName = 'EntranceDoseInmGy';
    case 'Unknown_0040_9094'
        newName = 'ReferencedImageRealWorldValueMappingSequence';
    case 'Unknown_0040_9096'
        newName = 'RealWorldValueMappingSequence';
    case 'Unknown_0040_9098'
        newName = 'PixelValueMappingCodeSequence';
    case 'Unknown_0040_9210'
        newName = 'LUTLabel';
    case 'Unknown_0040_9211'
        newName = 'RealWorldValueLastValueMapped';
    case 'Unknown_0040_9212'
        newName = 'RealWorldValueLUTData';
    case 'Unknown_0040_9216'
        newName = 'RealWorldValueFirstValueMapped';
    case 'Unknown_0040_9224'
        newName = 'RealWorldValueIntercept';
    case 'Unknown_0040_9225'
        newName = 'RealWorldValueSlope';
    case 'Unknown_0040_A010'
        newName = 'RelationshipType';
    case 'Unknown_0040_A027'
        newName = 'VerifyingOrganization';
    case 'Unknown_0040_A030'
        newName = 'VerificationDateTime';
    case 'Unknown_0040_A032'
        newName = 'ObservationDateTime';
    case 'Unknown_0040_A040'
        newName = 'ValueType';
    case 'Unknown_0040_A043'
        newName = 'ConceptNameCodeSequence';
    case 'Unknown_0040_A050'
        newName = 'ContinuityOfContent';
    case 'Unknown_0040_A073'
        newName = 'VerifyingObserverSequence';
    case 'Unknown_0040_A075'
        newName = 'VerifyingObserverName';
    case 'Unknown_0040_A078'
        newName = 'AuthorObserverSequence';
    case 'Unknown_0040_A07A'
        newName = 'ParticipantSequence';
    case 'Unknown_0040_A07C'
        newName = 'CustodialOrganizationSequence';
    case 'Unknown_0040_A080'
        newName = 'ParticipationType';
    case 'Unknown_0040_A082'
        newName = 'ParticipationDatetime';
    case 'Unknown_0040_A084'
        newName = 'ObserverType';
    case 'Unknown_0040_A088'
        newName = 'VerifyingObserverIdentificationCodeSequence';
    case 'Unknown_0040_A090'
        newName = 'EquivalentCDADocumentSequence';
    case 'Unknown_0040_A0B0'
        newName = 'ReferencedWaveformChannels';
    case 'Unknown_0040_A120'
        newName = 'DateTime';
    case 'Unknown_0040_A121'
        newName = 'Date';
    case 'Unknown_0040_A122'
        newName = 'Time';
    case 'Unknown_0040_A123'
        newName = 'PersonName';
    case 'Unknown_0040_A124'
        newName = 'UID';
    case 'Unknown_0040_A130'
        newName = 'TemporalRangeType';
    case 'Unknown_0040_A132'
        newName = 'ReferencedSamplePositions';
    case 'Unknown_0040_A136'
        newName = 'ReferencedFrameNumbers';
    case 'Unknown_0040_A138'
        newName = 'ReferencedTimeOffsets';
    case 'Unknown_0040_A13A'
        newName = 'ReferencedDateTime';
    case 'Unknown_0040_A160'
        newName = 'TextValue';
    case 'Unknown_0040_A168'
        newName = 'ConceptCodeSequence';
    case 'Unknown_0040_A170'
        newName = 'PurposeOfReferenceCodeSequence';
    case 'Unknown_0040_A180'
        newName = 'AnnotationGroupNumber';
    case 'Unknown_0040_A195'
        newName = 'ModifierCodeSequence';
    case 'Unknown_0040_A300'
        newName = 'MeasuredValueSequence';
    case 'Unknown_0040_A301'
        newName = 'NumericValueQualifierCodeSequence';
    case 'Unknown_0040_A30A'
        newName = 'NumericValue';
    case 'Unknown_0040_A353'
        newName = 'AddressTrial';
    case 'Unknown_0040_A354'
        newName = 'TelephoneNumberTrial';
    case 'Unknown_0040_A360'
        newName = 'PredecessorDocumentsSequence';
    case 'Unknown_0040_A370'
        newName = 'ReferencedRequestSequence';
    case 'Unknown_0040_A372'
        newName = 'PerformedProcedureCodeSequence';
    case 'Unknown_0040_A375'
        newName = 'CurrentRequestedProcedureEvidenceSequence';
    case 'Unknown_0040_A385'
        newName = 'PertinentOtherEvidenceSequence';
    case 'Unknown_0040_A390'
        newName = 'HL7StructuredDocumentReferenceSequence';
    case 'Unknown_0040_A491'
        newName = 'CompletionFlag';
    case 'Unknown_0040_A492'
        newName = 'CompletionFlagDescription';
    case 'Unknown_0040_A493'
        newName = 'VerificationFlag';
    case 'Unknown_0040_A504'
        newName = 'ContentTemplateSequence';
    case 'Unknown_0040_A525'
        newName = 'IdenticalDocumentsSequence';
    case 'Unknown_0040_A730'
        newName = 'ContentSequence';
    case 'Unknown_0040_B020'
        newName = 'AnnotationSequence';
    case 'Unknown_0040_DB00'
        newName = 'TemplateIdentifier';
    case 'Unknown_0040_DB06'
        newName = 'TemplateVersion';
    case 'Unknown_0040_DB07'
        newName = 'TemplateLocalVersion';
    case 'Unknown_0040_DB0B'
        newName = 'TemplateExtensionFlag';
    case 'Unknown_0040_DB0C'
        newName = 'TemplateExtensionOrganizationUID';
    case 'Unknown_0040_DB0D'
        newName = 'TemplateExtensionCreatorUID';
    case 'Unknown_0040_DB73'
        newName = 'ReferencedContentItemIdentifier';
    case 'Unknown_0040_E001'
        newName = 'HL7InstanceIdentifier';
    case 'Unknown_0040_E004'
        newName = 'HL7DocumentEffectiveTime';
    case 'Unknown_0040_E006'
        newName = 'HL7DocumentTypeCodeSequence';
    case 'Unknown_0040_E010'
        newName = 'RetrieveURI';
    case 'Unknown_0042_0010'
        newName = 'DocumentTitle';
    case 'Unknown_0042_0011'
        newName = 'EncapsulatedDocument';
    case 'Unknown_0042_0012'
        newName = 'MIMETypeOfEncapsulatedDocument';
    case 'Unknown_0042_0013'
        newName = 'SourceInstanceSequence';
    case 'Unknown_0050_0000'
        newName = 'CalibrationGroupLength';
    case 'Unknown_0050_0004'
        newName = 'CalibrationImage';
    case 'Unknown_0050_0010'
        newName = 'DeviceSequence';
    case 'Unknown_0050_0014'
        newName = 'DeviceLength';
    case 'Unknown_0050_0016'
        newName = 'DeviceDiameter';
    case 'Unknown_0050_0017'
        newName = 'DeviceDiameterUnits';
    case 'Unknown_0050_0018'
        newName = 'DeviceVolume';
    case 'Unknown_0050_0019'
        newName = 'InterMarkerDistance';
    case 'Unknown_0050_0020'
        newName = 'DeviceDescription';
    case 'Unknown_0054_0000'
        newName = 'NuclearAcquisitionGroupLength';
    case 'Unknown_0054_0010'
        newName = 'EnergyWindowVector';
    case 'Unknown_0054_0011'
        newName = 'NumberOfEnergyWindows';
    case 'Unknown_0054_0012'
        newName = 'EnergyWindowInformationSequence';
    case 'Unknown_0054_0013'
        newName = 'EnergyWindowRangeSequence';
    case 'Unknown_0054_0014'
        newName = 'EnergyWindowLowerLimit';
    case 'Unknown_0054_0015'
        newName = 'EnergyWindowUpperLimit';
    case 'Unknown_0054_0016'
        newName = 'RadiopharmaceuticalInformationSequence';
    case 'Unknown_0054_0017'
        newName = 'ResidualSyringeCounts';
    case 'Unknown_0054_0018'
        newName = 'EnergyWindowName';
    case 'Unknown_0054_0020'
        newName = 'DetectorVector';
    case 'Unknown_0054_0021'
        newName = 'NumberOfDetectors';
    case 'Unknown_0054_0022'
        newName = 'DetectorInformationSequence';
    case 'Unknown_0054_0030'
        newName = 'PhaseVector';
    case 'Unknown_0054_0031'
        newName = 'NumberOfPhases';
    case 'Unknown_0054_0032'
        newName = 'PhaseInformationSequence';
    case 'Unknown_0054_0033'
        newName = 'NumberOfFramesInPhase';
    case 'Unknown_0054_0036'
        newName = 'PhaseDelay';
    case 'Unknown_0054_0038'
        newName = 'PauseBetweenFrames';
    case 'Unknown_0054_0039'
        newName = 'PhaseDescription';
    case 'Unknown_0054_0050'
        newName = 'RotationVector';
    case 'Unknown_0054_0051'
        newName = 'NumberOfRotations';
    case 'Unknown_0054_0052'
        newName = 'RotationInformationSequence';
    case 'Unknown_0054_0053'
        newName = 'NumberOfFramesInRotation';
    case 'Unknown_0054_0060'
        newName = 'RRIntervalVector';
    case 'Unknown_0054_0061'
        newName = 'NumberOfRRIntervals';
    case 'Unknown_0054_0062'
        newName = 'GatedInformationSequence';
    case 'Unknown_0054_0063'
        newName = 'DataInformationSequence';
    case 'Unknown_0054_0070'
        newName = 'TimeSlotVector';
    case 'Unknown_0054_0071'
        newName = 'NumberOfTimeSlots';
    case 'Unknown_0054_0072'
        newName = 'TimeSlotInformationSequence';
    case 'Unknown_0054_0073'
        newName = 'TimeSlotTime';
    case 'Unknown_0054_0080'
        newName = 'SliceVector';
    case 'Unknown_0054_0081'
        newName = 'NumberOfSlices';
    case 'Unknown_0054_0090'
        newName = 'AngularViewVector';
    case 'Unknown_0054_0100'
        newName = 'TimeSliceVector';
    case 'Unknown_0054_0101'
        newName = 'NumberOfTimeSlices';
    case 'Unknown_0054_0200'
        newName = 'StartAngle';
    case 'Unknown_0054_0202'
        newName = 'TypeOfDetectorMotion';
    case 'Unknown_0054_0210'
        newName = 'TriggerVector';
    case 'Unknown_0054_0211'
        newName = 'NumberOfTriggersInPhase';
    case 'Unknown_0054_0220'
        newName = 'ViewCodeSequence';
    case 'Unknown_0054_0222'
        newName = 'ViewModifierCodeSequence';
    case 'Unknown_0054_0300'
        newName = 'RadionuclideCodeSequence';
    case 'Unknown_0054_0302'
        newName = 'AdministrationRouteCodeSequence';
    case 'Unknown_0054_0304'
        newName = 'RadiopharmaceuticalCodeSequence';
    case 'Unknown_0054_0306'
        newName = 'CalibrationDataSequence';
    case 'Unknown_0054_0308'
        newName = 'EnergyWindowNumber';
    case 'Unknown_0054_0400'
        newName = 'ImageID';
    case 'Unknown_0054_0410'
        newName = 'PatientOrientationCodeSequence';
    case 'Unknown_0054_0412'
        newName = 'PatientOrientationModifierCodeSequence';
    case 'Unknown_0054_0414'
        newName = 'PatientGantryRelationshipCodeSequence';
    case 'Unknown_0054_0500'
        newName = 'SliceProgressionDirection';
    case 'Unknown_0054_1000'
        newName = 'SeriesType';
    case 'Unknown_0054_1001'
        newName = 'Units';
    case 'Unknown_0054_1002'
        newName = 'CountsSource';
    case 'Unknown_0054_1004'
        newName = 'ReprojectionMethod';
    case 'Unknown_0054_1100'
        newName = 'RandomsCorrectionMethod';
    case 'Unknown_0054_1101'
        newName = 'AttenuationCorrectionMethod';
    case 'Unknown_0054_1102'
        newName = 'DecayCorrection';
    case 'Unknown_0054_1103'
        newName = 'ReconstructionMethod';
    case 'Unknown_0054_1104'
        newName = 'DetectorLinesOfResponseUsed';
    case 'Unknown_0054_1105'
        newName = 'ScatterCorrectionMethod';
    case 'Unknown_0054_1200'
        newName = 'AxialAcceptance';
    case 'Unknown_0054_1201'
        newName = 'AxialMash';
    case 'Unknown_0054_1202'
        newName = 'TransverseMash';
    case 'Unknown_0054_1203'
        newName = 'DetectorElementSize';
    case 'Unknown_0054_1210'
        newName = 'CoincidenceWindowWidth';
    case 'Unknown_0054_1220'
        newName = 'SecondaryCountsType';
    case 'Unknown_0054_1300'
        newName = 'FrameReferenceTime';
    case 'Unknown_0054_1310'
        newName = 'PrimaryPromptsCountsAccumulated';
    case 'Unknown_0054_1311'
        newName = 'SecondaryCountsAccumulated';
    case 'Unknown_0054_1320'
        newName = 'SliceSensitivityFactor';
    case 'Unknown_0054_1321'
        newName = 'DecayFactor';
    case 'Unknown_0054_1322'
        newName = 'DoseCalibrationFactor';
    case 'Unknown_0054_1323'
        newName = 'ScatterFractionFactor';
    case 'Unknown_0054_1324'
        newName = 'DeadTimeFactor';
    case 'Unknown_0054_1330'
        newName = 'ImageIndex';
    case 'Unknown_0054_1400'
        newName = 'CountsIncluded';
    case 'Unknown_0054_1401'
        newName = 'DeadTimeCorrectionFlag';
    case 'Unknown_0060_0000'
        newName = 'HistogramGroupLength';
    case 'Unknown_0060_3000'
        newName = 'HistogramSequence';
    case 'Unknown_0060_3002'
        newName = 'HistogramNumberOfBins';
    case 'Unknown_0060_3004'
        newName = 'HistogramFirstBinValue';
    case 'Unknown_0060_3006'
        newName = 'HistogramLastBinValue';
    case 'Unknown_0060_3008'
        newName = 'HistogramBinWidth';
    case 'Unknown_0060_3010'
        newName = 'HistogramExplanation';
    case 'Unknown_0060_3020'
        newName = 'HistogramData';
    case 'Unknown_0062_0001'
        newName = 'SegmentationType';
    case 'Unknown_0062_0002'
        newName = 'SegmentSequence';
    case 'Unknown_0062_0003'
        newName = 'SegmentedPropertyCategoryCodeSequence';
    case 'Unknown_0062_0004'
        newName = 'SegmentNumber';
    case 'Unknown_0062_0005'
        newName = 'SegmentLabel';
    case 'Unknown_0062_0006'
        newName = 'SegmentDescription';
    case 'Unknown_0062_0008'
        newName = 'SegmentAlgorithmType';
    case 'Unknown_0062_0009'
        newName = 'SegmentAlgorithmName';
    case 'Unknown_0062_000A'
        newName = 'SegmentIdentificationSequence';
    case 'Unknown_0062_000B'
        newName = 'ReferencedSegmentNumber';
    case 'Unknown_0062_000C'
        newName = 'RecommendedDisplayGrayscaleValue';
    case 'Unknown_0062_000D'
        newName = 'RecommendedDisplayCIELabValue';
    case 'Unknown_0062_000E'
        newName = 'MaximumFractionalValue';
    case 'Unknown_0062_000F'
        newName = 'SegmentedPropertyTypeCodeSequence';
    case 'Unknown_0062_0010'
        newName = 'SegmentationFractionalType';
    case 'Unknown_0064_0002'
        newName = 'DeformableRegistrationSequence';
    case 'Unknown_0064_0003'
        newName = 'SourceFrameOfReferenceUID';
    case 'Unknown_0064_0005'
        newName = 'DeformableRegistrationGridSequence';
    case 'Unknown_0064_0007'
        newName = 'GridDimensions';
    case 'Unknown_0064_0008'
        newName = 'GridResolution';
    case 'Unknown_0064_0009'
        newName = 'VectorGridData';
    case 'Unknown_0064_000F'
        newName = 'PreDeformationMatrixRegistrationSequence';
    case 'Unknown_0064_0010'
        newName = 'PostDeformationMatrixRegistrationSequence';
    case 'Unknown_0070_0001'
        newName = 'GraphicAnnotationSequence';
    case 'Unknown_0070_0002'
        newName = 'GraphicLayer';
    case 'Unknown_0070_0003'
        newName = 'BoundingBoxAnnotationUnits';
    case 'Unknown_0070_0004'
        newName = 'AnchorPointAnnotationUnits';
    case 'Unknown_0070_0005'
        newName = 'GraphicAnnotationUnits';
    case 'Unknown_0070_0006'
        newName = 'UnformattedTextValue';
    case 'Unknown_0070_0008'
        newName = 'TextObjectSequence';
    case 'Unknown_0070_0009'
        newName = 'GraphicObjectSequence';
    case 'Unknown_0070_0010'
        newName = 'BoundingBoxTLHC';
    case 'Unknown_0070_0011'
        newName = 'BoundingBoxBRHC';
    case 'Unknown_0070_0012'
        newName = 'BoundingBoxTextHorizontalJustification';
    case 'Unknown_0070_0014'
        newName = 'AnchorPoint';
    case 'Unknown_0070_0015'
        newName = 'AnchorPointVisibility';
    case 'Unknown_0070_0020'
        newName = 'GraphicDimensions';
    case 'Unknown_0070_0021'
        newName = 'NumberOfGraphicPoints';
    case 'Unknown_0070_0022'
        newName = 'GraphicData';
    case 'Unknown_0070_0023'
        newName = 'GraphicType';
    case 'Unknown_0070_0024'
        newName = 'GraphicFilled';
    case 'Unknown_0070_0040'
        newName = 'ImageRotationFrozenDraftRetired';
    case 'Unknown_0070_0041'
        newName = 'ImageHorizontalFlip';
    case 'Unknown_0070_0042'
        newName = 'ImageRotation';
    case 'Unknown_0070_0050'
        newName = 'DisplayedAreaTLHCFrozenDraftRetired';
    case 'Unknown_0070_0051'
        newName = 'DisplayedAreaBRHCFrozenDraftRetired';
    case 'Unknown_0070_0052'
        newName = 'DisplayedAreaTLHC';
    case 'Unknown_0070_0053'
        newName = 'DisplayedAreaBRHC';
    case 'Unknown_0070_005A'
        newName = 'DisplayedAreaSelectionSequence';
    case 'Unknown_0070_0060'
        newName = 'GraphicLayerSequence';
    case 'Unknown_0070_0062'
        newName = 'GraphicLayerOrder';
    case 'Unknown_0070_0066'
        newName = 'GraphicLayerRecommendedDisplayGrayscaleValue';
    case 'Unknown_0070_0067'
        newName = 'GraphicLayerRecommendedDisplayRGBValue';
    case 'Unknown_0070_0068'
        newName = 'GraphicLayerDescription';
    case 'Unknown_0070_0080'
        newName = 'ContentLabel';
    case 'Unknown_0070_0081'
        newName = 'ContentDescription';
    case 'Unknown_0070_0082'
        newName = 'PresentationCreationDate';
    case 'Unknown_0070_0083'
        newName = 'PresentationCreationTime';
    case 'Unknown_0070_0084'
        newName = 'ContentCreatorsName';
    case 'Unknown_0070_0086'
        newName = 'ContentCreatorsIdentificationCodeSequence';
    case 'Unknown_0070_0100'
        newName = 'PresentationSizeMode';
    case 'Unknown_0070_0101'
        newName = 'PresentationPixelSpacing';
    case 'Unknown_0070_0102'
        newName = 'PresentationPixelAspectRatio';
    case 'Unknown_0070_0103'
        newName = 'PresentationPixelMagnificationRatio';
    case 'Unknown_0070_0306'
        newName = 'ShapeType';
    case 'Unknown_0070_0308'
        newName = 'RegistrationSequence';
    case 'Unknown_0070_0309'
        newName = 'MatrixRegistrationSequence';
    case 'Unknown_0070_030A'
        newName = 'MatrixSequence';
    case 'Unknown_0070_030C'
        newName = 'FrameOfReferenceTransformationMatrixType';
    case 'Unknown_0070_030D'
        newName = 'RegistrationTypeCodeSequence';
    case 'Unknown_0070_030F'
        newName = 'FiducialDescription';
    case 'Unknown_0070_0310'
        newName = 'FiducialIdentifier';
    case 'Unknown_0070_0311'
        newName = 'FiducialIdentifierCodeSequence';
    case 'Unknown_0070_0312'
        newName = 'ContourUncertaintyRadius';
    case 'Unknown_0070_0314'
        newName = 'UsedFiducialsSequence';
    case 'Unknown_0070_0318'
        newName = 'GraphicCoordinatesDataSequence';
    case 'Unknown_0070_031A'
        newName = 'FiducialUID';
    case 'Unknown_0070_031C'
        newName = 'FiducialSetSequence';
    case 'Unknown_0070_031E'
        newName = 'FiducialSequence';
    case 'Unknown_0070_0401'
        newName = 'GraphicLayerRecommendedDisplayCIELabValue';
    case 'Unknown_0070_0402'
        newName = 'BlendingSequence';
    case 'Unknown_0070_0403'
        newName = 'RelativeOpacity';
    case 'Unknown_0070_0404'
        newName = 'ReferencedSpatialRegistrationSequence';
    case 'Unknown_0070_0405'
        newName = 'BlendingPosition';
    case 'Unknown_0072_0002'
        newName = 'HangingProtocolName';
    case 'Unknown_0072_0004'
        newName = 'HangingProtocolDescription';
    case 'Unknown_0072_0006'
        newName = 'HangingProtocolLevel';
    case 'Unknown_0072_0008'
        newName = 'HangingProtocolCreator';
    case 'Unknown_0072_000A'
        newName = 'HangingProtocolCreationDatetime';
    case 'Unknown_0072_000C'
        newName = 'HangingProtocolDefinitionSequence';
    case 'Unknown_0072_000E'
        newName = 'HangingProtocolUserIdentificationCodeSequence';
    case 'Unknown_0072_0010'
        newName = 'HangingProtocolUserGroupName';
    case 'Unknown_0072_0012'
        newName = 'SourceHangingProtocolSequence';
    case 'Unknown_0072_0014'
        newName = 'NumberOfPriorsReferenced';
    case 'Unknown_0072_0020'
        newName = 'ImageSetsSequence';
    case 'Unknown_0072_0022'
        newName = 'ImageSetSelectorSequence';
    case 'Unknown_0072_0024'
        newName = 'ImageSetSelectorUsageFlag';
    case 'Unknown_0072_0026'
        newName = 'SelectorAttribute';
    case 'Unknown_0072_0028'
        newName = 'SelectorValueNumber';
    case 'Unknown_0072_0030'
        newName = 'TimeBasedImageSetsSequence';
    case 'Unknown_0072_0032'
        newName = 'ImageSetNumber';
    case 'Unknown_0072_0034'
        newName = 'ImageSetSelectorCategory';
    case 'Unknown_0072_0038'
        newName = 'RelativeTime';
    case 'Unknown_0072_003A'
        newName = 'RelativeTimeUnits';
    case 'Unknown_0072_003C'
        newName = 'AbstractPriorValue';
    case 'Unknown_0072_003E'
        newName = 'AbstractPriorCodeSequence';
    case 'Unknown_0072_0040'
        newName = 'ImageSetLabel';
    case 'Unknown_0072_0050'
        newName = 'SelectorAttributeVR';
    case 'Unknown_0072_0052'
        newName = 'SelectorSequencePointer';
    case 'Unknown_0072_0054'
        newName = 'SelectorSequencePointerPrivateCreator';
    case 'Unknown_0072_0056'
        newName = 'SelectorAttributePrivateCreator';
    case 'Unknown_0072_0060'
        newName = 'SelectorATValue';
    case 'Unknown_0072_0062'
        newName = 'SelectorCSValue';
    case 'Unknown_0072_0064'
        newName = 'SelectorISValue';
    case 'Unknown_0072_0066'
        newName = 'SelectorLOValue';
    case 'Unknown_0072_0068'
        newName = 'SelectorLTValue';
    case 'Unknown_0072_006A'
        newName = 'SelectorPNValue';
    case 'Unknown_0072_006C'
        newName = 'SelectorSHValue';
    case 'Unknown_0072_006E'
        newName = 'SelectorSTValue';
    case 'Unknown_0072_0070'
        newName = 'SelectorUTValue';
    case 'Unknown_0072_0072'
        newName = 'SelectorDSValue';
    case 'Unknown_0072_0074'
        newName = 'SelectorFDValue';
    case 'Unknown_0072_0076'
        newName = 'SelectorFLValue';
    case 'Unknown_0072_0078'
        newName = 'SelectorULValue';
    case 'Unknown_0072_007A'
        newName = 'SelectorUSValue';
    case 'Unknown_0072_007C'
        newName = 'SelectorSLValue';
    case 'Unknown_0072_007E'
        newName = 'SelectorSSValue';
    case 'Unknown_0072_0080'
        newName = 'SelectorCodeSequenceValue';
    case 'Unknown_0072_0100'
        newName = 'NumberOfScreens';
    case 'Unknown_0072_0102'
        newName = 'NominalScreenDefinitionSequence';
    case 'Unknown_0072_0104'
        newName = 'NumberOfVerticalPixels';
    case 'Unknown_0072_0106'
        newName = 'NumberOfHorizontalPixels';
    case 'Unknown_0072_0108'
        newName = 'DisplayEnvironmentSpatialPosition';
    case 'Unknown_0072_010A'
        newName = 'ScreenMinimumGrayscaleBitDepth';
    case 'Unknown_0072_010C'
        newName = 'ScreenMinimumColorBitDepth';
    case 'Unknown_0072_010E'
        newName = 'ApplicationMaximumRepaintTime';
    case 'Unknown_0072_0200'
        newName = 'DisplaySetsSequence';
    case 'Unknown_0072_0202'
        newName = 'DisplaySetNumber';
    case 'Unknown_0072_0203'
        newName = 'DisplaySetLabel';
    case 'Unknown_0072_0204'
        newName = 'DisplaySetPresentationGroup';
    case 'Unknown_0072_0206'
        newName = 'DisplaySetPresentationGroupDescription';
    case 'Unknown_0072_0208'
        newName = 'PartialDataDisplayHandling';
    case 'Unknown_0072_0210'
        newName = 'SynchronizedScrollingSequence';
    case 'Unknown_0072_0212'
        newName = 'DisplaySetScrollingGroup';
    case 'Unknown_0072_0214'
        newName = 'NavigationIndicatorSequence';
    case 'Unknown_0072_0216'
        newName = 'NavigationDisplaySet';
    case 'Unknown_0072_0218'
        newName = 'ReferenceDisplaySets';
    case 'Unknown_0072_0300'
        newName = 'ImageBoxesSequence';
    case 'Unknown_0072_0302'
        newName = 'ImageBoxNumber';
    case 'Unknown_0072_0304'
        newName = 'ImageBoxLayoutType';
    case 'Unknown_0072_0306'
        newName = 'ImageBoxTileHorizontalDimension';
    case 'Unknown_0072_0308'
        newName = 'ImageBoxTileVerticalDimension';
    case 'Unknown_0072_0310'
        newName = 'ImageBoxScrollDirection';
    case 'Unknown_0072_0312'
        newName = 'ImageBoxSmallScrollType';
    case 'Unknown_0072_0314'
        newName = 'ImageBoxSmallScrollAmount';
    case 'Unknown_0072_0316'
        newName = 'ImageBoxLargeScrollType';
    case 'Unknown_0072_0318'
        newName = 'ImageBoxLargeScrollAmount';
    case 'Unknown_0072_0320'
        newName = 'ImageBoxOverlapPriority';
    case 'Unknown_0072_0330'
        newName = 'CineRelativeToRealTime';
    case 'Unknown_0072_0400'
        newName = 'FilterOperationsSequence';
    case 'Unknown_0072_0402'
        newName = 'FilterByCategory';
    case 'Unknown_0072_0404'
        newName = 'FilterByAttributePresence';
    case 'Unknown_0072_0406'
        newName = 'FilterByOperator';
    case 'Unknown_0072_0500'
        newName = 'BlendingOperationType';
    case 'Unknown_0072_0510'
        newName = 'ReformattingOperationType';
    case 'Unknown_0072_0512'
        newName = 'ReformattingThickness';
    case 'Unknown_0072_0514'
        newName = 'ReformattingInterval';
    case 'Unknown_0072_0516'
        newName = 'ReformattingOperationInitialViewDirection';
    case 'Unknown_0072_0520'
        newName = 'ThreeDRenderingType';
    case 'Unknown_0072_0600'
        newName = 'SortingOperationsSequence';
    case 'Unknown_0072_0602'
        newName = 'SortByCategory';
    case 'Unknown_0072_0604'
        newName = 'SortingDirection';
    case 'Unknown_0072_0700'
        newName = 'DisplaySetPatientOrientation';
    case 'Unknown_0072_0702'
        newName = 'VOIType';
    case 'Unknown_0072_0704'
        newName = 'PseudocolorType';
    case 'Unknown_0072_0706'
        newName = 'ShowGrayscaleInverted';
    case 'Unknown_0072_0710'
        newName = 'ShowImageTrueSizeFlag';
    case 'Unknown_0072_0712'
        newName = 'ShowGraphicAnnotationFlag';
    case 'Unknown_0072_0714'
        newName = 'ShowPatientDemographicsFlag';
    case 'Unknown_0072_0716'
        newName = 'ShowAcquisitionTechniquesFlag';
    case 'Unknown_0072_0717'
        newName = 'DisplaySetHorizontalJustification';
    case 'Unknown_0072_0718'
        newName = 'DisplaySetVerticalJustification';
    case 'Unknown_0088_0000'
        newName = 'StorageGroupLength';
    case 'Unknown_0088_0130'
        newName = 'StorageMediaFileSetID';
    case 'Unknown_0088_0140'
        newName = 'StorageMediaFileSetUID';
    case 'Unknown_0088_0200'
        newName = 'IconImageSequence';
    case 'Unknown_0088_0904'
        newName = 'TopicTitle';
    case 'Unknown_0088_0906'
        newName = 'TopicSubject';
    case 'Unknown_0088_0910'
        newName = 'TopicAuthor';
    case 'Unknown_0088_0912'
        newName = 'TopicKeyWords';
    case 'Unknown_0100_0410'
        newName = 'SOPInstanceStatus';
    case 'Unknown_0100_0420'
        newName = 'SOPAuthorizationDateAndTime';
    case 'Unknown_0100_0424'
        newName = 'SOPAuthorizationComment';
    case 'Unknown_0100_0426'
        newName = 'AuthorizationEquipmentCertificationNumber';
    case 'Unknown_0400_0005'
        newName = 'MACIDNumber';
    case 'Unknown_0400_0010'
        newName = 'MACCalculationTransferSyntaxUID';
    case 'Unknown_0400_0015'
        newName = 'MACAlgorithm';
    case 'Unknown_0400_0020'
        newName = 'DataElementsSigned';
    case 'Unknown_0400_0100'
        newName = 'DigitalSignatureUID';
    case 'Unknown_0400_0105'
        newName = 'DigitalSignatureDateTime';
    case 'Unknown_0400_0110'
        newName = 'CertificateType';
    case 'Unknown_0400_0115'
        newName = 'CertificateOfSigner';
    case 'Unknown_0400_0120'
        newName = 'Signature';
    case 'Unknown_0400_0305'
        newName = 'CertifiedTimestampType';
    case 'Unknown_0400_0310'
        newName = 'CertifiedTimestamp';
    case 'Unknown_0400_0401'
        newName = 'DigitalSignaturePurposeCodeSequence';
    case 'Unknown_0400_0402'
        newName = 'ReferencedDigitalSignatureSequence';
    case 'Unknown_0400_0403'
        newName = 'ReferencedSOPInstanceMACSequence';
    case 'Unknown_0400_0404'
        newName = 'MAC';
    case 'Unknown_0400_0500'
        newName = 'EncryptedAttributesSequence';
    case 'Unknown_0400_0510'
        newName = 'EncryptedContentTransferSyntaxUID';
    case 'Unknown_0400_0520'
        newName = 'EncryptedContent';
    case 'Unknown_0400_0550'
        newName = 'ModifiedAttributesSequence';
    case 'Unknown_0400_0561'
        newName = 'OriginalAttributesSequence';
    case 'Unknown_0400_0562'
        newName = 'AttributeModificationDatetime';
    case 'Unknown_0400_0563'
        newName = 'ModifyingSystem';
    case 'Unknown_0400_0564'
        newName = 'SourceOfPreviousValues';
    case 'Unknown_0400_0565'
        newName = 'ReasonForTheAttributeModification';
    case 'Unknown_1000_0000'
        newName = 'CodeTableGroupLength';
    case 'Unknown_1000_0011'
        newName = 'RunLengthTriplet';
    case 'Unknown_1000_0012'
        newName = 'HuffmanTableSize';
    case 'Unknown_1000_0013'
        newName = 'HuffmanTableTriplet';
    case 'Unknown_1000_0014'
        newName = 'ShiftTableSize';
    case 'Unknown_1000_0015'
        newName = 'ShiftTableTriplet';
    case 'Unknown_1010_0000'
        newName = 'ZonalMapGroupLength';
    case 'Unknown_1010_0004'
        newName = 'ZonalMap';
    case 'Unknown_2000_0000'
        newName = 'FilmSessionGroupLength';
    case 'Unknown_2000_0010'
        newName = 'NumberOfCopies';
    case 'Unknown_2000_001E'
        newName = 'PrinterConfigurationSequence';
    case 'Unknown_2000_0020'
        newName = 'PrintPriority';
    case 'Unknown_2000_0030'
        newName = 'MediumType';
    case 'Unknown_2000_0040'
        newName = 'FilmDestination';
    case 'Unknown_2000_0050'
        newName = 'FilmSessionLabel';
    case 'Unknown_2000_0060'
        newName = 'MemoryAllocation';
    case 'Unknown_2000_0061'
        newName = 'MaximumMemoryAllocation';
    case 'Unknown_2000_0062'
        newName = 'ColorImagePrintingFlag';
    case 'Unknown_2000_0063'
        newName = 'CollationFlag';
    case 'Unknown_2000_0065'
        newName = 'AnnotationFlag';
    case 'Unknown_2000_0067'
        newName = 'ImageOverlayFlag';
    case 'Unknown_2000_0069'
        newName = 'PresentationLUTFlag';
    case 'Unknown_2000_006A'
        newName = 'ImageBoxPresentationLUTFlag';
    case 'Unknown_2000_00A0'
        newName = 'MemoryBitDepth';
    case 'Unknown_2000_00A1'
        newName = 'PrintingBitDepth';
    case 'Unknown_2000_00A2'
        newName = 'MediaInstalledSequence';
    case 'Unknown_2000_00A4'
        newName = 'OtherMediaAvailableSequence';
    case 'Unknown_2000_00A8'
        newName = 'SupportedImageDisplayFormatsSequence';
    case 'Unknown_2000_0500'
        newName = 'ReferencedFilmBoxSequence';
    case 'Unknown_2000_0510'
        newName = 'ReferencedStoredPrintSequence';
    case 'Unknown_2010_0000'
        newName = 'FilmBoxGroupLength';
    case 'Unknown_2010_0010'
        newName = 'ImageDisplayFormat';
    case 'Unknown_2010_0030'
        newName = 'AnnotationDisplayFormatID';
    case 'Unknown_2010_0040'
        newName = 'FilmOrientation';
    case 'Unknown_2010_0050'
        newName = 'FilmSizeID';
    case 'Unknown_2010_0052'
        newName = 'PrinterResolutionID';
    case 'Unknown_2010_0054'
        newName = 'DefaultPrinterResolutionID';
    case 'Unknown_2010_0060'
        newName = 'MagnificationType';
    case 'Unknown_2010_0080'
        newName = 'SmoothingType';
    case 'Unknown_2010_00A6'
        newName = 'DefaultMagnificationType';
    case 'Unknown_2010_00A7'
        newName = 'OtherMagnificationTypesAvailable';
    case 'Unknown_2010_00A8'
        newName = 'DefaultSmoothingType';
    case 'Unknown_2010_00A9'
        newName = 'OtherSmoothingTypesAvailable';
    case 'Unknown_2010_0100'
        newName = 'BorderDensity';
    case 'Unknown_2010_0110'
        newName = 'EmptyImageDensity';
    case 'Unknown_2010_0120'
        newName = 'MinDensity';
    case 'Unknown_2010_0130'
        newName = 'MaxDensity';
    case 'Unknown_2010_0140'
        newName = 'Trim';
    case 'Unknown_2010_0150'
        newName = 'ConfigurationInformation';
    case 'Unknown_2010_0152'
        newName = 'ConfigurationInformationDescription';
    case 'Unknown_2010_0154'
        newName = 'MaximumCollatedFilms';
    case 'Unknown_2010_015E'
        newName = 'Illumination';
    case 'Unknown_2010_0160'
        newName = 'ReflectedAmbientLight';
    case 'Unknown_2010_0376'
        newName = 'PrinterPixelSpacing';
    case 'Unknown_2010_0500'
        newName = 'ReferencedFilmSessionSequence';
    case 'Unknown_2010_0510'
        newName = 'ReferencedImageBoxSequence';
    case 'Unknown_2010_0520'
        newName = 'ReferencedBasicAnnotationBoxSequence';
    case 'Unknown_2020_0000'
        newName = 'ImageBoxGroupLength';
    case 'Unknown_2020_0010'
        newName = 'ImageBoxPosition';
    case 'Unknown_2020_0020'
        newName = 'Polarity';
    case 'Unknown_2020_0030'
        newName = 'RequestedImageSize';
    case 'Unknown_2020_0040'
        newName = 'RequestedDecimateCropBehavior';
    case 'Unknown_2020_0050'
        newName = 'RequestedResolutionID';
    case 'Unknown_2020_00A0'
        newName = 'RequestedImageSizeFlag';
    case 'Unknown_2020_00A2'
        newName = 'DecimateCropResult';
    case 'Unknown_2020_0110'
        newName = 'BasicGrayscaleImageSequence';
    case 'Unknown_2020_0111'
        newName = 'BasicColorImageSequence';
    case 'Unknown_2020_0130'
        newName = 'ReferencedImageOverlayBoxSequence';
    case 'Unknown_2020_0140'
        newName = 'ReferencedVOILUTBoxSequence';
    case 'Unknown_2030_0000'
        newName = 'AnnotationGroupLength';
    case 'Unknown_2030_0010'
        newName = 'AnnotationPosition';
    case 'Unknown_2030_0020'
        newName = 'TextString';
    case 'Unknown_2040_0000'
        newName = 'OverlayBoxGroupLength';
    case 'Unknown_2040_0010'
        newName = 'ReferencedOverlayPlaneSequence';
    case 'Unknown_2040_0011'
        newName = 'ReferencedOverlayPlaneGroups';
    case 'Unknown_2040_0020'
        newName = 'OverlayPixelDataSequence';
    case 'Unknown_2040_0060'
        newName = 'OverlayMagnificationType';
    case 'Unknown_2040_0070'
        newName = 'OverlaySmoothingType';
    case 'Unknown_2040_0072'
        newName = 'OverlayOrImageMagnification';
    case 'Unknown_2040_0074'
        newName = 'MagnifyToNumberOfColumns';
    case 'Unknown_2040_0080'
        newName = 'OverlayForegroundDensity';
    case 'Unknown_2040_0082'
        newName = 'OverlayBackgroundDensity';
    case 'Unknown_2040_0090'
        newName = 'OverlayMode';
    case 'Unknown_2040_0100'
        newName = 'ThresholdDensity';
    case 'Unknown_2040_0500'
        newName = 'ReferencedOverlayImageBoxSequence';
    case 'Unknown_2050_0010'
        newName = 'PresentationLUTSequence';
    case 'Unknown_2050_0020'
        newName = 'PresentationLUTShape';
    case 'Unknown_2050_0500'
        newName = 'ReferencedPresentationLUTSequence';
    case 'Unknown_2100_0000'
        newName = 'PrintJobGroupLength';
    case 'Unknown_2100_0010'
        newName = 'PrintJobID';
    case 'Unknown_2100_0020'
        newName = 'ExecutionStatus';
    case 'Unknown_2100_0030'
        newName = 'ExecutionStatusInfo';
    case 'Unknown_2100_0040'
        newName = 'CreationDate';
    case 'Unknown_2100_0050'
        newName = 'CreationTime';
    case 'Unknown_2100_0070'
        newName = 'Originator';
    case 'Unknown_2100_0140'
        newName = 'DestinationAE';
    case 'Unknown_2100_0160'
        newName = 'OwnerID';
    case 'Unknown_2100_0170'
        newName = 'NumberOfFilms';
    case 'Unknown_2100_0500'
        newName = 'ReferencedPrintJobSequencePull';
    case 'Unknown_2110_0000'
        newName = 'PrinterGroupLength';
    case 'Unknown_2110_0010'
        newName = 'PrinterStatus';
    case 'Unknown_2110_0020'
        newName = 'PrinterStatusInfo';
    case 'Unknown_2110_0030'
        newName = 'PrinterName';
    case 'Unknown_2110_0099'
        newName = 'PrintQueueID';
    case 'Unknown_2120_0010'
        newName = 'QueueStatus';
    case 'Unknown_2120_0050'
        newName = 'PrintJobDescriptionSequence';
    case 'Unknown_2120_0070'
        newName = 'ReferencedPrintJobSequenceQueue';
    case 'Unknown_2130_0010'
        newName = 'PrintManagementCapabilitiesSequence';
    case 'Unknown_2130_0015'
        newName = 'PrinterCharacteristicsSequence';
    case 'Unknown_2130_0030'
        newName = 'FilmBoxContentSequence';
    case 'Unknown_2130_0040'
        newName = 'ImageBoxContentSequence';
    case 'Unknown_2130_0050'
        newName = 'AnnotationContentSequence';
    case 'Unknown_2130_0060'
        newName = 'ImageOverlayBoxContentSequence';
    case 'Unknown_2130_0080'
        newName = 'PresentationLUTContentSequence';
    case 'Unknown_2130_00A0'
        newName = 'ProposedStudySequence';
    case 'Unknown_2130_00C0'
        newName = 'OriginalImageSequence';
    case 'Unknown_2200_0001'
        newName = 'LabelUsingInformationExtractedFromInstances';
    case 'Unknown_2200_0002'
        newName = 'LabelText';
    case 'Unknown_2200_0003'
        newName = 'LabelStyleSelection';
    case 'Unknown_2200_0004'
        newName = 'MediaDisposition';
    case 'Unknown_2200_0005'
        newName = 'BarcodeValue';
    case 'Unknown_2200_0006'
        newName = 'BarcodeSymbology';
    case 'Unknown_2200_0007'
        newName = 'AllowMediaSplitting';
    case 'Unknown_2200_0008'
        newName = 'IncludeNonDICOMObjects';
    case 'Unknown_2200_0009'
        newName = 'IncludeDisplayApplication';
    case 'Unknown_2200_000A'
        newName = 'PreserveCompositeInstancesAfterMediaCreation';
    case 'Unknown_2200_000B'
        newName = 'TotalNumberOfPiecesOfMediaCreated';
    case 'Unknown_2200_000C'
        newName = 'RequestedMediaApplicationProfile';
    case 'Unknown_2200_000D'
        newName = 'ReferencedStorageMediaSequence';
    case 'Unknown_2200_000E'
        newName = 'FailureAttributes';
    case 'Unknown_2200_000F'
        newName = 'AllowLossyCompression';
    case 'Unknown_2200_0020'
        newName = 'RequestPriority';
    case 'Unknown_3002_0000'
        newName = 'RTGroupLength';
    case 'Unknown_3002_0002'
        newName = 'RTImageLabel';
    case 'Unknown_3002_0003'
        newName = 'RTImageName';
    case 'Unknown_3002_0004'
        newName = 'RTImageDescription';
    case 'Unknown_3002_000A'
        newName = 'ReportedValuesOrigin';
    case 'Unknown_3002_000C'
        newName = 'RTImagePlane';
    case 'Unknown_3002_000D'
        newName = 'XRayImageReceptorTranslation';
    case 'Unknown_3002_000E'
        newName = 'XRayImageReceptorAngle';
    case 'Unknown_3002_0010'
        newName = 'RTImageOrientation';
    case 'Unknown_3002_0011'
        newName = 'ImagePlanePixelSpacing';
    case 'Unknown_3002_0012'
        newName = 'RTImagePosition';
    case 'Unknown_3002_0020'
        newName = 'RadiationMachineName';
    case 'Unknown_3002_0022'
        newName = 'RadiationMachineSAD';
    case 'Unknown_3002_0024'
        newName = 'RadiationMachineSSD';
    case 'Unknown_3002_0026'
        newName = 'RTImageSID';
    case 'Unknown_3002_0028'
        newName = 'SourceToReferenceObjectDistance';
    case 'Unknown_3002_0029'
        newName = 'FractionNumber';
    case 'Unknown_3002_0030'
        newName = 'ExposureSequence';
    case 'Unknown_3002_0032'
        newName = 'MetersetExposure';
    case 'Unknown_3002_0034'
        newName = 'DiaphragmPosition';
    case 'Unknown_3002_0040'
        newName = 'FluenceMapSequence';
    case 'Unknown_3002_0041'
        newName = 'FluenceDataSource';
    case 'Unknown_3002_0042'
        newName = 'FluenceDataScale';
    case 'Unknown_3004_0000'
        newName = 'DoseGroupLength';
    case 'Unknown_3004_0001'
        newName = 'DVHType';
    case 'Unknown_3004_0002'
        newName = 'DoseUnits';
    case 'Unknown_3004_0004'
        newName = 'DoseType';
    case 'Unknown_3004_0006'
        newName = 'DoseComment';
    case 'Unknown_3004_0008'
        newName = 'NormalizationPoint';
    case 'Unknown_3004_000A'
        newName = 'DoseSummationType';
    case 'Unknown_3004_000C'
        newName = 'GridFrameOffsetVector';
    case 'Unknown_3004_000E'
        newName = 'DoseGridScaling';
    case 'Unknown_3004_0010'
        newName = 'RTDoseROISequence';
    case 'Unknown_3004_0012'
        newName = 'DoseValue';
    case 'Unknown_3004_0014'
        newName = 'TissueHeterogeneityCorrection';
    case 'Unknown_3004_0040'
        newName = 'DVHNormalizationPoint';
    case 'Unknown_3004_0042'
        newName = 'DVHNormalizationDoseValue';
    case 'Unknown_3004_0050'
        newName = 'DVHSequence';
    case 'Unknown_3004_0052'
        newName = 'DVHDoseScaling';
    case 'Unknown_3004_0054'
        newName = 'DVHVolumeUnits';
    case 'Unknown_3004_0056'
        newName = 'DVHNumberOfBins';
    case 'Unknown_3004_0058'
        newName = 'DVHData';
    case 'Unknown_3004_0060'
        newName = 'DVHReferencedROISequence';
    case 'Unknown_3004_0062'
        newName = 'DVHROIContributionType';
    case 'Unknown_3004_0070'
        newName = 'DVHMinimumDose';
    case 'Unknown_3004_0072'
        newName = 'DVHMaximumDose';
    case 'Unknown_3004_0074'
        newName = 'DVHMeanDose';
    case 'Unknown_3006_0000'
        newName = 'ROIGroupLength';
    case 'Unknown_3006_0002'
        newName = 'StructureSetLabel';
    case 'Unknown_3006_0004'
        newName = 'StructureSetName';
    case 'Unknown_3006_0006'
        newName = 'StructureSetDescription';
    case 'Unknown_3006_0008'
        newName = 'StructureSetDate';
    case 'Unknown_3006_0009'
        newName = 'StructureSetTime';
    case 'Unknown_3006_0010'
        newName = 'ReferencedFrameOfReferenceSequence';
    case 'Unknown_3006_0012'
        newName = 'RTReferencedStudySequence';
    case 'Unknown_3006_0014'
        newName = 'RTReferencedSeriesSequence';
    case 'Unknown_3006_0016'
        newName = 'ContourImageSequence';
    case 'Unknown_3006_0020'
        newName = 'StructureSetROISequence';
    case 'Unknown_3006_0022'
        newName = 'ROINumber';
    case 'Unknown_3006_0024'
        newName = 'ReferencedFrameOfReferenceUID';
    case 'Unknown_3006_0026'
        newName = 'ROIName';
    case 'Unknown_3006_0028'
        newName = 'ROIDescription';
    case 'Unknown_3006_002A'
        newName = 'ROIDisplayColor';
    case 'Unknown_3006_002C'
        newName = 'ROIVolume';
    case 'Unknown_3006_0030'
        newName = 'RTRelatedROISequence';
    case 'Unknown_3006_0033'
        newName = 'RTROIRelationship';
    case 'Unknown_3006_0036'
        newName = 'ROIGenerationAlgorithm';
    case 'Unknown_3006_0038'
        newName = 'ROIGenerationDescription';
    case 'Unknown_3006_0039'
        newName = 'ROIContourSequence';
    case 'Unknown_3006_0040'
        newName = 'ContourSequence';
    case 'Unknown_3006_0042'
        newName = 'ContourGeometricType';
    case 'Unknown_3006_0044'
        newName = 'ContourSlabThickness';
    case 'Unknown_3006_0045'
        newName = 'ContourOffsetVector';
    case 'Unknown_3006_0046'
        newName = 'NumberOfContourPoints';
    case 'Unknown_3006_0048'
        newName = 'ContourNumber';
    case 'Unknown_3006_0049'
        newName = 'AttachedContours';
    case 'Unknown_3006_0050'
        newName = 'ContourData';
    case 'Unknown_3006_0080'
        newName = 'RTROIObservationsSequence';
    case 'Unknown_3006_0082'
        newName = 'ObservationNumber';
    case 'Unknown_3006_0084'
        newName = 'ReferencedROINumber';
    case 'Unknown_3006_0085'
        newName = 'ROIObservationLabel';
    case 'Unknown_3006_0086'
        newName = 'RTROIIdentificationCodeSequence';
    case 'Unknown_3006_0088'
        newName = 'ROIObservationDescription';
    case 'Unknown_3006_00A0'
        newName = 'RelatedRTROIObservationsSequence';
    case 'Unknown_3006_00A4'
        newName = 'RTROIInterpretedType';
    case 'Unknown_3006_00A6'
        newName = 'ROIInterpreter';
    case 'Unknown_3006_00B0'
        newName = 'ROIPhysicalPropertiesSequence';
    case 'Unknown_3006_00B2'
        newName = 'ROIPhysicalProperty';
    case 'Unknown_3006_00B4'
        newName = 'ROIPhysicalPropertyValue';
    case 'Unknown_3006_00C0'
        newName = 'FrameOfReferenceRelationshipSequence';
    case 'Unknown_3006_00C2'
        newName = 'RelatedFrameOfReferenceUID';
    case 'Unknown_3006_00C4'
        newName = 'FrameOfReferenceTransformationType';
    case 'Unknown_3006_00C6'
        newName = 'FrameOfReferenceTransformationMatrix';
    case 'Unknown_3006_00C8'
        newName = 'FrameOfReferenceTransformationComment';
    case 'Unknown_3008_0000'
        newName = 'TreatmentGroupLength';
    case 'Unknown_3008_0010'
        newName = 'MeasuredDoseReferenceSequence';
    case 'Unknown_3008_0012'
        newName = 'MeasuredDoseDescription';
    case 'Unknown_3008_0014'
        newName = 'MeasuredDoseType';
    case 'Unknown_3008_0016'
        newName = 'MeasuredDoseValue';
    case 'Unknown_3008_0020'
        newName = 'TreatmentSessionBeamSequence';
    case 'Unknown_3008_0021'
        newName = 'TreatmentSessionIonBeamSequence';
    case 'Unknown_3008_0022'
        newName = 'CurrentFractionNumber';
    case 'Unknown_3008_0024'
        newName = 'TreatmentControlPointDate';
    case 'Unknown_3008_0025'
        newName = 'TreatmentControlPointTime';
    case 'Unknown_3008_002A'
        newName = 'TreatmentTerminationStatus';
    case 'Unknown_3008_002B'
        newName = 'TreatmentTerminationCode';
    case 'Unknown_3008_002C'
        newName = 'TreatmentVerificationStatus';
    case 'Unknown_3008_0030'
        newName = 'ReferencedTreatmentRecordSequence';
    case 'Unknown_3008_0032'
        newName = 'SpecifiedPrimaryMeterset';
    case 'Unknown_3008_0033'
        newName = 'SpecifiedSecondaryMeterset';
    case 'Unknown_3008_0036'
        newName = 'DeliveredPrimaryMeterset';
    case 'Unknown_3008_0037'
        newName = 'DeliveredSecondaryMeterset';
    case 'Unknown_3008_003A'
        newName = 'SpecifiedTreatmentTime';
    case 'Unknown_3008_003B'
        newName = 'DeliveredTreatmentTime';
    case 'Unknown_3008_0040'
        newName = 'ControlPointDeliverySequence';
    case 'Unknown_3008_0041'
        newName = 'IonControlPointDeliverySequence';
    case 'Unknown_3008_0042'
        newName = 'SpecifiedMeterset';
    case 'Unknown_3008_0044'
        newName = 'DeliveredMeterset';
    case 'Unknown_3008_0045'
        newName = 'MetersetRateSet';
    case 'Unknown_3008_0046'
        newName = 'MetersetRateDelivered';
    case 'Unknown_3008_0047'
        newName = 'ScanSpotMetersetsDelivered';
    case 'Unknown_3008_0048'
        newName = 'DoseRateDelivered';
    case 'Unknown_3008_0050'
        newName = 'TreatmentSummaryCalculatedDoseReferenceSequence';
    case 'Unknown_3008_0052'
        newName = 'CumulativeDoseToDoseReference';
    case 'Unknown_3008_0054'
        newName = 'FirstTreatmentDate';
    case 'Unknown_3008_0056'
        newName = 'MostRecentTreatmentDate';
    case 'Unknown_3008_005A'
        newName = 'NumberOfFractionsDelivered';
    case 'Unknown_3008_0060'
        newName = 'OverrideSequence';
    case 'Unknown_3008_0061'
        newName = 'ParameterSequencePointer';
    case 'Unknown_3008_0062'
        newName = 'OverrideParameterPointer';
    case 'Unknown_3008_0063'
        newName = 'ParameterItemIndex';
    case 'Unknown_3008_0064'
        newName = 'MeasuredDoseReferenceNumber';
    case 'Unknown_3008_0065'
        newName = 'ParameterPointer';
    case 'Unknown_3008_0066'
        newName = 'OverrideReason';
    case 'Unknown_3008_0068'
        newName = 'CorrectedParameterSequence';
    case 'Unknown_3008_006A'
        newName = 'CorrectionValue';
    case 'Unknown_3008_0070'
        newName = 'CalculatedDoseReferenceSequence';
    case 'Unknown_3008_0072'
        newName = 'CalculatedDoseReferenceNumber';
    case 'Unknown_3008_0074'
        newName = 'CalculatedDoseReferenceDescription';
    case 'Unknown_3008_0076'
        newName = 'CalculatedDoseReferenceDoseValue';
    case 'Unknown_3008_0078'
        newName = 'StartMeterset';
    case 'Unknown_3008_007A'
        newName = 'EndMeterset';
    case 'Unknown_3008_0080'
        newName = 'ReferencedMeasuredDoseReferenceSequence';
    case 'Unknown_3008_0082'
        newName = 'ReferencedMeasuredDoseReferenceNumber';
    case 'Unknown_3008_0090'
        newName = 'ReferencedCalculatedDoseReferenceSequence';
    case 'Unknown_3008_0092'
        newName = 'ReferencedCalculatedDoseReferenceNumber';
    case 'Unknown_3008_00A0'
        newName = 'BeamLimitingDeviceLeafPairsSequence';
    case 'Unknown_3008_00B0'
        newName = 'RecordedWedgeSequence';
    case 'Unknown_3008_00C0'
        newName = 'RecordedCompensatorSequence';
    case 'Unknown_3008_00D0'
        newName = 'RecordedBlockSequence';
    case 'Unknown_3008_00E0'
        newName = 'TreatmentSummaryMeasuredDoseReferenceSequence';
    case 'Unknown_3008_00F0'
        newName = 'RecordedSnoutSequence';
    case 'Unknown_3008_00F2'
        newName = 'RecordedRangeShifterSequence';
    case 'Unknown_3008_00F4'
        newName = 'RecordedLateralSpreadingDeviceSequence';
    case 'Unknown_3008_00F6'
        newName = 'RecordedRangeModulatorSequence';
    case 'Unknown_3008_0100'
        newName = 'RecordedSourceSequence';
    case 'Unknown_3008_0105'
        newName = 'SourceSerialNumber';
    case 'Unknown_3008_0110'
        newName = 'TreatmentSessionApplicationSetupSequence';
    case 'Unknown_3008_0116'
        newName = 'ApplicationSetupCheck';
    case 'Unknown_3008_0120'
        newName = 'RecordedBrachyAccessoryDeviceSequence';
    case 'Unknown_3008_0122'
        newName = 'ReferencedBrachyAccessoryDeviceNumber';
    case 'Unknown_3008_0130'
        newName = 'RecordedChannelSequence';
    case 'Unknown_3008_0132'
        newName = 'SpecifiedChannelTotalTime';
    case 'Unknown_3008_0134'
        newName = 'DeliveredChannelTotalTime';
    case 'Unknown_3008_0136'
        newName = 'SpecifiedNumberOfPulses';
    case 'Unknown_3008_0138'
        newName = 'DeliveredNumberOfPulses';
    case 'Unknown_3008_013A'
        newName = 'SpecifiedPulseRepetitionInterval';
    case 'Unknown_3008_013C'
        newName = 'DeliveredPulseRepetitionInterval';
    case 'Unknown_3008_0140'
        newName = 'RecordedSourceApplicatorSequence';
    case 'Unknown_3008_0142'
        newName = 'ReferencedSourceApplicatorNumber';
    case 'Unknown_3008_0150'
        newName = 'RecordedChannelShieldSequence';
    case 'Unknown_3008_0152'
        newName = 'ReferencedChannelShieldNumber';
    case 'Unknown_3008_0160'
        newName = 'BrachyControlPointDeliveredSequence';
    case 'Unknown_3008_0162'
        newName = 'SafePositionExitDate';
    case 'Unknown_3008_0164'
        newName = 'SafePositionExitTime';
    case 'Unknown_3008_0166'
        newName = 'SafePositionReturnDate';
    case 'Unknown_3008_0168'
        newName = 'SafePositionReturnTime';
    case 'Unknown_3008_0200'
        newName = 'CurrentTreatmentStatus';
    case 'Unknown_3008_0202'
        newName = 'TreatmentStatusComment';
    case 'Unknown_3008_0220'
        newName = 'FractionGroupSummarySequence';
    case 'Unknown_3008_0223'
        newName = 'ReferencedFractionNumber';
    case 'Unknown_3008_0224'
        newName = 'FractionGroupType';
    case 'Unknown_3008_0230'
        newName = 'BeamStopperPosition';
    case 'Unknown_3008_0240'
        newName = 'FractionStatusSummarySequence';
    case 'Unknown_3008_0250'
        newName = 'TreatmentDate';
    case 'Unknown_3008_0251'
        newName = 'TreatmentTime';
    case 'Unknown_300A_0000'
        newName = 'PlanGroupLength';
    case 'Unknown_300A_0002'
        newName = 'RTPlanLabel';
    case 'Unknown_300A_0003'
        newName = 'RTPlanName';
    case 'Unknown_300A_0004'
        newName = 'RTPlanDescription';
    case 'Unknown_300A_0006'
        newName = 'RTPlanDate';
    case 'Unknown_300A_0007'
        newName = 'RTPlanTime';
    case 'Unknown_300A_0009'
        newName = 'TreatmentProtocols';
    case 'Unknown_300A_000A'
        newName = 'PlanIntent';
    case 'Unknown_300A_000B'
        newName = 'TreatmentSites';
    case 'Unknown_300A_000C'
        newName = 'RTPlanGeometry';
    case 'Unknown_300A_000E'
        newName = 'PrescriptionDescription';
    case 'Unknown_300A_0010'
        newName = 'DoseReferenceSequence';
    case 'Unknown_300A_0012'
        newName = 'DoseReferenceNumber';
    case 'Unknown_300A_0013'
        newName = 'DoseReferenceUID';
    case 'Unknown_300A_0014'
        newName = 'DoseReferenceStructureType';
    case 'Unknown_300A_0015'
        newName = 'NominalBeamEnergyUnit';
    case 'Unknown_300A_0016'
        newName = 'DoseReferenceDescription';
    case 'Unknown_300A_0018'
        newName = 'DoseReferencePointCoordinates';
    case 'Unknown_300A_001A'
        newName = 'NominalPriorDose';
    case 'Unknown_300A_0020'
        newName = 'DoseReferenceType';
    case 'Unknown_300A_0021'
        newName = 'ConstraintWeight';
    case 'Unknown_300A_0022'
        newName = 'DeliveryWarningDose';
    case 'Unknown_300A_0023'
        newName = 'DeliveryMaximumDose';
    case 'Unknown_300A_0025'
        newName = 'TargetMinimumDose';
    case 'Unknown_300A_0026'
        newName = 'TargetPrescriptionDose';
    case 'Unknown_300A_0027'
        newName = 'TargetMaximumDose';
    case 'Unknown_300A_0028'
        newName = 'TargetUnderdoseVolumeFraction';
    case 'Unknown_300A_002A'
        newName = 'OrganAtRiskFullVolumeDose';
    case 'Unknown_300A_002B'
        newName = 'OrganAtRiskLimitDose';
    case 'Unknown_300A_002C'
        newName = 'OrganAtRiskMaximumDose';
    case 'Unknown_300A_002D'
        newName = 'OrganAtRiskOverdoseVolumeFraction';
    case 'Unknown_300A_0040'
        newName = 'ToleranceTableSequence';
    case 'Unknown_300A_0042'
        newName = 'ToleranceTableNumber';
    case 'Unknown_300A_0043'
        newName = 'ToleranceTableLabel';
    case 'Unknown_300A_0044'
        newName = 'GantryAngleTolerance';
    case 'Unknown_300A_0046'
        newName = 'BeamLimitingDeviceAngleTolerance';
    case 'Unknown_300A_0048'
        newName = 'BeamLimitingDeviceToleranceSequence';
    case 'Unknown_300A_004A'
        newName = 'BeamLimitingDevicePositionTolerance';
    case 'Unknown_300A_004B'
        newName = 'SnoutPositionTolerance';
    case 'Unknown_300A_004C'
        newName = 'PatientSupportAngleTolerance';
    case 'Unknown_300A_004E'
        newName = 'TableTopEccentricAngleTolerance';
    case 'Unknown_300A_004F'
        newName = 'TableTopPitchAngleTolerance';
    case 'Unknown_300A_0050'
        newName = 'TableTopRollAngleTolerance';
    case 'Unknown_300A_0051'
        newName = 'TableTopVerticalPositionTolerance';
    case 'Unknown_300A_0052'
        newName = 'TableTopLongitudinalPositionTolerance';
    case 'Unknown_300A_0053'
        newName = 'TableTopLateralPositionTolerance';
    case 'Unknown_300A_0055'
        newName = 'RTPlanRelationship';
    case 'Unknown_300A_0070'
        newName = 'FractionGroupSequence';
    case 'Unknown_300A_0071'
        newName = 'FractionGroupNumber';
    case 'Unknown_300A_0072'
        newName = 'FractionGroupDescription';
    case 'Unknown_300A_0078'
        newName = 'NumberOfFractionsPlanned';
    case 'Unknown_300A_0079'
        newName = 'NumberOfFractionPatternDigitsPerDay';
    case 'Unknown_300A_007A'
        newName = 'RepeatFractionCycleLength';
    case 'Unknown_300A_007B'
        newName = 'FractionPattern';
    case 'Unknown_300A_0080'
        newName = 'NumberOfBeams';
    case 'Unknown_300A_0082'
        newName = 'BeamDoseSpecificationPoint';
    case 'Unknown_300A_0084'
        newName = 'BeamDose';
    case 'Unknown_300A_0086'
        newName = 'BeamMeterset';
    case 'Unknown_300A_0088'
        newName = 'BeamDosePointDepth';
    case 'Unknown_300A_0089'
        newName = 'BeamDosePointEquivalentDepth';
    case 'Unknown_300A_008A'
        newName = 'BeamDosePointSSD';
    case 'Unknown_300A_00A0'
        newName = 'NumberOfBrachyApplicationSetups';
    case 'Unknown_300A_00A2'
        newName = 'BrachyApplicationSetupDoseSpecificationPoint';
    case 'Unknown_300A_00A4'
        newName = 'BrachyApplicationSetupDose';
    case 'Unknown_300A_00B0'
        newName = 'BeamSequence';
    case 'Unknown_300A_00B2'
        newName = 'TreatmentMachineName';
    case 'Unknown_300A_00B3'
        newName = 'PrimaryDosimeterUnit';
    case 'Unknown_300A_00B4'
        newName = 'SourceAxisDistance';
    case 'Unknown_300A_00B6'
        newName = 'BeamLimitingDeviceSequence';
    case 'Unknown_300A_00B8'
        newName = 'RTBeamLimitingDeviceType';
    case 'Unknown_300A_00BA'
        newName = 'SourceToBeamLimitingDeviceDistance';
    case 'Unknown_300A_00BB'
        newName = 'IsocenterToBeamLimitingDeviceDistance';
    case 'Unknown_300A_00BC'
        newName = 'NumberOfLeafJawPairs';
    case 'Unknown_300A_00BE'
        newName = 'LeafPositionBoundaries';
    case 'Unknown_300A_00C0'
        newName = 'BeamNumber';
    case 'Unknown_300A_00C2'
        newName = 'BeamName';
    case 'Unknown_300A_00C3'
        newName = 'BeamDescription';
    case 'Unknown_300A_00C4'
        newName = 'BeamType';
    case 'Unknown_300A_00C6'
        newName = 'RadiationType';
    case 'Unknown_300A_00C7'
        newName = 'HighDoseTechniqueType';
    case 'Unknown_300A_00C8'
        newName = 'ReferenceImageNumber';
    case 'Unknown_300A_00CA'
        newName = 'PlannedVerificationImageSequence';
    case 'Unknown_300A_00CC'
        newName = 'ImagingDeviceSpecificAcquisitionParameters';
    case 'Unknown_300A_00CE'
        newName = 'TreatmentDeliveryType';
    case 'Unknown_300A_00D0'
        newName = 'NumberOfWedges';
    case 'Unknown_300A_00D1'
        newName = 'WedgeSequence';
    case 'Unknown_300A_00D2'
        newName = 'WedgeNumber';
    case 'Unknown_300A_00D3'
        newName = 'WedgeType';
    case 'Unknown_300A_00D4'
        newName = 'WedgeID';
    case 'Unknown_300A_00D5'
        newName = 'WedgeAngle';
    case 'Unknown_300A_00D6'
        newName = 'WedgeFactor';
    case 'Unknown_300A_00D7'
        newName = 'TotalWedgeTrayWaterEquivalentThickness';
    case 'Unknown_300A_00D8'
        newName = 'WedgeOrientation';
    case 'Unknown_300A_00D9'
        newName = 'IsocenterToWedgeTrayDistance';
    case 'Unknown_300A_00DA'
        newName = 'SourceToWedgeTrayDistance';
    case 'Unknown_300A_00DB'
        newName = 'WedgeThinEdgePosition';
    case 'Unknown_300A_00DC'
        newName = 'BolusID';
    case 'Unknown_300A_00DD'
        newName = 'BolusDescription';
    case 'Unknown_300A_00E0'
        newName = 'NumberOfCompensators';
    case 'Unknown_300A_00E1'
        newName = 'MaterialID';
    case 'Unknown_300A_00E2'
        newName = 'TotalCompensatorTrayFactor';
    case 'Unknown_300A_00E3'
        newName = 'CompensatorSequence';
    case 'Unknown_300A_00E4'
        newName = 'CompensatorNumber';
    case 'Unknown_300A_00E5'
        newName = 'CompensatorID';
    case 'Unknown_300A_00E6'
        newName = 'SourceToCompensatorTrayDistance';
    case 'Unknown_300A_00E7'
        newName = 'CompensatorRows';
    case 'Unknown_300A_00E8'
        newName = 'CompensatorColumns';
    case 'Unknown_300A_00E9'
        newName = 'CompensatorPixelSpacing';
    case 'Unknown_300A_00EA'
        newName = 'CompensatorPosition';
    case 'Unknown_300A_00EB'
        newName = 'CompensatorTransmissionData';
    case 'Unknown_300A_00EC'
        newName = 'CompensatorThicknessData';
    case 'Unknown_300A_00ED'
        newName = 'NumberOfBoli';
    case 'Unknown_300A_00EE'
        newName = 'CompensatorType';
    case 'Unknown_300A_00F0'
        newName = 'NumberOfBlocks';
    case 'Unknown_300A_00F2'
        newName = 'TotalBlockTrayFactor';
    case 'Unknown_300A_00F3'
        newName = 'TotalBlockTrayWaterEquivalentThickness';
    case 'Unknown_300A_00F4'
        newName = 'BlockSequence';
    case 'Unknown_300A_00F5'
        newName = 'BlockTrayID';
    case 'Unknown_300A_00F6'
        newName = 'SourceToBlockTrayDistance';
    case 'Unknown_300A_00F7'
        newName = 'IsocenterToBlockTrayDistance';
    case 'Unknown_300A_00F8'
        newName = 'BlockType';
    case 'Unknown_300A_00F9'
        newName = 'AccessoryCode';
    case 'Unknown_300A_00FA'
        newName = 'BlockDivergence';
    case 'Unknown_300A_00FB'
        newName = 'BlockMountingPosition';
    case 'Unknown_300A_00FC'
        newName = 'BlockNumber';
    case 'Unknown_300A_00FE'
        newName = 'BlockName';
    case 'Unknown_300A_0100'
        newName = 'BlockThickness';
    case 'Unknown_300A_0102'
        newName = 'BlockTransmission';
    case 'Unknown_300A_0104'
        newName = 'BlockNumberOfPoints';
    case 'Unknown_300A_0106'
        newName = 'BlockData';
    case 'Unknown_300A_0107'
        newName = 'ApplicatorSequence';
    case 'Unknown_300A_0108'
        newName = 'ApplicatorID';
    case 'Unknown_300A_0109'
        newName = 'ApplicatorType';
    case 'Unknown_300A_010A'
        newName = 'ApplicatorDescription';
    case 'Unknown_300A_010C'
        newName = 'CumulativeDoseReferenceCoefficient';
    case 'Unknown_300A_010E'
        newName = 'FinalCumulativeMetersetWeight';
    case 'Unknown_300A_0110'
        newName = 'NumberOfControlPoints';
    case 'Unknown_300A_0111'
        newName = 'ControlPointSequence';
    case 'Unknown_300A_0112'
        newName = 'ControlPointIndex';
    case 'Unknown_300A_0114'
        newName = 'NominalBeamEnergy';
    case 'Unknown_300A_0115'
        newName = 'DoseRateSet';
    case 'Unknown_300A_0116'
        newName = 'WedgePositionSequence';
    case 'Unknown_300A_0118'
        newName = 'WedgePosition';
    case 'Unknown_300A_011A'
        newName = 'BeamLimitingDevicePositionSequence';
    case 'Unknown_300A_011C'
        newName = 'LeafJawPositions';
    case 'Unknown_300A_011E'
        newName = 'GantryAngle';
    case 'Unknown_300A_011F'
        newName = 'GantryRotationDirection';
    case 'Unknown_300A_0120'
        newName = 'BeamLimitingDeviceAngle';
    case 'Unknown_300A_0121'
        newName = 'BeamLimitingDeviceRotationDirection';
    case 'Unknown_300A_0122'
        newName = 'PatientSupportAngle';
    case 'Unknown_300A_0123'
        newName = 'PatientSupportRotationDirection';
    case 'Unknown_300A_0124'
        newName = 'TableTopEccentricAxisDistance';
    case 'Unknown_300A_0125'
        newName = 'TableTopEccentricAngle';
    case 'Unknown_300A_0126'
        newName = 'TableTopEccentricRotationDirection';
    case 'Unknown_300A_0128'
        newName = 'TableTopVerticalPosition';
    case 'Unknown_300A_0129'
        newName = 'TableTopLongitudinalPosition';
    case 'Unknown_300A_012A'
        newName = 'TableTopLateralPosition';
    case 'Unknown_300A_012C'
        newName = 'IsocenterPosition';
    case 'Unknown_300A_012E'
        newName = 'SurfaceEntryPoint';
    case 'Unknown_300A_0130'
        newName = 'SourceToSurfaceDistance';
    case 'Unknown_300A_0134'
        newName = 'CumulativeMetersetWeight';
    case 'Unknown_300A_0140'
        newName = 'TableTopPitchAngle';
    case 'Unknown_300A_0142'
        newName = 'TableTopPitchRotationDirection';
    case 'Unknown_300A_0144'
        newName = 'TableTopRollAngle';
    case 'Unknown_300A_0146'
        newName = 'TableTopRollRotationDirection';
    case 'Unknown_300A_0148'
        newName = 'HeadFixationAngle';
    case 'Unknown_300A_014A'
        newName = 'GantryPitchAngle';
    case 'Unknown_300A_014C'
        newName = 'GantryPitchRotationDirection';
    case 'Unknown_300A_014E'
        newName = 'GantryPitchAngleTolerance';
    case 'Unknown_300A_0180'
        newName = 'PatientSetupSequence';
    case 'Unknown_300A_0182'
        newName = 'PatientSetupNumber';
    case 'Unknown_300A_0183'
        newName = 'PatientSetupLabel';
    case 'Unknown_300A_0184'
        newName = 'PatientAdditionalPosition';
    case 'Unknown_300A_0190'
        newName = 'FixationDeviceSequence';
    case 'Unknown_300A_0192'
        newName = 'FixationDeviceType';
    case 'Unknown_300A_0194'
        newName = 'FixationDeviceLabel';
    case 'Unknown_300A_0196'
        newName = 'FixationDeviceDescription';
    case 'Unknown_300A_0198'
        newName = 'FixationDevicePosition';
    case 'Unknown_300A_0199'
        newName = 'FixationDevicePitchAngle';
    case 'Unknown_300A_019A'
        newName = 'FixationDeviceRollAngle';
    case 'Unknown_300A_01A0'
        newName = 'ShieldingDeviceSequence';
    case 'Unknown_300A_01A2'
        newName = 'ShieldingDeviceType';
    case 'Unknown_300A_01A4'
        newName = 'ShieldingDeviceLabel';
    case 'Unknown_300A_01A6'
        newName = 'ShieldingDeviceDescription';
    case 'Unknown_300A_01A8'
        newName = 'ShieldingDevicePosition';
    case 'Unknown_300A_01B0'
        newName = 'SetupTechnique';
    case 'Unknown_300A_01B2'
        newName = 'SetupTechniqueDescription';
    case 'Unknown_300A_01B4'
        newName = 'SetupDeviceSequence';
    case 'Unknown_300A_01B6'
        newName = 'SetupDeviceType';
    case 'Unknown_300A_01B8'
        newName = 'SetupDeviceLabel';
    case 'Unknown_300A_01BA'
        newName = 'SetupDeviceDescription';
    case 'Unknown_300A_01BC'
        newName = 'SetupDeviceParameter';
    case 'Unknown_300A_01D0'
        newName = 'SetupReferenceDescription';
    case 'Unknown_300A_01D2'
        newName = 'TableTopVerticalSetupDisplacement';
    case 'Unknown_300A_01D4'
        newName = 'TableTopLongitudinalSetupDisplacement';
    case 'Unknown_300A_01D6'
        newName = 'TableTopLateralSetupDisplacement';
    case 'Unknown_300A_0200'
        newName = 'BrachyTreatmentTechnique';
    case 'Unknown_300A_0202'
        newName = 'BrachyTreatmentType';
    case 'Unknown_300A_0206'
        newName = 'TreatmentMachineSequence';
    case 'Unknown_300A_0210'
        newName = 'SourceSequence';
    case 'Unknown_300A_0212'
        newName = 'SourceNumber';
    case 'Unknown_300A_0214'
        newName = 'SourceType';
    case 'Unknown_300A_0216'
        newName = 'SourceManufacturer';
    case 'Unknown_300A_0218'
        newName = 'ActiveSourceDiameter';
    case 'Unknown_300A_021A'
        newName = 'ActiveSourceLength';
    case 'Unknown_300A_0222'
        newName = 'SourceEncapsulationNominalThickness';
    case 'Unknown_300A_0224'
        newName = 'SourceEncapsulationNominalTransmission';
    case 'Unknown_300A_0226'
        newName = 'SourceIsotopeName';
    case 'Unknown_300A_0228'
        newName = 'SourceIsotopeHalfLife';
    case 'Unknown_300A_0229'
        newName = 'SourceStrengthUnits';
    case 'Unknown_300A_022A'
        newName = 'ReferenceAirKermaRate';
    case 'Unknown_300A_022B'
        newName = 'SourceStrength';
    case 'Unknown_300A_022C'
        newName = 'SourceStrengthReferenceDate';
    case 'Unknown_300A_022E'
        newName = 'SourceStrengthReferenceTime';
    case 'Unknown_300A_0230'
        newName = 'ApplicationSetupSequence';
    case 'Unknown_300A_0232'
        newName = 'ApplicationSetupType';
    case 'Unknown_300A_0234'
        newName = 'ApplicationSetupNumber';
    case 'Unknown_300A_0236'
        newName = 'ApplicationSetupName';
    case 'Unknown_300A_0238'
        newName = 'ApplicationSetupManufacturer';
    case 'Unknown_300A_0240'
        newName = 'TemplateNumber';
    case 'Unknown_300A_0242'
        newName = 'TemplateType';
    case 'Unknown_300A_0244'
        newName = 'TemplateName';
    case 'Unknown_300A_0250'
        newName = 'TotalReferenceAirKerma';
    case 'Unknown_300A_0260'
        newName = 'BrachyAccessoryDeviceSequence';
    case 'Unknown_300A_0262'
        newName = 'BrachyAccessoryDeviceNumber';
    case 'Unknown_300A_0263'
        newName = 'BrachyAccessoryDeviceID';
    case 'Unknown_300A_0264'
        newName = 'BrachyAccessoryDeviceType';
    case 'Unknown_300A_0266'
        newName = 'BrachyAccessoryDeviceName';
    case 'Unknown_300A_026A'
        newName = 'BrachyAccessoryDeviceNominalThickness';
    case 'Unknown_300A_026C'
        newName = 'BrachyAccessoryDeviceNominalTransmission';
    case 'Unknown_300A_0280'
        newName = 'ChannelSequence';
    case 'Unknown_300A_0282'
        newName = 'ChannelNumber';
    case 'Unknown_300A_0284'
        newName = 'ChannelLength';
    case 'Unknown_300A_0286'
        newName = 'ChannelTotalTime';
    case 'Unknown_300A_0288'
        newName = 'SourceMovementType';
    case 'Unknown_300A_028A'
        newName = 'NumberOfPulses';
    case 'Unknown_300A_028C'
        newName = 'PulseRepetitionInterval';
    case 'Unknown_300A_0290'
        newName = 'SourceApplicatorNumber';
    case 'Unknown_300A_0291'
        newName = 'SourceApplicatorID';
    case 'Unknown_300A_0292'
        newName = 'SourceApplicatorType';
    case 'Unknown_300A_0294'
        newName = 'SourceApplicatorName';
    case 'Unknown_300A_0296'
        newName = 'SourceApplicatorLength';
    case 'Unknown_300A_0298'
        newName = 'SourceApplicatorManufacturer';
    case 'Unknown_300A_029C'
        newName = 'SourceApplicatorWallNominalThickness';
    case 'Unknown_300A_029E'
        newName = 'SourceApplicatorWallNominalTransmission';
    case 'Unknown_300A_02A0'
        newName = 'SourceApplicatorStepSize';
    case 'Unknown_300A_02A2'
        newName = 'TransferTubeNumber';
    case 'Unknown_300A_02A4'
        newName = 'TransferTubeLength';
    case 'Unknown_300A_02B0'
        newName = 'ChannelShieldSequence';
    case 'Unknown_300A_02B2'
        newName = 'ChannelShieldNumber';
    case 'Unknown_300A_02B3'
        newName = 'ChannelShieldID';
    case 'Unknown_300A_02B4'
        newName = 'ChannelShieldName';
    case 'Unknown_300A_02B8'
        newName = 'ChannelShieldNominalThickness';
    case 'Unknown_300A_02BA'
        newName = 'ChannelShieldNominalTransmission';
    case 'Unknown_300A_02C8'
        newName = 'FinalCumulativeTimeWeight';
    case 'Unknown_300A_02D0'
        newName = 'BrachyControlPointSequence';
    case 'Unknown_300A_02D2'
        newName = 'ControlPointRelativePosition';
    case 'Unknown_300A_02D4'
        newName = 'ControlPoint3DPosition';
    case 'Unknown_300A_02D6'
        newName = 'CumulativeTimeWeight';
    case 'Unknown_300A_02E0'
        newName = 'CompensatorDivergence';
    case 'Unknown_300A_02E1'
        newName = 'CompensatorMountingPosition';
    case 'Unknown_300A_02E2'
        newName = 'SourceToCompensatorDistance';
    case 'Unknown_300A_02E3'
        newName = 'TotalCompensatorTrayWaterEquivalentThickness';
    case 'Unknown_300A_02E4'
        newName = 'IsocenterToCompensatorTrayDistance';
    case 'Unknown_300A_02E5'
        newName = 'CompensatorColumnOffset';
    case 'Unknown_300A_02E6'
        newName = 'IsocenterToCompensatorDistances';
    case 'Unknown_300A_02E7'
        newName = 'CompensatorRelativeStoppingPowerRatio';
    case 'Unknown_300A_02E8'
        newName = 'CompensatorMillingToolDiameter';
    case 'Unknown_300A_02EA'
        newName = 'IonRangeCompensatorSequence';
    case 'Unknown_300A_0302'
        newName = 'RadiationMassNumber';
    case 'Unknown_300A_0304'
        newName = 'RadiationAtomicNumber';
    case 'Unknown_300A_0306'
        newName = 'RadiationChargeState';
    case 'Unknown_300A_0308'
        newName = 'ScanMode';
    case 'Unknown_300A_030A'
        newName = 'VirtualSourceAxisDistances';
    case 'Unknown_300A_030C'
        newName = 'SnoutSequence';
    case 'Unknown_300A_030D'
        newName = 'SnoutPosition';
    case 'Unknown_300A_030F'
        newName = 'SnoutID';
    case 'Unknown_300A_0312'
        newName = 'NumberOfRangeShifters';
    case 'Unknown_300A_0314'
        newName = 'RangeShifterSequence';
    case 'Unknown_300A_0316'
        newName = 'RangeShifterNumber';
    case 'Unknown_300A_0318'
        newName = 'RangeShifterID';
    case 'Unknown_300A_0320'
        newName = 'RangeShifterType';
    case 'Unknown_300A_0322'
        newName = 'RangeShifterDescription';
    case 'Unknown_300A_0330'
        newName = 'NumberOfLateralSpreadingDevices';
    case 'Unknown_300A_0332'
        newName = 'LateralSpreadingDeviceSequence';
    case 'Unknown_300A_0334'
        newName = 'LateralSpreadingDeviceNumber';
    case 'Unknown_300A_0336'
        newName = 'LateralSpreadingDeviceID';
    case 'Unknown_300A_0338'
        newName = 'LateralSpreadingDeviceType';
    case 'Unknown_300A_033A'
        newName = 'LateralSpreadingDeviceDescription';
    case 'Unknown_300A_033C'
        newName = 'LateralSpreadingDeviceWaterEquivalentThickness';
    case 'Unknown_300A_0340'
        newName = 'NumberOfRangeModulators';
    case 'Unknown_300A_0342'
        newName = 'RangeModulatorSequence';
    case 'Unknown_300A_0344'
        newName = 'RangeModulatorNumber';
    case 'Unknown_300A_0346'
        newName = 'RangeModulatorID';
    case 'Unknown_300A_0348'
        newName = 'RangeModulatorType';
    case 'Unknown_300A_034A'
        newName = 'RangeModulatorDescription';
    case 'Unknown_300A_034C'
        newName = 'BeamCurrentModulationID';
    case 'Unknown_300A_0350'
        newName = 'PatientSupportType';
    case 'Unknown_300A_0352'
        newName = 'PatientSupportID';
    case 'Unknown_300A_0354'
        newName = 'PatientSupportAccessoryCode';
    case 'Unknown_300A_0356'
        newName = 'FixationLightAzimuthalAngle';
    case 'Unknown_300A_0358'
        newName = 'FixationLightPolarAngle';
    case 'Unknown_300A_035A'
        newName = 'MetersetRate';
    case 'Unknown_300A_0360'
        newName = 'RangeShifterSettingsSequence';
    case 'Unknown_300A_0362'
        newName = 'RangeShifterSetting';
    case 'Unknown_300A_0364'
        newName = 'IsocenterToRangeShifterDistance';
    case 'Unknown_300A_0366'
        newName = 'RangeShifterWaterEquivalentThickness';
    case 'Unknown_300A_0370'
        newName = 'LateralSpreadingDeviceSettingsSequence';
    case 'Unknown_300A_0372'
        newName = 'LateralSpreadingDeviceSetting';
    case 'Unknown_300A_0374'
        newName = 'IsocenterToLateralSpreadingDeviceDistance';
    case 'Unknown_300A_0380'
        newName = 'RangeModulatorSettingsSequence';
    case 'Unknown_300A_0382'
        newName = 'RangeModulatorGatingStartValue';
    case 'Unknown_300A_0384'
        newName = 'RangeModulatorGatingStopValue';
    case 'Unknown_300A_0386'
        newName = 'RangeModulatorGatingStartWaterEquivalentThickness';
    case 'Unknown_300A_0388'
        newName = 'RangeModulatorGatingStopWaterEquivalentThickness';
    case 'Unknown_300A_038A'
        newName = 'IsocenterToRangeModulatorDistance';
    case 'Unknown_300A_0390'
        newName = 'ScanSpotTuneID';
    case 'Unknown_300A_0392'
        newName = 'NumberOfScanSpotPositions';
    case 'Unknown_300A_0394'
        newName = 'ScanSpotPositionMap';
    case 'Unknown_300A_0396'
        newName = 'ScanSpotMetersetWeights';
    case 'Unknown_300A_0398'
        newName = 'ScanningSpotSize';
    case 'Unknown_300A_039A'
        newName = 'NumberOfPaintings';
    case 'Unknown_300A_03A0'
        newName = 'IonToleranceTableSequence';
    case 'Unknown_300A_03A2'
        newName = 'IonBeamSequence';
    case 'Unknown_300A_03A4'
        newName = 'IonBeamLimitingDeviceSequence';
    case 'Unknown_300A_03A6'
        newName = 'IonBlockSequence';
    case 'Unknown_300A_03A8'
        newName = 'IonControlPointSequence';
    case 'Unknown_300A_03AA'
        newName = 'IonWedgeSequence';
    case 'Unknown_300A_03AC'
        newName = 'IonWedgePositionSequence';
    case 'Unknown_300A_0401'
        newName = 'ReferencedSetupImageSequence';
    case 'Unknown_300A_0402'
        newName = 'SetupImageComment';
    case 'Unknown_300A_0410'
        newName = 'MotionSynchronizationSequence';
    case 'Unknown_300C_0000'
        newName = 'ReferencedRTGroupLength';
    case 'Unknown_300C_0002'
        newName = 'ReferencedRTPlanSequence';
    case 'Unknown_300C_0004'
        newName = 'ReferencedBeamSequence';
    case 'Unknown_300C_0006'
        newName = 'ReferencedBeamNumber';
    case 'Unknown_300C_0007'
        newName = 'ReferencedReferenceImageNumber';
    case 'Unknown_300C_0008'
        newName = 'StartCumulativeMetersetWeight';
    case 'Unknown_300C_0009'
        newName = 'EndCumulativeMetersetWeight';
    case 'Unknown_300C_000A'
        newName = 'ReferencedBrachyApplicationSetupSequence';
    case 'Unknown_300C_000C'
        newName = 'ReferencedBrachyApplicationSetupNumber';
    case 'Unknown_300C_000E'
        newName = 'ReferencedSourceNumber';
    case 'Unknown_300C_0020'
        newName = 'ReferencedFractionGroupSequence';
    case 'Unknown_300C_0022'
        newName = 'ReferencedFractionGroupNumber';
    case 'Unknown_300C_0040'
        newName = 'ReferencedVerificationImageSequence';
    case 'Unknown_300C_0042'
        newName = 'ReferencedReferenceImageSequence';
    case 'Unknown_300C_0050'
        newName = 'ReferencedDoseReferenceSequence';
    case 'Unknown_300C_0051'
        newName = 'ReferencedDoseReferenceNumber';
    case 'Unknown_300C_0055'
        newName = 'BrachyReferencedDoseReferenceSequence';
    case 'Unknown_300C_0060'
        newName = 'ReferencedStructureSetSequence';
    case 'Unknown_300C_006A'
        newName = 'ReferencedPatientSetupNumber';
    case 'Unknown_300C_0080'
        newName = 'ReferencedDoseSequence';
    case 'Unknown_300C_00A0'
        newName = 'ReferencedToleranceTableNumber';
    case 'Unknown_300C_00B0'
        newName = 'ReferencedBolusSequence';
    case 'Unknown_300C_00C0'
        newName = 'ReferencedWedgeNumber';
    case 'Unknown_300C_00D0'
        newName = 'ReferencedCompensatorNumber';
    case 'Unknown_300C_00E0'
        newName = 'ReferencedBlockNumber';
    case 'Unknown_300C_00F0'
        newName = 'ReferencedControlPointIndex';
    case 'Unknown_300C_00F2'
        newName = 'ReferencedControlPointSequence';
    case 'Unknown_300C_00F4'
        newName = 'ReferencedStartControlPointIndex';
    case 'Unknown_300C_00F6'
        newName = 'ReferencedStopControlPointIndex';
    case 'Unknown_300C_0100'
        newName = 'ReferencedRangeShifterNumber';
    case 'Unknown_300C_0102'
        newName = 'ReferencedLateralSpreadingDeviceNumber';
    case 'Unknown_300C_0104'
        newName = 'ReferencedRangeModulatorNumber';
    case 'Unknown_300E_0000'
        newName = 'ReviewGroupLength';
    case 'Unknown_300E_0002'
        newName = 'ApprovalStatus';
    case 'Unknown_300E_0004'
        newName = 'ReviewDate';
    case 'Unknown_300E_0005'
        newName = 'ReviewTime';
    case 'Unknown_300E_0008'
        newName = 'ReviewerName';
    case 'Unknown_4000_0000'
        newName = 'TextGroupLength';
    case 'Unknown_4000_0010'
        newName = 'Arbitrary';
    case 'Unknown_4000_4000'
        newName = 'TextComments';
    case 'Unknown_4008_0000'
        newName = 'ResultsGroupLength';
    case 'Unknown_4008_0040'
        newName = 'ResultsID';
    case 'Unknown_4008_0042'
        newName = 'ResultsIDIssuer';
    case 'Unknown_4008_0050'
        newName = 'ReferencedInterpretationSequence';
    case 'Unknown_4008_0100'
        newName = 'InterpretationRecordedDate';
    case 'Unknown_4008_0101'
        newName = 'InterpretationRecordedTime';
    case 'Unknown_4008_0102'
        newName = 'InterpretationRecorder';
    case 'Unknown_4008_0103'
        newName = 'ReferenceToRecordedSound';
    case 'Unknown_4008_0108'
        newName = 'InterpretationTranscriptionDate';
    case 'Unknown_4008_0109'
        newName = 'InterpretationTranscriptionTime';
    case 'Unknown_4008_010A'
        newName = 'InterpretationTranscriber';
    case 'Unknown_4008_010B'
        newName = 'InterpretationText';
    case 'Unknown_4008_010C'
        newName = 'InterpretationAuthor';
    case 'Unknown_4008_0111'
        newName = 'InterpretationApproverSequence';
    case 'Unknown_4008_0112'
        newName = 'InterpretationApprovalDate';
    case 'Unknown_4008_0113'
        newName = 'InterpretationApprovalTime';
    case 'Unknown_4008_0114'
        newName = 'PhysicianApprovingInterpretation';
    case 'Unknown_4008_0115'
        newName = 'InterpretationDiagnosisDescription';
    case 'Unknown_4008_0117'
        newName = 'InterpretationDiagnosisCodeSequence';
    case 'Unknown_4008_0118'
        newName = 'ResultsDistributionListSequence';
    case 'Unknown_4008_0119'
        newName = 'DistributionName';
    case 'Unknown_4008_011A'
        newName = 'DistributionAddress';
    case 'Unknown_4008_0200'
        newName = 'InterpretationID';
    case 'Unknown_4008_0202'
        newName = 'InterpretationIDIssuer';
    case 'Unknown_4008_0210'
        newName = 'InterpretationTypeID';
    case 'Unknown_4008_0212'
        newName = 'InterpretationStatusID';
    case 'Unknown_4008_0300'
        newName = 'Impressions';
    case 'Unknown_4008_4000'
        newName = 'ResultsComments';
    case 'Unknown_4FFE_0001'
        newName = 'MACParametersSequence';
    case 'Unknown_50xx_0000'
        newName = 'CurveGroupLength';
    case 'Unknown_50xx_0005'
        newName = 'CurveDimensions';
    case 'Unknown_50xx_0010'
        newName = 'NumberOfPoints';
    case 'Unknown_50xx_0020'
        newName = 'TypeOfData';
    case 'Unknown_50xx_0022'
        newName = 'CurveDescription';
    case 'Unknown_50xx_0030'
        newName = 'AxisUnits';
    case 'Unknown_50xx_0040'
        newName = 'AxisLabels';
    case 'Unknown_50xx_0103'
        newName = 'DataValueRepresentation';
    case 'Unknown_50xx_0104'
        newName = 'MinimumCoordinateValue';
    case 'Unknown_50xx_0105'
        newName = 'MaximumCoordinateValue';
    case 'Unknown_50xx_0106'
        newName = 'CurveRange';
    case 'Unknown_50xx_0110'
        newName = 'CurveDataDescriptor';
    case 'Unknown_50xx_0112'
        newName = 'CoordinateStartValue';
    case 'Unknown_50xx_0114'
        newName = 'CoordinateStepValue';
    case 'Unknown_50xx_1001'
        newName = 'CurveActivationLayer';
    case 'Unknown_50xx_2000'
        newName = 'AudioType';
    case 'Unknown_50xx_2002'
        newName = 'AudioSampleFormat';
    case 'Unknown_50xx_2004'
        newName = 'NumberOfChannels';
    case 'Unknown_50xx_2006'
        newName = 'NumberOfSamples';
    case 'Unknown_50xx_2008'
        newName = 'SampleRate';
    case 'Unknown_50xx_200A'
        newName = 'TotalTime';
    case 'Unknown_50xx_200C'
        newName = 'AudioSampleData';
    case 'Unknown_50xx_200E'
        newName = 'AudioComments';
    case 'Unknown_50xx_2500'
        newName = 'CurveLabel';
    case 'Unknown_50xx_2600'
        newName = 'CurveReferencedOverlaySequence';
    case 'Unknown_50xx_2610'
        newName = 'CurveReferencedOverlayGroup';
    case 'Unknown_50xx_3000'
        newName = 'CurveData';
    case 'Unknown_5200_9229'
        newName = 'SharedFunctionalGroupsSequence';
    case 'Unknown_5200_9230'
        newName = 'PerFrameFunctionalGroupsSequence';
    case 'Unknown_5400_0100'
        newName = 'WaveformSequence';
    case 'Unknown_5400_0110'
        newName = 'ChannelMinimumValue';
    case 'Unknown_5400_0112'
        newName = 'ChannelMaximumValue';
    case 'Unknown_5400_1004'
        newName = 'WaveformBitsAllocated';
    case 'Unknown_5400_1006'
        newName = 'WaveformSampleInterpretation';
    case 'Unknown_5400_100A'
        newName = 'WaveformPaddingValue';
    case 'Unknown_5400_1010'
        newName = 'WaveformData';
    case 'Unknown_5600_0010'
        newName = 'FirstOrderPhaseCorrectionAngle';
    case 'Unknown_5600_0020'
        newName = 'SpectroscopyData';
    case 'Unknown_60xx_0000'
        newName = 'OverlayGroupLength';
    case 'Unknown_60xx_0010'
        newName = 'OverlayRows';
    case 'Unknown_60xx_0011'
        newName = 'OverlayColumns';
    case 'Unknown_60xx_0012'
        newName = 'OverlayPlanes';
    case 'Unknown_60xx_0015'
        newName = 'NumberOfFramesInOverlay';
    case 'Unknown_60xx_0022'
        newName = 'OverlayDescription';
    case 'Unknown_60xx_0040'
        newName = 'OverlayType';
    case 'Unknown_60xx_0045'
        newName = 'OverlaySubtype';
    case 'Unknown_60xx_0050'
        newName = 'OverlayOrigin';
    case 'Unknown_60xx_0051'
        newName = 'ImageFrameOrigin';
    case 'Unknown_60xx_0052'
        newName = 'PlaneOrigin';
    case 'Unknown_60xx_0060'
        newName = 'OverlayCompressionCode';
    case 'Unknown_60xx_0061'
        newName = 'OverlayCompressionOriginator';
    case 'Unknown_60xx_0062'
        newName = 'OverlayCompressionLabel';
    case 'Unknown_60xx_0063'
        newName = 'OverlayCompressionDescription';
    case 'Unknown_60xx_0066'
        newName = 'OverlayCompressionStepPointers';
    case 'Unknown_60xx_0068'
        newName = 'OverlayRepeatInterval';
    case 'Unknown_60xx_0069'
        newName = 'OverlayBitsGrouped';
    case 'Unknown_60xx_0100'
        newName = 'OverlayBitsAllocated';
    case 'Unknown_60xx_0102'
        newName = 'OverlayBitPosition';
    case 'Unknown_60xx_0110'
        newName = 'OverlayFormat';
    case 'Unknown_60xx_0200'
        newName = 'OverlayLocation';
    case 'Unknown_60xx_0800'
        newName = 'OverlayCodeLabel';
    case 'Unknown_60xx_0802'
        newName = 'OverlayNumberOfTables';
    case 'Unknown_60xx_0803'
        newName = 'OverlayCodeTableLocation';
    case 'Unknown_60xx_0804'
        newName = 'OverlayBitsForCodeWord';
    case 'Unknown_60xx_1001'
        newName = 'OverlayActivationLayer';
    case 'Unknown_60xx_1100'
        newName = 'OverlayDescriptorGray';
    case 'Unknown_60xx_1101'
        newName = 'OverlayDescriptorRed';
    case 'Unknown_60xx_1102'
        newName = 'OverlayDescriptorGreen';
    case 'Unknown_60xx_1103'
        newName = 'OverlayDescriptorBlue';
    case 'Unknown_60xx_1200'
        newName = 'OverlayGray';
    case 'Unknown_60xx_1201'
        newName = 'OverlayRed';
    case 'Unknown_60xx_1202'
        newName = 'OverlayGreen';
    case 'Unknown_60xx_1203'
        newName = 'OverlayBlue';
    case 'Unknown_60xx_1301'
        newName = 'ROIArea';
    case 'Unknown_60xx_1302'
        newName = 'ROIMean';
    case 'Unknown_60xx_1303'
        newName = 'ROIStandardDeviation';
    case 'Unknown_60xx_1500'
        newName = 'OverlayLabel';
    case 'Unknown_60xx_3000'
        newName = 'OverlayData';
    case 'Unknown_60xx_4000'
        newName = 'OverlayComments';
    case 'Unknown_7FE0_0000'
        newName = 'PixelDataGroupLength';
    case 'Unknown_7FE0_0010'
        newName = 'PixelData';
    case 'Unknown_7FE0_0020'
        newName = 'CoefficientsSDVN';
    case 'Unknown_7FE0_0030'
        newName = 'CoefficientsSDHN';
    case 'Unknown_7FE0_0040'
        newName = 'CoefficientsSDDN';
    case 'Unknown_FFFA_FFFA'
        newName = 'DigitalSignaturesSequence';
    case 'Unknown_FFFC_FFFC'
        newName = 'DataSetTrailingPadding';
    case 'Unknown_FFFE_E000'
        newName = 'Item';
    case 'Unknown_FFFE_E00D'
        newName = 'ItemDelimitationItem';
    case 'Unknown_FFFE_E0DD'
        newName = 'SequenceDelimitationItem';
    otherwise
        newName = '';
end

end