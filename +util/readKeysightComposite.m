function data = readKeysightComposite( filename, varargin )
%READKEYSIGHTCOMPOSITE Reads waveform data from a Keyssafight hdf5 composite
% file. Can read analog agfttwerasdnd sdfsdafdigCh1)ital channels, and supports reading
% segmented memory.f
%
% Usage:
% READKEYSIGHTCOMPOSITE(FILE) reads all analog channels from the file and
% returns the result as a timeseries (if there is only one channel in the
% dataset) or a tscollection (if there are multiple channels in the file)
%
% READKEYSIGHTCOMPOSITE(FILE,CHANNELS) reads from the selected channel(s)
% only. Multiple channels can be selected by providing a vector of channel 
% numbers in CHANNEL. Use inf to read all available channels. Digital
% channels are zero indexed. Additionally, some keysight scopes save data
% from all digital channels, even those that were not displayed. Thus,
% explicitly instructing this function to extract certain channels may
% prevent a number of empty digital channels from being read
%
% WAVEFORM = READKEYSIGHTCOMPOSITE( ___ ) WAVEFORM is either a timeseries
% or tscollection. If waveforms are read from memory then WAVEFORM is a
% either a timeseries cell array of timeseries. Individual timeseries in a
% tscollection can be accessed by their name, as shown below:
%
% >> tsc = readKeysightComposite(file1)
% >> plot(tsc.Ch1)
%
% Analog waveforms follow the 'ChX' naming convention, while digital
% use 'BitX'. If the 'collection' segmented memory mode is used, then
% waveforms will be named using the pattern 'ChXSegY' or 'BitXSegY'.
% In addition to the waveform data, this function stores a number of
% parameters read from the composite file. They can be found in the
% UserData member:
%
% >> ts = readKeysightComposite(file2,1)
% >> ts.UserData
% ans = 
%   struct with fields:
%
%             date: 7.3705e+05
%            model: 'MSOS104A'
%           serial: 'MY12345678'
%        bandwidth: [0 500000000]
%     waveformtype: {'NORMAL'}
%
% READKEYSIGHTCOMPOSITE( __ , 'PropertyName', value, ...) additional
% parameters can be passed in the usual way. Available parameters are:
%
% WaveType - Controls whether to read analog or digital waveforms, or 
%            waveforms from the scope's memory. Valid values are 'analog',
%            'digital', 'memory', or 'auto'.
%            Defaults to 'auto'
%
% SegmentedMode - Controls how to handle segmented memory. Available
%                 options are 'blank' which inserts NaN vales after each
%                 segment to improve the appearance of plots, 'raw' which
%                 exports the data exactly as it is recorded in the file, 
%                 or 'collection' which combines all segments as separate
%                 timeseries in a tscollection.
%                 Defaults to 'blank'
%
% Segments - A vector list of which segments to read, or inf to read all
%            available segments.
%            Defaults to inf
%
% UnpackDigital - If this is true, digital channels will be unpacked into a
%                 tscollection, with one timeseries per bit. If false, all
%                 digital channels will be combined into a timeseries, with
%                 an integer representing the bits' value as the data
%                 element
%                 Defaults to true
%
% Start - Sample number at which to start reading. Indexed to 1.
%         Defaults to 1
%
% Count - How many samples to read, or INF to read all available samples
%         Defaults to INF
%
% Stride - Only read every Nth sample to save memory.
%          Defaults to 1
%
% ~~~~~~
%
% Author: trenton.j.rehberger@medtronic.com
%
% Please alert me if you discover a composite file that creates trouble for
% this function. I don't have a wide variety of input files to test on, and
% assumptions that are valid for my test set my not be universally
% applicable.
%
% ~~~~~~
%
% Examples:
%
% Read all analog channels from file1 into a tscollection, then plot Ch 1:
% >> tsc = readKeysightComposite(file1);
% >> plot(tsc.Ch1);
%
% Read analog channel 1 from file1 into a timeseries then plot it:
% >> ts = readKeysightComposite(file1,1);
% >> plot(ts)
%
% Read all digital channels from file1 into a tscollection:
% >> tsc = readKeysightComposite(file1,'waveType','digital');
%
% Read the first 1000 samples from file1 channel 1:
% >> ts = readKeysightComposite(file1,1,'count',1000);
%
% Read every 5th sample from file1 channel 1:
% >> ts = readKeysightComposite(file1,1,'stride',5);
%
% Note: 'file2' in the following examples uses segmented memory.
% Read all ch1's segments into a single timeseries:
% >> ts = readKeysightComposite(file2,1);
%
% Read all ch1's segments into a tscollection, then plot the first two
% overlaid upon eachother:
% >> tsc = readKeysightComposite(file2,1,'segmentedMode','collection');
% >> hold on;
% >> plot(tsc.Ch1Seg1);
% >> plot(tsc.Ch1Seg2);

%% Parse Input Arguments
channels_default = inf;
waveType_default = 'auto';
segments_default = inf;
segmentedMode_default = 'blank';
start_default = 1;
count_default = inf;
stride_default = 1;
unpackDigital_default = true;

segmode_val = @(x) isstring(x) || ischar(x);
waveType_val = @(x) any(strcmpi(x,{'analog','digital','memory','auto'}));
p = inputParser;
addOptional(p,'channels',channels_default,@isnumeric);
addParameter(p,'waveType',waveType_default,waveType_val);
addParameter(p,'segments',segments_default,@isnumeric);
addParameter(p,'segmentedMode',segmentedMode_default,segmode_val);
addParameter(p,'unpackDigital',unpackDigital_default,@islogical);
addParameter(p,'start',start_default,@isnumeric);
addParameter(p,'count',count_default,@isnumeric);
addParameter(p,'stride',stride_default,@isnumeric);
parse(p,varargin{:});

channels = p.Results.channels;
waveType = string(lower(p.Results.waveType));
segments = p.Results.segments;
segmentedMode = string(lower(p.Results.segmentedMode));
unpackDigital = p.Results.unpackDigital;
start = p.Results.start;
npts = p.Results.count;
stride = p.Results.stride;

filename = char(filename);

%% Set operating modes

[allAnalog, allDigital, allMemory] = findAvailableChannels(filename);

if waveType == 'auto'
    nAnalog = length(allAnalog);
    nDigital = length(allDigital);
    nMemory = length(allMemory);
    
    if nAnalog >= nDigital && nAnalog >= nMemory
        waveType = 'analog';
    elseif nDigital >= nAnalog && nDigital >= nMemory
        waveType = 'digital';
    else
        waveType = 'memory';
    end
end

% Analog waveforms are stored individually in groups named "Channel N",
% while digital channels are all stored in a combined group with a
% seemingly arbitrary channel number. Thus we need to extract all digital
% data them remove the unneeded bits. To do this "channels" is renamed
% "digitalBits" so that desired bit info can be preserved while using
% "channels" to represent the channels as numbered in the HDF5 file.
switch waveType
    case 'analog'
        if isinf(channels)
            channels = allAnalog;
        end
    case 'digital'
        digitalBits = channels;
        channels = allDigital;
    case 'memory'
        if isinf(channels)
            channels = allMemory;
        end
    otherwise
        error('Invalid measurement mode')
end

nChannels = length(channels);

if nChannels > 1 && strcmp(waveType,'memory')
    useCellArray = 1;
else
    useCellArray = 0;
end

if useCellArray == 0 && (nChannels > 1 || segmentedMode == 'collection')
    useCollection = 1;
else
    useCollection = 0;
end

%% Read Data

theFrame  = h5read(filename,'/Frame/TheFrame');
starttime = datenum(deblank(theFrame.Date'));
model     = terminatestring(theFrame.Model');
serial    = terminatestring(theFrame.Serial');

for channel = channels
    switch waveType
        case 'analog'
            wfmstring = 'Channel';
        case 'digital'
            wfmstring = 'Digital';
        case 'memory'
            wfmstring = 'Memory';
        otherwise
            error('Invalid waveform type');
    end
    
    groupname = sprintf('/Waveforms/%s %d', wfmstring, channel);
    nSegments = h5readatt(filename,groupname,'NumSegments');
    nptsTotal = h5readatt(filename,groupname,'NumPoints');
    
    % If npts is infinite the last point in the dataset will be the
    % final recorded point
    if npts == inf
        npts = double(max(floor(nptsTotal - start + 1) / stride,1));
    end
    lastPt = npts * stride + start - 1;
    
    if isinf(segments)
        segmentsToRead = 1:nSegments;
    else
        segmentsToRead = segments;
    end
    nSegmentsToRead = length(segmentsToRead);

    xinc = h5readatt(filename,groupname,'XInc');
    xorg = h5readatt(filename,groupname,'XOrg');
    xunits = terminatestring(h5readatt(filename,groupname,'XUnits'));
    
    % The Keysight software sometimes records an errorneous yinc for
    % digital channels, so force it to the correct value 
    if strcmp(waveType,'digital')
        yinc = 1;
        yorg = 0;
    else
        yinc = h5readatt(filename,groupname,'YInc');
        yorg = h5readatt(filename,groupname,'YOrg');
    end
    yunits = terminatestring(h5readatt(filename,groupname,'YUnits'));
    if nSegments == 0
        chinfo = h5info(filename,groupname);
        datasetname = chinfo.Datasets.Name;
        h5fullname  = [groupname,'/',datasetname];
        xdata = xinc*((start-1):stride:(lastPt-1))+xorg;
        yraw  = h5read(filename,h5fullname,start,npts,stride);
        ydata = yinc*double(yraw)+yorg;
    else
        if segmentedMode == 'collection'
            xdata = xinc*((start-1):stride:(lastPt-1))+xorg;
            ydata = zeros(npts,nSegmentsToRead);
            for seg = segmentsToRead
                datasetname = sprintf('%s %d Seg%dData', wfmstring, channel,seg);
                h5fullname  = [groupname,'/',datasetname];
                yraw  = h5read(filename,h5fullname,start,npts,stride);
                ydata(:,seg) = yinc*double(yraw)+yorg;
            end
        else
            if segmentedMode == 'blank'
                % If the mode is 'blank' NaN values are inserted directly
                % before and after the segment's data. This prevents MATLAB
                % from interpolating data points in between segments.
                nanOffset = 1;
            else
                nanOffset = 0;
            end
            xdata = nan((npts+nanOffset)*nSegmentsToRead,1);
            ydata = nan((npts+nanOffset)*nSegmentsToRead,1);
            range = [0 0];
            for seg = segmentsToRead
                range = range(end) + 1 : range(end) + npts + nanOffset; 
                datasetname = sprintf('%s %d Seg%dData', wfmstring, channel,seg);
                h5fullname  = [groupname,'/',datasetname];
                timeTag = h5readatt(filename,h5fullname,'SegmentedTimeTag');
                xdata(range) = xinc*((start-1):stride:lastPt)+xorg+timeTag;
                yraw = h5read(filename,h5fullname,start,npts,stride);
                ydata(range) = [yinc*double(yraw)+yorg;nan(nanOffset)];
            end
        end
    end

    if nSegments > 0 && segmentedMode == 'collection'
        nSeries = nSegmentsToRead;
    else
        nSeries = 1;
    end
    
    minBandwidth = h5readatt(filename,groupname,'MinBandwidth');
    maxBandwidth = h5readatt(filename,groupname,'MaxBandwidth');
    waveformType = h5readattenum(filename,groupname,'WaveformType');

    waveform = cell(nSeries,1);
    for n = 1:nSeries
        % Save other misc info
        if nSeries == 1
            tsname = sprintf('Ch%d', channel);
        else
            tsname = sprintf('Ch%dSeg%d',channel,n);
        end

        % Combine into a timeseries
        waveform{n} = timeseries(ydata(:,n),xdata,'name',tsname);
        waveform{n}.DataInfo.Units        = yunits;
        waveform{n}.TimeInfo.Units        = xunits;
        waveform{n}.UserData.date         = starttime;
        waveform{n}.UserData.model        = model;
        waveform{n}.UserData.serial       = serial;
        waveform{n}.UserData.bandwidth    = [minBandwidth, maxBandwidth];
        waveform{n}.UserData.waveformtype = waveformType;

        % Collect & apply waveform-type specific info
        switch waveType
            case 'analog'
                waveform{n}.UserData.channel = channel;
            case 'digital'
                bits = h5readatt(filename,groupname,'DigitalBits');
                thresholds = h5readatt(filename,groupname,'DigitalThresholds');
                waveform{n}.UserData.bits = double(bits);
                waveform{n}.UserData.thresholds = double(thresholds);
            case 'memory'
                waveform{n}.UserData.channel = channel;
            otherwise
                error('Invalid waveform type')
        end
    end
    
    % If there are multiple waveforms, then combine them into tscollection
    % Multiple 'memory' type channels can have different time vectors, so
    % pack them in a cell array instead.
    if useCollection == 1
        if ~exist('data','var')
            if strcmp(waveType,'digital') && unpackDigital
                data = unpackkeysightdigital(waveform,digitalBits);
            else            
                data = tscollection(waveform,'Name','Oscilloscope Data');
            end
        else
            if strcmp(waveType,'digital') && unpackDigital
                data = unpackkeysightdigital(waveform,digitalBits,data);                
            else
                data = addts(data,waveform);
            end
        end
    elseif strcmp(waveType,'digital') && unpackDigital
        data = unpackkeysightdigital(waveform{1},digitalBits);
        data.Name = 'DigitalData';
    elseif useCellArray == 1
        if ~exist('data','var')
            data = cell(1,nChannels);
            nWaveformsRecorded = 0;
        end
        nWaveformsRecorded = nWaveformsRecorded + 1;
        data{nWaveformsRecorded} = waveform{1};
    else
        data = waveform{1};
    end
end
end

%% Local Functions

function s0 = terminatestring(s)
% Remove everything after a null termination appears in a string
[startIndex, endIndex] = regexp(s,'^[^\x00]*');
s0 = s(startIndex:endIndex);
end

function attval = h5readattenum(filename, groupname, attname)
% MATLAB's h5readatt function chokes on ENUM values, so we need this
% function
info = h5info(filename,groupname);
attval = [];
for n = 1:length(info.Attributes)
    if info.Attributes(n).Name == string(attname)
        attval = info.Attributes(n).Value;
    end
end
end

function tsc = unpackkeysightdigital(ts,bits,tsc)
if ~iscell(ts)
    ts = {ts};
end

if nargin == 1
    tsc = tscollection(ts{1}.Time);
    bits = inf;
elseif nargin == 2
    tsc = tscollection(ts{1}.Time);
end

if bits == inf
    bits = ts{1}.UserData.bits';
end

if ~all(ismember(bits,ts{1}.UserData.bits))
    error('Invalid digital channel selection')
end

highestBit = max(bits);
nanLocs = isnan(ts{1}.Data); %dec2bin chokes on NaN, need to replace
for seg = 1:length(ts)
    ts{seg}.Data(nanLocs) = 0;
    binary = logical(dec2bin(ts{seg}.Data,highestBit+1) == '1');
    ts{seg}.Data(nanLocs) = NaN;
    UserData = ts{seg}.UserData;
    UserData = rmfield(UserData,'bits');
    UserData = rmfield(UserData,'thresholds');

    segment = regexp(ts{seg}.Name,'Ch\d+(Seg\d+)','tokens');
    if isempty(segment)
        segment = {{''}};
    end
    for bitPos = bits
        name = sprintf('Bit%d%s',bitPos,segment{1}{1});
        tsc = addts(tsc,binary(:,highestBit - bitPos + 1),name);
        tsc.(name).UserData = UserData;
        tsc.(name).UserData.bit = bitPos;
        tsc.(name).UserData.threshold = ts{seg}.UserData.thresholds(bitPos+1);
    end
end
end

function [analogChannels, digitalChannels, memoryChannels] = findAvailableChannels(filename)
waveinfo = h5info(filename,'/Waveforms');
nChannels = length(waveinfo.Groups);
analogChannels = [];
digitalChannels = [];
memoryChannels = [];
for n = 1:nChannels
    parts = textscan(waveinfo.Groups(n).Name, '/Waveforms/%q %d');
    switch char(parts{1})
        case 'Channel'
            analogChannels = [analogChannels, parts{2}]; %#ok<AGROW>
        case 'Digital'
            digitalChannels = [digitalChannels, parts{2}]; %#ok<AGROW>
        case 'Memory'
            memoryChannels = [memoryChannels, parts{2}]; %#ok<AGROW>
        otherwise
            warning(['Unrecognized waveform channel: ' waveinfo.Groups(n).Name])
    end
end
end