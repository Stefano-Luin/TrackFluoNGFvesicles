function [result,AcqTime] = ReadInfo(id,SerIn,SerFin)
% Copyright (C) 2012 - 2022 Carmine di Rienzo (contact Stefano Luin s.luin@sns.it)

tic
autoloadBioFormats = 1;

% Toggle the stitchFiles flag to control grouping of similarly
% named files into a single dataset based on file numbering.
stitchFiles = 0;

% -- Main function - no need to edit anything past this point --

% load the Bio-Formats library into the MATLAB environment
if autoloadBioFormats
    path = fullfile(ffpath('loci_tools.jar'), 'loci_tools.jar');
    javaaddpath(path);
end

% set LuraWave license code, if available
if exist('lurawaveLicense')
    path = fullfile(fileparts(mfilename('fullpath')), 'lwf_jsdk2.6.jar');
    javaaddpath(path);
    java.lang.System.setProperty('lurawave.license', lurawaveLicense);
end

% check MATLAB version, since typecast function requires MATLAB 7.1+
canTypecast = versionCheck(version, 7, 1);

% check Bio-Formats version, since makeDataArray2D function requires trunk
bioFormatsVersion = char(loci.formats.FormatTools.VERSION);
isBioFormatsTrunk = versionCheck(bioFormatsVersion, 5, 0);

r = loci.formats.ChannelFiller();
r = loci.formats.ChannelSeparator(r);
if stitchFiles
    r = loci.formats.FileStitcher(r);
end


% i1=0;
r.setId(id);

r = bfGetReader(id, stitchFiles);
% Ome1=r.getMetadataStore();

meta1 = r.getSeriesMetadata();
keys = loci.formats.MetadataTools.keys(meta1);
numSeries1 = r.getSeriesCount();
numSeries=SerFin-SerIn+1;
if numSeries>numSeries1
    numSeries=numSeries1;
    SerIn=1;
    SerFin=numSeries;
end
result = cell(8+size(keys,1), 2,numSeries);
% result = cell(8, 2,numSeries);

AcqTime = cell(numSeries,1);
for s1 = 1:numSeries
    s=s1+SerIn-1;
    r.setSeries(s - 1);
    meta=r.getMetadata();
    meta1 = r.getSeriesMetadata();
    keys = loci.formats.MetadataTools.keys(meta1);
    Ome=r.getMetadataStore();
    AcqTime{s1}=Ome.getImageAcquisitionDate(0).getValue();
    result{1,1,s1}='width'; result{1,2,s1} = r.getSizeX();
    result{2,1,s1}='height'; result{2,2,s1} = r.getSizeY();
    result{3,1,s1}='pixelType'; result{3,2,s1} = r.getPixelType();
    result{4,1,s1}='numImages'; result{4,2,s1} = r.getImageCount();
    result{5,1,s1}=' metadataList'; result{5,2,s1} = meta;
    result{6,1,s1}=' NumSeries'; result{6,2,s1} = numSeries1;
%     result{7,1,s1}=' CreationDate'; result{7,2,s1} = r.getImageAcquisitionDate(s - 1);
    result{7,1,s1}='depth'; result{7,2,s1} = r.getSizeZ();
    result{8,1,s1}='TimeStep'; result{8,2,s1} = r.getSizeT();
    for Ikeys=1:size(keys,1)
        result{8+Ikeys,1,s1}=keys(Ikeys); result{8+Ikeys,2,s1} = get(meta1,keys(Ikeys));

    end
end
r.close();
toc




% -- Helper functions --

function [result] = versionCheck(v, maj, min)

tokens = regexp(v, '[^\d]*(\d+)[^\d]+(\d+).*', 'tokens');
majToken = tokens{1}(1);
minToken = tokens{1}(2);
major = str2num(majToken{1});
minor = str2num(minToken{1});
result = major > maj || (major == maj && minor >= min);
