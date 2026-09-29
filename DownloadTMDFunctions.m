function pathTMD = DownloadTMDFunctions(pathMain)

if nargin < 1 || isempty(pathMain)
    pathMain = fileparts(mfilename('fullpath'));
end

pathTMD = fullfile(pathMain, 'TMD_functions');
pathTMDPredict = fullfile(pathTMD, 'tmd_predict.m');

if exist(pathTMDPredict, 'file')
    fprintf('TMD functions found:\n%s\n', pathTMD);
    return;
end

if exist(pathTMD, 'dir')
    error('TMD_functions exists, but tmd_predict.m is missing: %s', pathTMD);
end

tmdCommit = 'fd89bd61e921dd81d95eda227ac23a995a6ca650';
downloadURL = ['https://github.com/chadagreene/' ...
    'Tide-Model-Driver/archive/' tmdCommit '.zip'];

downloadFolder = tempname;
mkdir(downloadFolder);
zipFile = fullfile(downloadFolder, 'Tide-Model-Driver.zip');
sourceFolder = fullfile(downloadFolder, ...
    ['Tide-Model-Driver-' tmdCommit]);

fprintf('Downloading TMD functions...\n');
websave(zipFile, downloadURL);

fprintf('Extracting TMD functions...\n');
unzip(zipFile, downloadFolder);

if ~exist(fullfile(sourceFolder, 'tmd_predict.m'), 'file')
    error('tmd_predict.m was not found after extracting %s.', zipFile);
end

movefile(sourceFolder, pathTMD);
delete(zipFile);
rmdir(downloadFolder, 's');

if ~exist(pathTMDPredict, 'file')
    error('TMD functions were not installed correctly: %s', pathTMD);
end

fprintf('TMD functions saved:\n%s\n', pathTMD);

end
