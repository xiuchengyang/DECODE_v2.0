function tideModel = DownloadEOT20Model(pathMain)

if nargin < 1 || isempty(pathMain)
    pathMain = fileparts(mfilename('fullpath'));
end

modelFolder = fullfile(pathMain, 'TMD_functions', ...
    'TideModelData', 'EOT20_ocean.nc');
tideModel = fullfile(modelFolder, 'EOT20_ocean.nc');

if exist(tideModel, 'file')
    ncinfo(tideModel);
    fprintf('EOT20 ocean tide model found:\n%s\n', tideModel);
    return;
end

if ~exist(modelFolder, 'dir')
    mkdir(modelFolder);
end

downloadURL = ...
    'https://www.chadagreene.com/tide_data/EOT20_ocean.nc.zip';
zipFile = fullfile(modelFolder, 'EOT20_ocean.nc.zip');

fprintf('Downloading EOT20 ocean tide model...\n');
websave(zipFile, downloadURL);

fprintf('Extracting EOT20 ocean tide model...\n');
unzip(zipFile, modelFolder);
delete(zipFile);

if ~exist(tideModel, 'file')
    modelFiles = dir(fullfile(modelFolder, '**', 'EOT20_ocean.nc'));

    if numel(modelFiles) ~= 1
        error('EOT20_ocean.nc was not found after extracting %s.', zipFile);
    end

    downloadedModel = fullfile(modelFiles(1).folder, modelFiles(1).name);
    movefile(downloadedModel, tideModel);
end

ncinfo(tideModel);
fprintf('EOT20 ocean tide model saved:\n%s\n', tideModel);

end
