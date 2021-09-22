filepath = uigetdir();
files = dir(filepath);
files(1:2) = [];
for ii = 1:length(files)
    temp = dir(fullfile(filepath,files(ii).name,'*.raw'));
    if ~isempty(temp)
        my_raw2tiff(fullfile(filepath,files(ii).name));
    end
end
        
