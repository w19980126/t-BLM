function my_raw2tiff(filepath)

filelist = dir([filepath '\*.raw']);
[folder_structure,current_folder] = fileparts(filepath);

if ~isempty(filelist)
    if exist([folder_structure '\TIFF']) ==0
        mkdir([folder_structure '\TIFF']);
    end
    if exist([folder_structure '\TIFF\' current_folder  ]) == 0
        mkdir([folder_structure '\TIFF\' current_folder  ]);
    end
end

eight_bit = 0; 
xres = 640;
yres = 480;
for file = 1:length(filelist)
    fid = fopen([filepath '\' filelist(file).name]);
    A = fread(fid, 'uint8=>uint8');
    fclose(fid);
    E = double(A(1:2:end));
    F = double(A(2:2:end));
    G = 64*E+F/4;
    if eight_bit == 1
        intensity = reshape(E, [xres yres]);
    elseif eight_bit == 0
        intensity = reshape(G, [xres yres])';
    end

    imwrite(uint16(intensity), [folder_structure '\TIFF\' current_folder '\' filelist(file).name '.tiff'], 'Compression', 'none');
    
end

end