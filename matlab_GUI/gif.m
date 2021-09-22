%% ×ö³Égif

TiffRoute = uigetdir();
if ~isempty(dir(fullfile(TiffRoute,'\*.tiff')))
    TiffDir = dir(fullfile(TiffRoute,'\*.tiff'));
elseif ~isempty(dir(fullfile(TiffRoute,'\*.tif')))
    TiffDir = dir(fullfile(TiffRoute,'\*.tif'));
end

h = fspecial('average', 3)
for ii = 1:length(TiffDir)
    
    I = imread(fullfile(TiffRoute,TiffDir(ii).name));
    I = imfilter(I,h,'replicate');
    imagesc(imadjust(I))
    axis off
    
    %set(gcf,'units','normalized','position',[0.2 0.2 0.5 0.5]);
    % set(0,'defaultfigurecolor','w');
    drawnow;
    
    F = getframe(gcf);
    I = frame2im(F);
    [I,map] = rgb2ind(I,256);
    
    if ii == 1
        imwrite(I,map,'Scratch.gif','gif','Loopcount',Inf,'DelayTime',0.1);
    else
        imwrite(I,map,'Scratch.gif','gif','WriteMode','append','DelayTime',0.1);
    end
    
end

I0 = imfilter(I,h,'replicate');
subplot(121)
imagesc(I)
subplot(122)
imagesc(I0)
    
    
