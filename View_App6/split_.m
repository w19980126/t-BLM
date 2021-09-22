function names = split_(folder,symbol)

% 本函数是因为旧版MATLAB里面split函数没法使用，因此盗版了一个拆分函数，功能比较
% 局限，仅适用于将folder里面的分隔符拆分掉

l = length(folder);
loc = find(folder == symbol);
names = cell(length(loc) + 1,1);
for ii = 1:length(names)
    if ii == 1
        names{1} = folder(1:loc(1) - 1);
    elseif ii == length(names)
        names{ii} = folder(loc(ii - 1) + 1:end);
    else
        names{ii} = folder(loc(ii - 1) + 1:loc(ii) - 1);
    end
end

end