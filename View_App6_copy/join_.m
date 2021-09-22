function C = join_(A,~)

% 这个函数是因为旧版MATLAB的join不好使，所以仿照写一个join_
% A 是split后的names数组，B是符号，'\'之类的，注意使用单引号

l = length(A);
C = cell(1,1);
for ii = 1:l
    C{1} = fullfile(C{1},A{ii});
end

end


