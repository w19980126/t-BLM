%% EIS analysis
% import frequency,Z_real and Z_imag from .csv data to .txt data for Zview
% analysing

filepath = uigetdir();
files = dir(fullfile(filepath,'*.csv'));
for ii = 1:length(files)
    
    temp = importdata(fullfile(filepath,files(ii).name));
    if isfield(temp,'colheaders')
        if length(temp.colheaders) == 5
            Freq = temp.data(:,1);
            Z_real = temp.data(:,2);
            Z_imag = temp.data(:,3);
            saveroute = fullfile(filepath,[files(ii).name '.txt']);
            Result_Mat = [Freq Z_real Z_imag];
            save(saveroute,'Result_Mat','-ascii');
        end
    end
    
end

%% LSV数据

filepath = [];
filepath = uigetdir(filepath);
datapath = filepath;
files = dir(fullfile(filepath,'*.csv'));
kk = 0;
for ii = 1:length(files)
    
    temp = importdata(fullfile(filepath,files(ii).name));
    if isfield(temp,'colheaders')
        if length(temp.colheaders) == 3
            kk = kk+1;
            LSV.potential(:,kk) = temp.data(:,1);
            LSV.current(:,kk) = temp.data(:,2);
        end
    end
    
end

figure
hold on
for ii = 1:size(LSV.potential,2)
    plot(LSV.potential(:,ii),LSV.current(:,ii),'linewidth',3);
end
legend
ha = findobj(gcf,'type','axes');
set(ha,'fontsize',15,'fontweight','bold','titlefontweight','bold');
title('食人鱼溶液处理17min后的金片');
xlabel('Potential (V vs. Ag/AgCl)');
ylabel('Current (A)')
box on

%% EIS作图

filepath = uigetdir();
files = dir(fullfile(filepath,'*.txt'));
figure
hold on
for ii = 1:length(files)
    temp = importdata(fullfile(filepath,files(ii).name));
    x = temp(:,2);
    y = -temp(:,3);
    plot(x,y,'linewidth',2)
end
hold off
box on
xlim([0,20000])




