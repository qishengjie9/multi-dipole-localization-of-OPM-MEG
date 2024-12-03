%% 绘制相关性的箱线图
clear all
close all
file_data = '5S_corr1';
data = load(['.\data_results\poserr_500result_',file_data]);
dataset = data.poserr_500result_5S_corr1;
%% 数据准备
x = 1:4;
method = 3;
switch method
    case 1
        DLE_range  = [0 25];
        RMSE_range = [0 2.1];
    case 2
        DLE_range  = [0 25];
        RMSE_range = [0 2.1];
    case 3
        DLE_range  = [0 31];
        RMSE_range = [0 2.1];
end
        
% DLE calculation
DLE_dataset = squeeze(dataset(1,method,:,:,:));
DLE_mean = zeros(4,6);
DLE_std  = zeros(4,6);

RMSE_dataset = squeeze(dataset(2,method,:,:,:));
RMSE_mean = zeros(4,6);
RMSE_std  = zeros(4,6);
for i = 1:4 %for four SNR
    DLE_mean(i,:) = mean(squeeze(DLE_dataset(i,:,:)),1);
    DLE_std(i,:)  = std(squeeze(DLE_dataset(i,:,:)));

    RMSE_mean(i,:) = mean(squeeze(RMSE_dataset(i,:,:)),1);
    RMSE_std(i,:)  = std(squeeze(RMSE_dataset(i,:,:)));
end
 
%% 绘制并保存DLE柱状图
GO = bar(x,DLE_mean,1,'EdgeColor','k');
hold on
% 添加误差棒
[M,N] = size(DLE_mean);
xpos = zeros(M,N);
for i = 1:N
    xpos(:,i) = GO(1,i).XEndPoints'; % v2019b
end
AVG = DLE_std/sqrt(500);
hE = errorbar(xpos, DLE_mean, AVG);

% 柱状图赋色
GO(1).FaceColor = [229, 157, 160]/255;
GO(2).FaceColor = [141, 204, 229]/255;
GO(3).FaceColor = [244, 202, 169]/255;
GO(4).FaceColor = [179, 212, 163]/255;
GO(5).FaceColor = [160, 186, 149]/255;
GO(6).FaceColor = [230, 205, 207]/255;
% 误差棒属性
set(hE, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1.2)
set(gca, 'Box', 'off', ...                                         % 边框
         'XGrid', 'off', 'YGrid', 'on', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.01 .01], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'Ylim' , DLE_range,...
         'Xticklabel',{'-10 dB' '0 dB' '10 dB' '20 dB'})% X坐标轴刻度标签
%'Ylim' , [0 7], ...  
%         'Ylim' , [0 7], ... %
% 坐标区调整
% set(gca, 'Box', 'off', ...                                         % 边框
%          'XGrid', 'off', 'YGrid', 'on', ...                        % 网格
%          'TickDir', 'out', 'TickLength', [.01 .01], ...            % 刻度
%          'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
%          'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1],...           % 坐标轴颜色
%          'YTick', 0:0.1:1,...                                      % 刻度位置、间隔
%          'Ylim' , [0 0.5], ...                                     % 坐标轴范围
%          'Xticklabel',{'0 dB' '10 dB' '20dB'},...% X坐标轴刻度标签
%          'Yticklabel',{[0:0.1:1]})                                 % Y坐标轴刻度标签
% Legend 设置
% for gain error
% hLegend = legend([GO(1),GO(2),GO(3),GO(4),GO(5),GO(6)], ...
%                  'No gain error', '2% gain error', '4% gain error','6% gain error', '8% gain error','10% gain error', ...
%                  'Location', 'northeast');
% for crosstalk
% hLegend = legend([GO(1),GO(2),GO(3),GO(4),GO(5),GO(6)], ...
%                  'No crosstalk', '1% crosstalk', '2% crosstalk','3% crosstalk', '4% crosstalk','5% crosstalk', ...
%                  'Location', 'northeast');
%for position error
hLegend = legend([GO(1),GO(2),GO(3),GO(4),GO(5),GO(6)], ...
                 'No error', '1 mm error', '2 mm error','3 mm error', '4 mm error','5 mm error', ...
                 'Location', 'southeast');
%for angular error
% hLegend = legend([GO(1),GO(2),GO(3),GO(4),GO(5),GO(6)], ...
%                  'No error', '2° error', '4° error','6° error', '8° error','10° error', ...
%                  'Location', 'northeast');
hLegend.NumColumns = 2;
% 字体和字号
set(gca, 'FontName', 'Times New Roman')
%set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 14,'FontWeight' , 'bold')
% 背景颜色
set(gcf,'Color',[1 1 1])
hLegend.FontSize = 10;
% 图片输出
% 获取当前图形的句柄
fig = gcf;
% 设置图形的背景颜色为白色（或任何其他颜色），以避免白边
set(fig, 'Color', 'w');
% 设置坐标轴的背景颜色为白色（或任何其他颜色），以避免白边
ax = gca;
set(ax, 'Color', 'w');
% po = get(hLegend,'Position');
% set(hLegend,'position',[po(1),po(2)+0.025,po(3),po(4)]);
% 保存图形为图像文件，设置分辨率和无白边
filename = 'fig_results\single axis\angular error\output_image.png';
print(fig, filename, '-dpng', '-r300', '-painters', '-loose');
% 读取保存的图像
img = imread(filename);
% 调整图像大小为原来的一半
img_resized = imresize(img, 0.5);
% 保存调整后的图像
imwrite(img_resized, ['fig_results\single axis\pose error\',sprintf('DLE_%d_',method),file_data,'.png']);
close all
%% 绘制并保存RMSE柱状图
GO = bar(x,RMSE_mean,1,'EdgeColor','k');
hold on
% 添加误差棒
[M,N] = size(RMSE_mean);
xpos = zeros(M,N);
for i = 1:N
    xpos(:,i) = GO(1,i).XEndPoints'; % v2019b
end
AVG = RMSE_std/sqrt(500);
hE = errorbar(xpos, RMSE_mean, AVG);

% 柱状图赋色
GO(1).FaceColor = [229, 157, 160]/255;
GO(2).FaceColor = [141, 204, 229]/255;
GO(3).FaceColor = [244, 202, 169]/255;
GO(4).FaceColor = [179, 212, 163]/255;
GO(5).FaceColor = [160, 186, 149]/255;
GO(6).FaceColor = [230, 205, 207]/255;
% 误差棒属性
set(hE, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1.2)
set(gca, 'Box', 'off', ...                                         % 边框
         'XGrid', 'off', 'YGrid', 'on', ...                        % 网格
         'TickDir', 'out', 'TickLength', [.01 .01], ...            % 刻度
         'XMinorTick', 'off', 'YMinorTick', 'off', ...             % 小刻度
         'Ylim' , RMSE_range,...
         'Xticklabel',{'-10 dB' '0 dB' '10 dB' '20 dB'})% X坐标轴刻度标签
% for gain error   
% hLegend = legend([GO(1),GO(2),GO(3),GO(4),GO(5),GO(6)], ...
%                  'No gain error', '2% gain error', '4% gain error','6% gain error', '8% gain error','10% gain error', ...
%                  'Location', 'northeast');
% for crosstalk
% hLegend = legend([GO(1),GO(2),GO(3),GO(4),GO(5),GO(6)], ...
%                  'No crosstalk', '1% crosstalk', '2% crosstalk','3% crosstalk', '4% crosstalk','5% crosstalk',...
%                  'Location', 'northeast');
%for position error
hLegend = legend([GO(1),GO(2),GO(3),GO(4),GO(5),GO(6)], ...
                 'No error', '1 mm error', '2 mm error','3 mm error', '4 mm error','5 mm error', ...
                 'Location', 'southeast');                         
%for angular error
% hLegend = legend([GO(1),GO(2),GO(3),GO(4),GO(5),GO(6)], ...
%                  'No error', '2° error', '4° error','6° error', '8° error','10° error', ...
%                  'Location', 'northeast');
hLegend.NumColumns = 2;
% 字体和字号
set(gca, 'FontName', 'Times New Roman')
%set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 14,'FontWeight' , 'bold')
% 背景颜色
set(gcf,'Color',[1 1 1])
hLegend.FontSize = 10;
% po = get(hLegend,'Position');
% set(hLegend,'position',[po(1),po(2)+0.02,po(3),po(4)]);
% 图片输出
% 获取当前图形的句柄
fig = gcf;
set(fig, 'Color', 'w');
ax = gca;
set(ax, 'Color', 'w');
% 保存图形为图像文件，设置分辨率和无白边
filename = 'fig_results\single axis\angular error\output_image.png';
print(fig, filename, '-dpng', '-r300', '-painters', '-loose');
% 读取保存的图像
img = imread(filename);
% 调整图像大小为原来的一半
img_resized = imresize(img, 0.5);
% 保存调整后的图像
imwrite(img_resized, ['fig_results\single axis\pose error\',sprintf('RMSE_%d_',method),file_data,'.png']);
% P = hLegend.Position;
% hLegend.Position = P + [0.015 0.03 0 0];