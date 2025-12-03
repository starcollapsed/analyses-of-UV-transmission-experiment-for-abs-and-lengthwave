clear
clc

%需要自己修改的参数
curvature_threshold = 0.3;
dataname = 'ZnO_nano.txt';

%读数(记得改成自己文件名）
data = readtable(dataname, 'Delimiter', ' \t');
wl = data.nm;
Abs = data.Abs;

E = 1.6e-19;%元电荷
h = 6.63*10^-34;%普朗克常数
hv = (h*3*10^17/E)./wl;
a = (Abs.*hv).^2;

%作图
subplot(2,1,1)
plot(wl,Abs,'-.')
title('截线法与Tauc plot法')
xlabel('wavelength (nm) ')
ylabel('abs')
xlim([min(wl),max(wl)])

subplot(2,1,2)
plot(hv,a,'-')
xlabel('engergy (eV)')
ylabel('(\alphah\nu)^2')
xlim([0,max(hv)])

%截线法
grad = gradient(Abs,wl);
i=1;
while grad(i)<=0
    i=i+1;
end
[grad1,l1] = min(grad(i:length(grad)));

%Tauc plot法
dy = gradient(a,hv);
d2y = gradient(dy,hv);
denominator = (1+dy.^2).^(3/2); %去除分母为0的情况
denominator(denominator == 0) = eps;
k = abs(d2y)./denominator;

is_linear = k < curvature_threshold; %提取符合线性区域
linear_segments = {};
current_segment = [];

for i = 1:length(is_linear)
    if is_linear(i)
        current_segment = [current_segment, i];
    elseif ~isempty(current_segment)
        if length(current_segment) >= 10  % 最小长度要求
            linear_segments{end+1} = current_segment;
        end
        current_segment = [];
    end
end
if ~isempty(current_segment) && length(current_segment) >= 10 %防止以线性结尾而不被计入
    linear_segments{end+1} = current_segment;
end

if ~isempty(linear_segments) %找出最佳线性区域

    segment_lengths = cellfun(@length, linear_segments); % 按长度排序（降序），找出长度前三的段
    [sorted_lengths, length_order] = sort(segment_lengths, 'descend');
    
    num_top_segments = min(3, length(linear_segments));% 选择长度前三的段（如果不足三个，则选择全部）
    top_length_segments = length_order(1:num_top_segments);
    
    avg_curvatures_top = zeros(1, num_top_segments); % 计算前三段的平均曲率
    for j = 1:num_top_segments
        seg_idx = top_length_segments(j);
        segment_indices = linear_segments{seg_idx};
        avg_curvatures_top(j) = mean(k(segment_indices));
    end
    
    % 在这前三段中，找出平均曲率最低的一段
    [min_curvature, min_idx] = min(avg_curvatures_top);
    best_seg_idx = top_length_segments(min_idx);
    best_segment = linear_segments{best_seg_idx};

    x_linear = hv(best_segment);%拟合
    y_linear = a(best_segment);
    [coefficients, ~] = polyfit(x_linear, y_linear, 1);
    y_fit = polyval(coefficients, x_linear);
    x_poly = linspace(0,hv(length(hv)),10*length(hv));
    y_poly = polyval(coefficients, x_poly);

    subplot(2,1,2)
    hold on
    plot(x_poly, y_poly,'--')
    legend('实验数据','线性拟合方程')
    hold off

else
    warning('未找到足够长的线性段');

end

%显示结果
disp(['截线法结果：带隙为',num2str(h*3*10^8*10^9/(l1*E)),'eV，对应波长为',num2str(l1),'nm'])
disp(['Tauc plot法结果：带隙为',num2str( -coefficients(2) / coefficients(1)),'eV'])
fprintf('最佳线性区间: %.3f-%.3f eV\n', x_linear(1), x_linear(end));