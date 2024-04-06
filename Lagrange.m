%6889 
[m,n]=size(one);
% 使用 unique 函数找到唯一值列表
% 示例数据

onecopy=one;
% 使用unique函数找到唯一值列表
unique_values = unique(one(:,5));

% 显示唯一值列表
disp('列向量的所有取值：');
disp(unique_values);

for i=1:m
    if strcmp(onecopy(i,5),'轻度污染')
        onecopy(i,5)=4;
    end
    if strcmp(onecopy(i,5),'中度污染')
        onecopy(i,5)=3;
    end
    if  strcmp(onecopy(i,5),'重度污染')
        onecopy(i,5)=2;
    end
    if strcmp(onecopy(i,5),'严重污染')
        onecopy(i,5)=1;
    end
    if strcmp(onecopy(i,5),'')
        onecopy(i,5)=0;
    end
    if strcmp(onecopy(i,5),'无')
        onecopy(i,5)=0;
    end
    if  strcmp(onecopy(i,5),'良')
        onecopy(i,5)=5;
    end
    if  strcmp(onecopy(i,5),'优')
        onecopy(i,5)=6;
    end
end
for i=1:m
    for j=1:n
        if strcmp(onecopy(i,j),'')
            onecopy(i,j)=0;
        end
    end
    if strcmp(onecopy(i,8),'1')
        onecopy(i,8)=0;
    end
end
oneshu=str2double(onecopy);
oneshu(oneshu == 0) = NaN;

% 使用多重插补法填补缺失值
num_imputations = 5; % 设置插补次数
matrix_imputed = NaN(m, n, num_imputations); % 存储插补后的矩阵


% 计算每一列的均值
mean_values = mean(oneshu, 'omitnan'); % 计算每一列的均值，忽略 NaN 值

for i = 1:num_imputations
    % 使用内置函数填补缺失值（此处仅简单地用均值填补）
   matrix_imputed(:, :, i) = fillmissing(oneshu, 'constant', mean_values);
end

% 计算插补后的均值和标准差
mean_imputed = mean(matrix_imputed, 3);
std_imputed = std(matrix_imputed, 0, 3, 'omitnan');

% 显示插补后的均值和标准差
disp('插补后的均值：');
disp(mean_imputed);
disp('插补后的标准差：');
disp(std_imputed);

twocope=str2double(two);
[m2,n2]=size(two)
twocope(twocope==999990)=NaN;
twocope(twocope==999999)=NaN;

twomiddle=twocope(:,1:5)
%_______________________________________________________________拉格朗日函数matlab——————————————————————————————————————————————————————————————————————
%输入x数据点x坐标（后期是去掉缺失点），y数据点y坐标（后期是去掉缺失点），xo（缺失点x坐标），lie插补列的序号，原矩阵twocope，
% 假设原始数据向量
lie=7
data = twocope(:,lie);


% 找到缺失值的位置
missing_idx = find(isnan(data));

%输入数据点的坐标
x = 1:m2
y = data
%  

valid_indices = ~isnan(y);  % 找到 y 中不是 NaN 的元素的索引
x = x(valid_indices); % 根据索引从 x 中选取相应的元素
y = y(valid_indices); % 根据索引从 y 中选取相应的元素
plot(x,y,'o')
%输入所要估值的插值点的x坐标向量
x0 = missing_idx
%  
% y0=Lagrange(x0,x,y);


% 找到缺失值的位置
missing_idx = x0';

% 获取已知数据点的坐标和值
x_known = x;
y_known = y';

% 拉格朗日插值
for i = 1:length(missing_idx)
    % 当前缺失值的索引
    idx = missing_idx(i);
    
    % 当前缺失值的 x 坐标
    x0 = idx;
    
    % 计算拉格朗日插值
    y0 = Lagrange(x_known, y_known, x0);
    
    % 将插值结果填补到缺失值位置,插补的列位
    twocope(idx,lie) = y0;
end


% 拉格朗日插值函数
% 拉格朗日插值函数（使用前后十个数进行插值）
function y0 = Lagrange(x, y, x0)
    n = length(x);
    y0 = 0;
    
    % 确定插值范围
    range = 10;
    % 寻找最近的两个点的索引
    [~, idx] = sort(abs(x - x0));
    closest_indices = idx(1:range);
    
    for i = closest_indices
        % 计算每个基函数的值
        L = 1;
        for j = closest_indices
            if j ~= i
                L = L * (x0 - x(j)) / (x(i) - x(j));
            end
        end
        % 累加基函数乘以对应的 y 值
        y0 = y0 + L * y(i);
    end
end
