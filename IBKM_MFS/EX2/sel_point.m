function coll = sel_point(points, nb)
    % selectRandomPoints 从给定的三维点集中随机选择若干个点
    %
    % 输入参数:
    %   points - N x 3 的矩阵，每一行代表一个三维点
    %   nb     - 要随机选择的点的数量
    %
    % 输出参数:
    %   coll   - 包含随机选择点的 nb x 3 矩阵

    % 获取总点数
    N = size(points, 1);
    
    % 检查输入参数的有效性
    if nb > N
        error('选择的点数 (%d) 不能超过总点数 (%d)', nb, N);
    elseif nb <= 0
        error('选择的点数必须是一个正整数');
    end
    
    % 使用 randsample 随机选择索引
    selectedIndices = randsample(N, nb);
    
    % 根据选定的索引提取点
    coll = points(selectedIndices, :);
end