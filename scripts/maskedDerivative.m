function D = maskedDerivative(F,G,valid,dim,h,edgeX,edgeY)
% 对 ny×nx×3 的 F 求一阶导数，邻点必须属于同一晶粒。
% dim=2 对列求导（x 方向）；dim=1 对行求导（y 方向）。
    % 将输出导数初始化为 NaN，表示尚不能计算。
    D = nan(size(F));
    % 读入图像的行数 ny 和列数 nx。
    [ny,nx,~] = size(F);
    % 遍历每一行像素。
    for i = 1:ny
        % 遍历该行的每一列像素。
        for j = 1:nx
            % 无效像素或无晶粒号像素直接跳过。
            if ~valid(i,j) || G(i,j)<=0 || ~all(isfinite(reshape(F(i,j,:),1,[]))), continue; end
            % 当 dim=2 时，沿列方向即 x 方向计算差分。
            if dim == 2
                % 检查左邻点是否存在、有效、属于同一晶粒且三个场分量都有限。
                hasM = j>1  && edgeX(i,j-1) && valid(i,j-1) && G(i,j-1)==G(i,j) && all(isfinite(reshape(F(i,j-1,:),1,[])));
                % 以同样条件检查右邻点。
                hasP = j<nx && edgeX(i,j) && valid(i,j+1) && G(i,j+1)==G(i,j) && all(isfinite(reshape(F(i,j+1,:),1,[])));
                % 左右同晶粒邻点都可用时采用中心差分。
                if hasM && hasP
                    % 用左右邻点之差除以 2h。
                    D(i,j,:) = (F(i,j+1,:)-F(i,j-1,:))/(2*h);
                % 只有右邻点可用时改用前向差分。
                elseif hasP
                    % 用右邻点与当前点之差除以 h。
                    D(i,j,:) = (F(i,j+1,:)-F(i,j,:))/h;
                % 只有左邻点可用时改用后向差分。
                elseif hasM
                    % 用当前点与左邻点之差除以 h。
                    D(i,j,:) = (F(i,j,:)-F(i,j-1,:))/h;
                % 结束 x 方向邻点选择；无合适邻点则保留 NaN。
                end
            % 当 dim=1 时，沿行方向即 y 方向计算差分。
            elseif dim == 1
                % 检查上邻点是否存在、有效、同晶粒且场分量有限。
                hasM = i>1  && edgeY(i-1,j) && valid(i-1,j) && G(i-1,j)==G(i,j) && all(isfinite(reshape(F(i-1,j,:),1,[])));
                % 以同样条件检查下邻点。
                hasP = i<ny && edgeY(i,j) && valid(i+1,j) && G(i+1,j)==G(i,j) && all(isfinite(reshape(F(i+1,j,:),1,[])));
                % 上下同晶粒邻点均可用时采用中心差分。
                if hasM && hasP
                    % 用下邻点与上邻点之差除以 2h。
                    D(i,j,:) = (F(i+1,j,:)-F(i-1,j,:))/(2*h);
                % 只有下邻点可用时采用前向差分。
                elseif hasP
                    % 用下邻点与当前点之差除以 h。
                    D(i,j,:) = (F(i+1,j,:)-F(i,j,:))/h;
                % 只有上邻点可用时采用后向差分。
                elseif hasM
                    % 用当前点与上邻点之差除以 h。
                    D(i,j,:) = (F(i,j,:)-F(i-1,j,:))/h;
                % 结束 y 方向邻点选择；无可用邻点则保留 NaN。
                end
            % 当 dim 既不是 1 也不是 2 时进入非法参数分支。
            else
                % 说明合法方向参数只有 1 与 2。
                error('dim must be 1 (y) or 2 (x).');
            % 结束方向判断。
            end
        % 结束列循环。
        end
    % 结束行循环。
    end
% 结束同晶粒差分函数。
end

% 定义晶界绘图函数：输入晶粒号网格、坐标向量和两个步长。
