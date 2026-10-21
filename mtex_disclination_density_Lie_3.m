%% MTEX：基于李代数方法计算二维 EBSD 向错密度
% -------------------------------------------------------------------------
% 本脚本按下述论文中的主要数值步骤计算：
% C. Du et al., Acta Materialia 257 (2023) 119176,
% "A Lie-algebra-based description of disclination densities ..."
%
% 输出变量：
%   theta13、theta23、theta33：二维 EBSD 可计算的向错密度张量分量；
%                             单位为 rad/(输入长度单位)^2。
%   rhoTheta：上述三个分量的平方和开根号。
%   CSV、MAT 和 PNG 文件保存在原始 EBSD 文件旁边。
%
% 适用范围与解释：
%   1. 要求二维 EBSD 坐标落在轴对齐的规则矩形网格上；六角采样、
%      显著畸变的网格和散点数据会报错，不会被自动重采样。
%   2. EBSD 测得晶格取向，不直接测得独立的弹性旋转场；以下结果
%      是按照论文中取向场处理方式得到的模型估计值。
%   3. 此版本直接使用文件中保存的取向代表元，未做晶体对称等价
%      取向的空间连续性优化。代表元突然跳变可产生伪影，需要检查图像。
%   4. 晶粒内部优先用中心差分；晶界使用指向同一晶粒内部的单侧差分。
%   5. 脚本采用的李代数修正公式：
%       kappa_x = d_x(omega) + 1/2 omega x d_x(omega)
%       kappa_y = d_y(omega) + 1/2 omega x d_y(omega)
%       theta(:,3) = D~_x(kappa_y) - D~_y(kappa_x)
%       D~_j(v) = d_j(v) + 1/2 omega x d_j(v)
%      其中 x 表示叉乘；omega 为弧度，坐标单位取决于输入文件，
%      通常为微米。
%
% 运行要求：MATLAB 已安装 MTEX，且 MTEX 所在路径已加入 MATLAB 路径。
% 这里只进行了静态检查和输入坐标验证，未在 MATLAB/MTEX 中完整运行；
% 定量发表前应进一步用参考数据核查计算方法和结果。
% -------------------------------------------------------------------------

% 清理工作区变量、关闭现有图窗，并清空命令窗口。
clearvars; close all; clc;

%% 用户设置
% 可以填写 MTEX 支持的 EBSD 文件路径，例如 .ctf、.ang、.osc 或 .h5。
% EBSD.load 通常根据扩展名选择导入接口。
% 输入文件路径留空时，运行过程中会弹出文件选择窗口。
inputFile = '';

% 只分析一个物相：取向矩阵转换要求数据具有同一晶体对称性。
% 此处需填写 MTEX 显示的准确物相名称。
% 只分析 MTEX 中名称为 Magnesium 的相，避免混合晶体对称性的取向矩阵报错。
phaseName = 'Magnesium';

% Oxford .h5oina 可能同时保存清理后的 Data Processing 图和原始 EBSD 图。
% 原始 EBSD 图包含明确的测点 x、y 坐标。
% 对于 Oxford h5oina 文件，选择原始 EBSD 数据集。
h5oinaDataSet = 'EBSD';

% 相邻测点的取向差达到此角度时，晶粒重建可能将它们分在不同晶粒。
% 用 3° 的相邻取向差作为重建晶粒的晶界阈值。
grainBoundaryAngleDeg = 3;

% 脚本默认不删除小晶粒，以免改变原始的小尺度结构。

% true 表示额外保存 EBSD 取向与晶界诊断图。
% 决定是否另存一张取向诊断图。
makeIPFDiagnostic = true;

%% 选择并导入 EBSD 数据
% 如果没有预先填入路径，则进入交互式选文件流程。
if isempty(inputFile)
    % 弹出文件选择框；以下多行列出可选扩展名、类型名称与窗口标题。
    [fn, fp] = uigetfile( ...
        {'*.ctf;*.ang;*.osc;*.h5;*.h5oina;*.oh5;*.dream3d;*.crc;*.cpr', ...
         'MTEX-supported EBSD files (*.ctf, *.ang, *.osc, *.h5, etc.)'; ...
         '*.*', 'All files'}, ...
        'Select an EBSD file');
    % 检查用户是否在文件对话框中按了“取消”。
    if isequal(fn,0)
        % 未选文件时立刻报错，阻止后续计算。
        error('No EBSD file selected.');
    % 结束取消选择的判断。
    end
    % 将所选文件夹和文件名拼成完整路径。
    inputFile = fullfile(fp,fn);
% 结束输入路径为空时的选文件分支。
end

% 核实输入文件确实存在。
if ~exist(inputFile,'file')
    % 路径无效时在错误信息中列出该路径。
    error('Input file does not exist: %s',inputFile);
% 结束文件存在性检查。
end

% 在命令窗口打印正在导入的文件路径。
fprintf('Loading EBSD file:\n  %s\n',inputFile);
% 尝试执行 MTEX 导入；若失败则由后面的 catch 提供说明。
try
    % 从路径中提取文件扩展名，前两个输出用波浪号舍弃。
    [~,~,inputExt] = fileparts(inputFile);
    % 判断是否为 Oxford 的 .h5oina 格式。
    if strcmpi(inputExt,'.h5oina')
        % 从文件中明确读取用户指定的 EBSD 数据集。
        ebsdAll = EBSD.load(inputFile,'dataSet',h5oinaDataSet);
    % 对于其他文件格式，进入 MTEX 自动导入分支。
    else
        % 调用 MTEX 的 EBSD.load 自动识别并读取文件。
        ebsdAll = EBSD.load(inputFile);
    % 结束格式分支。
    end
% 捕获导入时抛出的异常，并把异常对象命名为 ME。
catch ME
    % 把 MTEX 原始错误附在提示后面；这条 error 语句跨三行。
    error(['MTEX could not import this file automatically. For text/CSV/Excel ' ...
        'files, you may need to specify ColumnNames, phase columns, and crystal ' ...
        'symmetry in EBSD.load. Original error: %s'], ME.message);
% 结束导入异常处理。
end

% 去掉未索引测点，只用已索引测点进行晶粒重建。
% 求导时还会限制差分邻点处在同一晶粒。
% 尝试筛选已经成功索引的 EBSD 测点。
try
    % MTEX 的 indexed 选择器会排除未索引点。
    ebsd = ebsdAll('indexed');
% 如果此版本 MTEX 无法筛选 indexed，则进入异常分支。
catch
    % 提示检查 MTEX 版本和 EBSD 数据。
    error('Could not select indexed EBSD measurements. Check your MTEX version and input data.');
% 结束 indexed 筛选的异常处理。
end
% 若 phaseName 不为空，则进一步筛选目标物相。
if ~isempty(phaseName)
    % 尝试按物相名称从已索引测点中选择单一物相。
    try
        % 仅保留目标物相的测点和相应取向。
        ebsd = ebsd(phaseName);
    % 捕获物相名称不存在或筛选失败的异常。
    catch ME
        % 报告目标物相名称和 MTEX 的原始错误；下一行传入参数。
        error('Could not select phase "%s". Check the mineral name in the MTEX phase list. Details: %s', ...
            phaseName, ME.message);
    % 结束物相筛选异常处理。
    end
% 结束是否选择物相的判断。
end
% 检查筛选后是否仍有足够的测点。
if numel(ebsd) < 9
    % 少于九个测点时停止二维导数计算。
    error('Too few indexed measurements for a 2-D derivative calculation.');
% 结束测点数量检查。
end

%% 重建晶粒并取得各测点的晶粒编号
% 打印晶粒重建阈值，便于核对输入参数。
fprintf('Reconstructing grains (threshold %.2f deg)...\n',grainBoundaryAngleDeg);
% 尝试根据所选物相的取向重建晶粒。
try
    % calcGrains 返回晶粒集合 grains，以及带 grainId 的 EBSD 测点。
    [grains, ebsd] = calcGrains(ebsd,'angle',grainBoundaryAngleDeg*degree);
% 捕获晶粒重建失败的异常。
catch ME
    % 说明可能是晶体对称性或 MTEX 设置问题。
    error('calcGrains failed. Check the imported crystal symmetry / MTEX setup. Details: %s',ME.message);
% 结束晶粒重建异常处理。
end
% 优先尝试从 ebsd 的直接属性读取晶粒编号。
try
    % 将晶粒编号转成列向量 double，便于随后按测点编号索引。
    grainId = double(ebsd.grainId(:));
% 如果当前 MTEX 版本不支持直接访问，使用备用属性路径。
catch
    % 尝试从 ebsd.prop 中取得晶粒编号。
    try
        % 读取并转换该版本 MTEX 保存的 grainId。
        grainId = double(ebsd.prop.grainId(:));
    % 捕获两种 grainId 读取方式均失败的异常。
    catch ME
        % 输出最后一次读取失败的原因。
        error('Could not read grainId after calcGrains. Details: %s',ME.message);
    % 结束备用读取的异常处理。
    end
% 结束晶粒编号读取的异常处理。
end

%% 根据测点坐标构造轴对齐的规则矩形网格
% 读取每个目标物相测点的 x 坐标并转为 double 列向量。
x = double(ebsd.x(:));
% 读取对应的 y 坐标并转为 double 列向量。
y = double(ebsd.y(:));
% 保存筛选后测点的总数 N。
N = numel(ebsd);
% 检查 x、y 坐标的数量是否都等于测点数。
if numel(x) ~= N || numel(y) ~= N
    % 坐标长度不一致时停止处理，以免错位。
    error('EBSD x/y coordinate vectors are inconsistent with the measurement count.');
% 结束坐标数量检查。
end

% 此处容差仅用于合并浮点表示造成的微小坐标波动。
% 取 x 跨度、y 跨度和 1 的最大值，作为坐标容差的尺度。
coordScale = max([max(x)-min(x), max(y)-min(y), 1]);
% 设置只用于合并浮点抖动的坐标量化容差。
coordTol = coordScale*1e-8;
% 将每个 x 坐标映射到量化整数键。
xKey = round((x-min(x))/coordTol);
% 以相同方式量化 y 坐标。
yKey = round((y-min(y))/coordTol);
% 提取有序的 x 列位置，并换算回坐标单位。
xVals = unique(xKey,'sorted')*coordTol + min(x);
% 提取有序的 y 行位置，并换算回坐标单位。
yVals = unique(yKey,'sorted')*coordTol + min(y);
% 强制 xVals 和 yVals 为列向量。
xVals = xVals(:); yVals = yVals(:);

% 判断 x 或 y 方向是否少于三个不同采样位置。
if numel(xVals)<3 || numel(yVals)<3
    % 位置过少则无法进行稳定的二维差分。
    error('At least 3 unique x and y coordinates are required.');
% 结束采样位置数量检查。
end

% 用相邻距离的典型值估计整数网格索引。
% 若直接取量化间距的最小值，极小的舍入误差会跨许多列累积，
% 使本来规则的 EBSD 图被误判为非规则网格。
% 计算相邻 x 列距离，排除低于容差阈值的数值噪声。
dxCandidates = diff(xVals); dxCandidates = dxCandidates(dxCandidates > coordTol*10);
% 计算相邻 y 行距离，并排除数值噪声。
dyCandidates = diff(yVals); dyCandidates = dyCandidates(dyCandidates > coordTol*10);
% 检查任一方向是否找不到有效距离候选值。
if isempty(dxCandidates) || isempty(dyCandidates)
    % 无法估计步长时报告错误。
    error('Could not determine the scan pitch from x/y coordinates.');
% 结束候选距离检查。
end
% 用候选距离的中位数初始化 x、y 步长，避免最小值的舍入偏差。
dxSeed = median(dxCandidates); dySeed = median(dyCandidates);
% 检查初始步长是否有限且严格大于零。
if ~isfinite(dxSeed) || ~isfinite(dySeed) || dxSeed<=0 || dySeed<=0
    % 无效步长导致后续网格映射不可用。
    error('Could not determine the scan pitch from x/y coordinates.');
% 结束初始步长检查。
end

% 用实测最小坐标作为网格原点。
x0 = min(x); y0 = min(y);
% 根据初始 x 步长把测点量化为列号，并计算总列数。
nx = max(round((x-x0)/dxSeed))+1;
% 根据初始 y 步长把测点量化为行号，并计算总行数。
ny = max(round((y-y0)/dySeed))+1;
% 根据真实坐标端点修正步长，消除前面量化造成的累计误差。
% 用 x 坐标的真实总跨度除以列间隔数，校正最终 dx。
dx = (max(x)-x0)/(nx-1);
% 用 y 坐标的真实总跨度除以行间隔数，校正最终 dy。
dy = (max(y)-y0)/(ny-1);
% 按校正后的 dx 计算 MATLAB 从 1 开始的列索引 ix。
ix = round((x-x0)/dx)+1;
% 按校正后的 dy 计算 MATLAB 从 1 开始的行索引 iy。
iy = round((y-y0)/dy)+1;

% 每个矩形网格测点都必须足够接近由原点和步长构造的理想网格点。
% 计算每个测点 x 坐标相对理想矩形网格的绝对残差。
errX = abs(x-(x0+(ix-1)*dx));
% 计算每个测点 y 坐标相对理想矩形网格的绝对残差。
errY = abs(y-(y0+(iy-1)*dy));
% 如果任何坐标残差超过容差或步长的千分之一，则拒绝该网格。
if any(errX > max(coordTol*20,dx*1e-3)) || any(errY > max(coordTol*20,dy*1e-3))
    % 说明目前的导数程序只接受轴对齐的矩形采样点；错误文本跨五行。
    error(['This scan is not a regular axis-aligned square/rectangular grid. ' ...
        'The current script deliberately does not resample hexagonal or ' ...
        'distorted grids. Please export/reconstruct a regular grid or adapt ' ...
        'the derivative routine for your scan geometry. ' ...
        'Maximum coordinate residuals: x=%.6g, y=%.6g.'],max(errX),max(errY));
% 结束网格残差检查。
end
% 再次确认行列索引都落在生成的网格范围内。
if any(ix<1 | ix>nx | iy<1 | iy>ny)
    % 索引越界说明内部映射发生错误。
    error('Internal grid mapping failed.');
% 结束索引范围检查。
end

% 拒绝多个测点映射到同一像素，避免后一个取向无声覆盖前一个。
% 将每个测点的行、列索引转换为二维数组的线性索引。
linIdx = sub2ind([ny,nx],iy,ix);
% 检查是否有两个测点落入同一网格像素。
if numel(unique(linIdx)) ~= numel(linIdx)
    % 存在重复像素就报错，避免后写的取向覆盖先写的取向。
    error('Duplicate EBSD measurements share the same x/y coordinate.');
% 结束重复测点检查。
end

% 打印网格行列数以及由坐标估计得到的 dx、dy。
fprintf('Grid: %d columns x %d rows; dx = %.6g, dy = %.6g (input length unit).\n',nx,ny,dx,dy);

%% 将旋转矩阵经主值对数映射为旋转向量 omega
% MTEX 的 matrix(orientation) 预期返回 3×3×N 的数组。
% 尝试把 MTEX 取向对象转换成旋转矩阵数组。
try
    % 预期得到维度为 3×3×N 的取向矩阵。
    Rall = matrix(ebsd.orientations);
% 捕获取向矩阵转换失败的异常。
catch ME
    % 连同 MTEX 错误一起报告，方便定位物相或取向问题。
    error('Could not extract orientation matrices using matrix(ebsd.orientations): %s',ME.message);
% 结束矩阵转换异常处理。
end
% 逐一核查矩阵数组的三维长度。
if size(Rall,1) ~= 3 || size(Rall,2) ~= 3 || size(Rall,3) ~= N
    % 不是 3×3×N 时停止，以免错误地给测点分配取向。
    error('Unexpected MTEX matrix(orientations) dimensions. Expected 3x3xN.');
% 结束矩阵维度核查。
end

% 为每个网格像素预分配三个旋转向量分量，默认 NaN 表示无值。
omegaGrid = nan(ny,nx,3);
% 预分配晶粒编号图；零表示没有可用晶粒。
grainGrid = zeros(ny,nx);
% 预分配有效点掩膜，初始全部为 false。
validGrid = false(ny,nx);
% 逐点把 MTEX 测量列表映射到二维网格。
for p = 1:N
    % 跳过没有有限正数晶粒编号的测点。
    if ~isfinite(grainId(p)) || grainId(p) <= 0
        % 继续处理下一个测点。
        continue;
    % 结束无效晶粒编号判断。
    end
    % 取出第 p 个测点的 3×3 旋转矩阵。
    Rp = double(Rall(:,:,p));
    % 检测矩阵是否含 NaN 或 Inf。
    if any(~isfinite(Rp(:)))
        % 矩阵不完整时跳过此点。
        continue;
    % 结束矩阵有效性判断。
    end
    % 将略偏离 SO(3) 的数值矩阵投影回最近的合法旋转矩阵。
    % 用奇异值分解取得最近的正交矩阵因子。
    [U,~,V] = svd(Rp);
    % 构造行列式为 +1 的旋转矩阵，修正数值误差及反射。
    Rp = U*diag([1 1 det(U*V')])*V';
    % 求旋转矩阵的主值对数，并以三分量旋转向量 w 表示。
    w = localRotmLog(Rp);
    % 将三个旋转分量写入对应的网格行列。
    omegaGrid(iy(p),ix(p),:) = reshape(w,1,1,3);
    % 将这个像素所属的晶粒编号写入晶粒图。
    grainGrid(iy(p),ix(p)) = grainId(p);
    % 标记该网格位置为有效测点。
    validGrid(iy(p),ix(p)) = true;
% 结束逐测点循环。
end

% 检查有效网格测点是否仍至少有九个。
if nnz(validGrid) < 9
    % 有效点太少时停止向错密度计算。
    error('Fewer than 9 valid grid locations could be mapped.');
% 结束有效点检查。
end

%% 计算李代数修正的一阶导数与晶格曲率
% 两侧同晶粒邻点均可用时用中心差分；晶界处用指向晶粒内部的单侧差分。
% 对旋转向量场沿 x 列方向求同晶粒内的差分，步长为 dx。
dOmegaX = maskedDerivative(omegaGrid,grainGrid,validGrid,2,dx); % d/dx (columns)
% 对旋转向量场沿 y 行方向求同晶粒内的差分，步长为 dy。
dOmegaY = maskedDerivative(omegaGrid,grainGrid,validGrid,1,dy); % d/dy (rows)

% 三维旋转向量形式下，李括号 [a,b] 对应向量叉乘 a×b。
% 计算 x 向李代数修正曲率：旋转向量梯度加半个叉乘项。
kappaX = dOmegaX + 0.5*crossField(omegaGrid,dOmegaX);
% 以同样公式计算 y 向李代数修正曲率。
kappaY = dOmegaY + 0.5*crossField(omegaGrid,dOmegaY);

%% 计算曲率的李修正旋度，即二维可得的 theta(:,3)
% 按脚本采用的张量旋度约定：theta_i3 = D~_x(kappa_i,y) - D~_y(kappa_i,x)。
% 沿 x 方向对 y 向曲率求受晶粒掩膜约束的导数。
dKappaY_dx = maskedDerivative(kappaY,grainGrid,validGrid,2,dx);
% 沿 y 方向对 x 向曲率求受晶粒掩膜约束的导数。
dKappaX_dy = maskedDerivative(kappaX,grainGrid,validGrid,1,dy);
% 两个曲率导数相减，构成普通二维旋度。
ordinaryCurl = dKappaY_dx - dKappaX_dy;
% 再加上旋转向量与旋度的半个叉乘，得到脚本采用的李修正向错密度列。
theta = ordinaryCurl + 0.5*crossField(omegaGrid,ordinaryCurl);

% 三个分量的单位均为 rad/(输入长度单位)^2。
% 取向错密度张量第三列的第 1 个分量 theta13。
theta13 = theta(:,:,1);
% 取第三列的第 2 个分量 theta23。
theta23 = theta(:,:,2);
% 取第三列的第 3 个分量 theta33。
theta33 = theta(:,:,3);
% 计算三个可测分量的欧氏范数，作为二维标量密度图。
rhoTheta = sqrt(theta13.^2 + theta23.^2 + theta33.^2);
% 把未索引或无效网格位置的标量结果显式设为 NaN。
rhoTheta(~validGrid) = NaN;
% 把未索引位置的 theta13 设为 NaN。
theta13(~validGrid) = NaN;
% 把未索引位置的 theta23 设为 NaN。
theta23(~validGrid) = NaN;
% 把未索引位置的 theta33 设为 NaN。
theta33(~validGrid) = NaN;

% 差分邻点不足而留下 NaN 的像素不会计入有效结果。
% 只有标量结果为有限数值时，才标记为有效计算结果。
validResult = isfinite(rhoTheta);
% 打印可计算密度的像素数以及映射成功的测点数。
fprintf('Finite scalar-density pixels: %d / %d indexed-grid locations.\n',nnz(validResult),nnz(validGrid));
% 如果整张图没有一个有限结果，就进入错误分支。
if nnz(validResult) == 0
    % 说明可能是晶粒太小、晶粒编号过碎或网格条件不满足；错误文本跨两行。
    error(['No finite density values were produced. The map may be too small, ' ...
        'grain IDs may be fragmented, or the scan grid may not be regular.']);
% 结束结果存在性检查。
end

%% 保存数值结果
% 取得原始文件的目录和不带扩展名的文件名。
[outDir,baseName,~] = fileparts(inputFile);
% 如果只给了相对文件名，则把当前工作目录作为输出目录。
if isempty(outDir), outDir = pwd; end
% 设定 CSV 表格的保存路径和文件名。
csvFile = fullfile(outDir,[baseName '_LieDisclinationDensity.csv']);
% 设定 MAT 数据文件的保存路径和文件名。
matFile = fullfile(outDir,[baseName '_LieDisclinationDensity.mat']);
% 设定四宫格 PNG 图的保存路径和文件名。
figFile = fullfile(outDir,[baseName '_LieDisclinationDensity.png']);

% 按网格原点及步长生成每个像素的二维 x、y 坐标。
[Xgrid,Ygrid] = meshgrid(min(x)+(0:nx-1)*dx, min(y)+(0:ny-1)*dy);
% 建立结果表；以下三行依次放入坐标、晶粒号、掩膜、三个分量及范数，并命名表列。
T = table(Xgrid(:),Ygrid(:),grainGrid(:),validGrid(:), ...
    theta13(:),theta23(:),theta33(:),rhoTheta(:), ...
    'VariableNames',{'x','y','grainId','indexed','theta13','theta23','theta33','rhoTheta'});
% 把结果表写入 CSV 文件。
writetable(T,csvFile);
% 把主要结果和中间变量写入 MAT 文件；这条 save 语句跨三行。
save(matFile,'theta13','theta23','theta33','rhoTheta','omegaGrid', ...
    'kappaX','kappaY','grainGrid','validGrid','validResult','dx','dy', ...
    'grainBoundaryAngleDeg','inputFile','xVals','yVals');

%% 绘制并保存向错密度图与取向诊断图
% 创建白底的主结果图窗。
fig = figure('Color','w','Name','Lie-algebra disclination density');
% 把主图排成两行两列，并缩小子图间距。
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
% 计算绘图所需的每个网格列中心 x 坐标。
xPlot = min(x)+(0:nx-1)*dx;
% 计算绘图所需的每个网格行中心 y 坐标。
yPlot = min(y)+(0:ny-1)*dy;

% 切换到第一个子图。
nexttile;
% 显示 theta13 的二维数值图。
imagesc(xPlot,yPlot,theta13);
% 保持当前坐标轴，叠加从 grainGrid 推导的晶界，再结束叠加状态。
hold on; plotRasterBoundaries(grainGrid,xPlot,yPlot,dx,dy); hold off;
% 使用等比例坐标，令 y 朝上，并添加色标。
axis image; set(gca,'YDir','normal'); colorbar;
% 设置 theta13 的标题和 x、y 轴标签。
title('\theta_{13}'); xlabel('x'); ylabel('y');

% 切换到第二个子图。
nexttile;
% 显示 theta23 的二维数值图。
imagesc(xPlot,yPlot,theta23);
% 在 theta23 上叠加同一组晶界线。
hold on; plotRasterBoundaries(grainGrid,xPlot,yPlot,dx,dy); hold off;
% 设置等比例、正确的 y 方向及色标。
axis image; set(gca,'YDir','normal'); colorbar;
% 设置 theta23 的标题和坐标轴标签。
title('\theta_{23}'); xlabel('x'); ylabel('y');

% 切换到第三个子图。
nexttile;
% 显示 theta33 的二维数值图。
imagesc(xPlot,yPlot,theta33);
% 叠加与数值网格一致的晶界线。
hold on; plotRasterBoundaries(grainGrid,xPlot,yPlot,dx,dy); hold off;
% 设置等比例、y 方向及色标。
axis image; set(gca,'YDir','normal'); colorbar;
% 设置 theta33 的标题与坐标轴标签。
title('\theta_{33}'); xlabel('x'); ylabel('y');

% 切换到第四个子图。
nexttile;
% 显示三个分量的范数 rhoTheta。
imagesc(xPlot,yPlot,rhoTheta);
% 在标量密度图上叠加晶界。
hold on; plotRasterBoundaries(grainGrid,xPlot,yPlot,dx,dy); hold off;
% 设置等比例、y 方向与色标。
axis image; set(gca,'YDir','normal'); colorbar;
% 设置标量密度图标题和坐标轴标签。
title('\rho_{\theta}^{2D}'); xlabel('x'); ylabel('y');

% 在整张主图上写入 dx、dy 信息。
sgtitle(sprintf('Lie-algebra disclination density | dx=%.4g, dy=%.4g',dx,dy));
% 把主图保存为 PNG 文件。
saveas(fig,figFile);

% 如果用户启用了取向诊断图，开始绘制该图。
if makeIPFDiagnostic
    % 尝试进行 MTEX 绘图，避免可选图失败时中断数值结果。
    try
        % 创建单独的白底取向诊断图窗。
        figIPF = figure('Color','w','Name','EBSD orientation diagnostic');
        % 让 MTEX 按 EBSD 取向绘制相图或取向着色图。
        plot(ebsd,ebsd.orientations);
        % 在取向图上叠加与密度图相同的像素晶界。
        hold on; plotRasterBoundaries(grainGrid,xPlot,yPlot,dx,dy); hold off;
        % 设置取向诊断图标题。
        title('Imported EBSD orientations and reconstructed grain boundaries');
        % 将诊断图保存在输入文件旁边。
        saveas(figIPF,fullfile(outDir,[baseName '_EBSD_grain_diagnostic.png']));
    % 捕获可选 MTEX 图生成失败的异常。
    catch ME
        % 只发出警告；此时前面的数值结果已保存。
        warning('Could not create optional MTEX IPF diagnostic: %s',ME.message);
    % 结束取向图绘制异常处理。
    end
% 结束可选取向图判断。
end

% 在命令窗口输出 CSV、MAT 与主图的保存路径。
fprintf('\nDone. Outputs saved to:\n  %s\n  %s\n  %s\n',csvFile,matFile,figFile);
% 提醒检查对称性取向代表元跳变和晶界伪影；提示文字跨两行。
fprintf(['Check the plots for symmetry-representation jumps and boundary artifacts. ' ...
    'Do not interpret the scalar norm alone as the Frank-vector direction.\n']);

%% 脚本局部函数
% 定义局部函数：输入 3×3 旋转矩阵 R，输出主值旋转向量 w。
function w = localRotmLog(R)
% 将 SO(3) 旋转矩阵取主值对数，输出单位为弧度的三分量旋转向量。
    % 用 trace(R)=1+2cos(theta) 求旋转角的余弦。
    c = (trace(R)-1)/2;
    % 把余弦限制在 [-1,1] 内，防止浮点误差使 acos 产生复数。
    c = max(-1,min(1,c));
    % 反余弦求得主值旋转角，取值为 0 到 pi。
    th = acos(c);
    % 取 R 的反对称部分对应的三维轴向向量。
    v = [R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
    % 对于接近零角度的旋转，使用矩阵对数的一阶近似。
    if th < 1e-8
        % 此时旋转向量约等于反对称部分的一半。
        w = 0.5*v;
    % 角度没有逼近 pi 时进入常规公式分支。
    elseif (pi-th) > 1e-5
        % 用 theta/(2sin theta) 把反对称部分换算为旋转向量。
        w = (th/(2*sin(th)))*v;
    % 对于接近 pi 的旋转，使用避免除以极小 sin(theta) 的分支。
    else
        % 接近 pi 时 sin(theta) 很小，改用对称部分稳定求旋转轴。
        % 构造接近轴向外积的对称矩阵。
        A = (R+eye(3))/2;
        % 找出最适合恢复轴向的最大对角元位置。
        [~,k] = max(diag(A));
        % 初始化三分量旋转轴。
        axisVec = zeros(3,1);
        % 用所选对角元的平方根恢复对应轴分量。
        axisVec(k) = sqrt(max(A(k,k),0));
        % 若该分量仍非常小，则采用特征分解备用方案。
        if axisVec(k) < 1e-8
            % 对 R 的对称部分求特征值和特征向量。
            [V,D] = eig((R+R')/2);
            % 找到特征值最接近 1 的旋转轴方向。
            [~,ii] = min(abs(diag(D)-1));
            % 把相应实特征向量当作旋转轴。
            axisVec = real(V(:,ii));
        % 轴分量足够大时按矩阵非对角元恢复其他分量。
        else
            % 循环求出三个轴分量。
            for j = 1:3
                % 只计算主轴分量以外的位置。
                if j ~= k
                    % 由对称矩阵关系恢复第 j 个轴分量。
                    axisVec(j) = (R(j,k)+R(k,j))/(4*axisVec(k));
                % 结束是否为其他分量的判断。
                end
            % 结束逐分量循环。
            end
        % 结束旋转轴备用恢复方法分支。
        end
        % 计算恢复出的轴向量长度。
        axisNorm = norm(axisVec);
        % 若长度接近机器精度，则不能稳定归一化。
        if axisNorm < eps
            % 使用单位 x 方向作为退化情形的备用轴。
            axisVec = [1;0;0];
        % 轴长度足够时进入正常分支。
        else
            % 将旋转轴归一化为单位向量。
            axisVec = axisVec/axisNorm;
        % 结束轴向量长度判断。
        end
        % 反对称部分有足够信息时，用它确定旋转轴的正负方向。
        % 若反对称部分可辨别符号，且轴与其方向相反，就纠正方向。
        if norm(v) > 1e-8 && dot(axisVec,v) < 0
            % 翻转旋转轴的符号。
            axisVec = -axisVec;
        % 结束轴方向判断。
        end
        % 用主值角乘单位轴得到旋转向量。
        w = th*axisVec;
    % 结束接近零角、一般角及接近 pi 的三种情况。
    end
% 结束旋转矩阵对数局部函数。
end

% 定义同晶粒差分函数：输入三分量场 F、晶粒号 G、有效点掩膜、方向 dim 和步长 h。
function D = maskedDerivative(F,G,valid,dim,h)
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
            if ~valid(i,j) || G(i,j)<=0, continue; end
            % 当 dim=2 时，沿列方向即 x 方向计算差分。
            if dim == 2
                % 检查左邻点是否存在、有效、属于同一晶粒且三个场分量都有限。
                hasM = j>1  && valid(i,j-1) && G(i,j-1)==G(i,j) && all(isfinite(reshape(F(i,j-1,:),1,[])));
                % 以同样条件检查右邻点。
                hasP = j<nx && valid(i,j+1) && G(i,j+1)==G(i,j) && all(isfinite(reshape(F(i,j+1,:),1,[])));
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
                hasM = i>1  && valid(i-1,j) && G(i-1,j)==G(i,j) && all(isfinite(reshape(F(i-1,j,:),1,[])));
                % 以同样条件检查下邻点。
                hasP = i<ny && valid(i+1,j) && G(i+1,j)==G(i,j) && all(isfinite(reshape(F(i+1,j,:),1,[])));
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
function plotRasterBoundaries(G,xCoords,yCoords,dx,dy)
% 只绘制两个相邻有效像素所属晶粒编号不同的交界。
% 晶界线段来自密度图所用的同一个矩形像素网格。
    % 标出左右相邻且都已索引、晶粒编号不同的像素对。
    diffX = G(:,1:end-1) ~= G(:,2:end) & G(:,1:end-1) > 0 & G(:,2:end) > 0;
    % 取得这些竖直晶界对应的行、列编号。
    [rowX,colX] = find(diffX);
    % 把竖直线段放在两个像素中心之间的 x 中点。
    xMid = xCoords(colX) + dx/2;
    % 取每条竖直线段所在的行中心 y 坐标。
    yMid = yCoords(rowX);
    % 生成线段两端 x 值，使用 NaN 隔开不同线段。
    xVert = [xMid(:), xMid(:), nan(numel(xMid),1)].';
    % 生成竖直线段上下端的 y 值，并用 NaN 分隔。
    yVert = [yMid(:)-dy/2, yMid(:)+dy/2, nan(numel(yMid),1)].';

    % 标出上下相邻且都已索引、晶粒编号不同的像素对。
    diffY = G(1:end-1,:) ~= G(2:end,:) & G(1:end-1,:) > 0 & G(2:end,:) > 0;
    % 取得这些水平晶界对应的行、列编号。
    [rowY,colY] = find(diffY);
    % 取得每条水平线段所在列的 x 中心。
    xMid = xCoords(colY);
    % 把水平线段放在两个像素中心之间的 y 中点。
    yMid = yCoords(rowY) + dy/2;
    % 生成水平线段左右端的 x 值及 NaN 分隔符。
    xHorz = [xMid(:)-dx/2, xMid(:)+dx/2, nan(numel(xMid),1)].';
    % 生成水平线段两端 y 值及 NaN 分隔符。
    yHorz = [yMid(:), yMid(:), nan(numel(yMid),1)].';

    % 一次画出竖直、水平晶界线段；下一行设置黑色和线宽。
    line([xVert(:); xHorz(:)],[yVert(:); yHorz(:)], ...
        'Color','k','LineWidth',0.35);
% 结束晶界绘图函数。
end

% 定义三维向量场叉乘函数，输入和输出都是 ny×nx×3 数组。
function C = crossField(A,B)
% 对两个 ny×nx×3 向量场在每个像素处逐点计算三维叉乘。
    % 预分配结果数组 C，大小与 A 相同。
    C = nan(size(A));
    % 拆出 A 在三个轴向的二维分量。
    ax=A(:,:,1); ay=A(:,:,2); az=A(:,:,3);
    % 拆出 B 在三个轴向的二维分量。
    bx=B(:,:,1); by=B(:,:,2); bz=B(:,:,3);
    % 按向量叉乘公式计算第一分量 ay*bz-az*by。
    C(:,:,1) = ay.*bz-az.*by;
    % 按向量叉乘公式计算第二分量 az*bx-ax*bz。
    C(:,:,2) = az.*bx-ax.*bz;
    % 按向量叉乘公式计算第三分量 ax*by-ay*bx。
    C(:,:,3) = ax.*by-ay.*bx;
% 结束逐像素叉乘函数。
end
