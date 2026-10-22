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
%   3. 在矩阵对数和差分前检查索引质量，并在各晶粒的四邻接连通域内
%      仅选择晶体对称等价代表元；不平均取向、不合并母晶/孪晶。
%      无法连续提升的回路边及主值对数分支切口会被标记并禁止差分。
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
clearvars; clc;
% 将本脚本旁的数值辅助函数加入当前会话，不修改持久 MATLAB 路径。
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

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

% 索引质量阈值示例：MAD 单位通常为度，CI 无量纲，须按仪器导出定义核实。
% 字段匹配不区分大小写；不存在的字段会警告，不会伪装成质量已合格。
% 不给 BC/IQ 设置通用阈值，因其强烈依赖仪器、相和采集条件。
qualityRules = struct('field',{'mad','ci'},'minimum',{[],0.1},'maximum',{1.0,[]});
% 同晶粒的相邻点若约化取向差过大，也不传播或跨边差分，默认不超过重建阈值。
maxContinuityAngleDeg = 3; % 仅报告大物理增量，不作为删边阈值
% 离散图优化参数；容差针对浮点误差，不是物理取向差阈值。
representativeMaxIterations = 50; % 有限状态下降迭代上限，未收敛会标记
representativeEnergyTolerance = 1e-12; % rad^2，避免舍入误差引起反复切换
representativeAmbiguityTolerance = 1e-10; % rad^2，近简并候选需人工核查
representativeInvarianceTolerance = 1e-7; % rad，acos 在零附近的舍入放大
% 孤立异常：至少三个同晶粒邻点的约化角中位数超过 2 度时标记为可疑。
% 仅供诊断；不能凭此判据区分错误索引和真实局部取向梯度。
localOutlierAngleDeg = 2;
minLocalNeighbors = 3;
excludeLocalOutliers = false; % 必须为 false；代表元模块禁止按角度剔除
% 如已有独立识别的晶界/孪晶界，可给出保护边掩膜 MAT 文件（不合并孪晶）。
% 文件包含 blockedX(ny,nx-1)、blockedY(ny-1,nx)、xPlot、yPlot，坐标需与本图一致。
protectedBoundaryFile = '';
% 主值旋转向量跨边差大于该值（弧度）时，标记对数切口而非强行展开。
maxLogVectorJumpRad = pi;

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
% 明确列出输入相及未处理相；本脚本按目标相分别运行，不混用对称群。
inputPhaseIds=unique(ebsd.phaseId);
for pid=inputPhaseIds(:)'
    cs=ebsd.CSList{pid};
    fprintf('输入晶相 ID=%d, 名称=%s, 点数=%d\n',pid,cs.mineral,nnz(ebsd.phaseId==pid));
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
selectedPhaseIds=unique(ebsd.phaseId);
if numel(selectedPhaseIds)~=1
    error('请设置 phaseName 选择单一目标相，然后对其他相分别运行，禁止混用对称群。');
end
excludedPhaseIds=setdiff(inputPhaseIds,selectedPhaseIds);
if ~isempty(excludedPhaseIds)
    warning('本次仅计算 %s；未处理相 ID=%s，请修改 phaseName 后分别运行。', ...
        phaseName,mat2str(excludedPhaseIds(:)'));
end
% 在晶粒重建之前验证坐标，避免用非有限位置建立邻接关系。
if any(~isfinite(ebsd.x(:))) || any(~isfinite(ebsd.y(:)))
    error('EBSD coordinates contain NaN or Inf; inspect the imported scan.');
end
assert(isscalar(grainBoundaryAngleDeg) && isfinite(grainBoundaryAngleDeg) ...
    && grainBoundaryAngleDeg>0 && grainBoundaryAngleDeg<180, 'Invalid grain threshold.');
assert(isscalar(maxContinuityAngleDeg) && isfinite(maxContinuityAngleDeg) ...
    && maxContinuityAngleDeg>0 && maxContinuityAngleDeg<180, 'Invalid continuity threshold.');
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
    % 本地 MTEX 6.1.1 源码：第二输出是 grainId 数组，写入 EBSD 属性。
    [grains, ebsd.grainId] = calcGrains(ebsd,'angle',grainBoundaryAngleDeg*degree,'minPixel',1);
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

%% 索引质量检查与晶粒内对称等价取向连续化（先于对数和差分）
[qualityGood,seedScore,qualityReport] = checkIndexQuality(ebsd.prop,N,qualityRules);
Rall = double(matrix(ebsd.orientations));
if size(Rall,1)~=3 || size(Rall,2)~=3 || size(Rall,3)~=N
    error('Unexpected orientation matrix dimensions; expected 3x3xN.');
end
% 只使用晶体正旋转群。禁止加入试样对称操作或把孪晶当作对称操作。
ori = ebsd.orientations;
properCS = ori.CS.properGroup;
symOps = double(matrix(properCS.rot));
% 运行时再用当前加载的 MTEX 对象核实复合侧，避免依赖文档版本猜测。
for k=1:size(symOps,3)
    equivalentMatrix=double(matrix(ori(1)*properCS.rot(k)));
    if norm(equivalentMatrix-Rall(:,:,1)*symOps(:,:,k),'fro')>1e-10
        error('当前 MTEX 的矩阵复合与 R*S 不一致，停止以免使用错误坐标约定。');
    end
end
runtimeInfo=struct('ebsdClass',class(ebsd),'orientationClass',class(ori), ...
    'phaseId',selectedPhaseIds,'mineral',ori.CS.mineral, ...
    'properOperationCount',size(symOps,3),'orientationImplementation',which('orientation'));
fprintf('EBSD 类型=%s，取向类型=%s，晶相=%s，正旋转群阶数=%d\n', ...
    runtimeInfo.ebsdClass,runtimeInfo.orientationClass,runtimeInfo.mineral, ...
    runtimeInfo.properOperationCount);

options.maxEdgeAngle = maxContinuityAngleDeg*degree;
options.maxIterations=representativeMaxIterations;
options.energyTolerance=representativeEnergyTolerance;
options.ambiguityTolerance=representativeAmbiguityTolerance;
options.invarianceTolerance=representativeInvarianceTolerance;
options.phaseId=ebsd.phaseId;
options.localOutlierAngle = localOutlierAngleDeg*degree;
options.minLocalNeighbors = minLocalNeighbors;
options.excludeLocalOutliers = excludeLocalOutliers;
options.blockedX = false(ny,nx-1); options.blockedY = false(ny-1,nx);
xPlot = x0+(0:nx-1)*dx; yPlot = y0+(0:ny-1)*dy;
if ~isempty(protectedBoundaryFile)
    boundary = load(protectedBoundaryFile,'blockedX','blockedY','xPlot','yPlot');
    required = {'blockedX','blockedY','xPlot','yPlot'};
    if ~all(isfield(boundary,required))
        error('Protected-boundary file must contain blockedX, blockedY, xPlot and yPlot.');
    end
    if numel(boundary.xPlot)~=nx || numel(boundary.yPlot)~=ny ...
            || any(abs(boundary.xPlot(:)-xPlot(:))>coordTol*20) ...
            || any(abs(boundary.yPlot(:)-yPlot(:))>coordTol*20) ...
            || any(~isfinite([boundary.xPlot(:);boundary.yPlot(:)]))
        error('Protected-boundary coordinates do not match this EBSD grid.');
    end
    if ~all(ismember(boundary.blockedX(:),[0,1])) ...
            || ~all(ismember(boundary.blockedY(:),[0,1]))
        error('Protected-boundary masks must contain only 0 or 1.');
    end
    options.blockedX = logical(boundary.blockedX);
    options.blockedY = logical(boundary.blockedY);
end
field = prepareSymmetryContinuousField(Rall,grainId,linIdx,[ny,nx], ...
    symOps,qualityGood,seedScore,options);
grainGrid = field.grainGrid; validGrid = field.validGrid;
% 把原始扫描中的相信息映射到已经验证的网格，仅作诊断，不插值、不改坐标。
field.phaseGrid=zeros(ny,nx);
ax=double(ebsdAll.x(:)); ay=double(ebsdAll.y(:));
ai=round((ax-x0)/dx)+1; aj=round((ay-y0)/dy)+1;
onGrid=isfinite(ax) & isfinite(ay) & ai>=1 & ai<=nx & aj>=1 & aj<=ny ...
    & abs(ax-(x0+(ai-1)*dx))<=max(coordTol*20,dx*1e-3) ...
    & abs(ay-(y0+(aj-1)*dy))<=max(coordTol*20,dy*1e-3);
originalPhaseIds=double(ebsdAll.phaseId(:));
field.phaseGrid(sub2ind([ny,nx],aj(onGrid),ai(onGrid)))=originalPhaseIds(onGrid);
field.phaseBoundaryX=field.phaseGrid(:,1:end-1)~=field.phaseGrid(:,2:end);
field.phaseBoundaryY=field.phaseGrid(1:end-1,:)~=field.phaseGrid(2:end,:);
field.grainBoundaryX=grainGrid(:,1:end-1)~=grainGrid(:,2:end);
field.grainBoundaryY=grainGrid(1:end-1,:)~=grainGrid(2:end,:);

if nnz(validGrid)<9
    error('Fewer than 9 usable points remain after quality/rotation checks.');
end
% 保留原始代表元和修正代表元；两者用完全相同的质量掩膜与导数边。
omegaRaw = nan(ny,nx,3); omegaGrid = nan(ny,nx,3);
for p = find(field.validPoint)'
    omegaRaw(iy(p),ix(p),:) = reshape(localRotmLog(field.Rraw(:,:,p)),1,1,3);
    omegaGrid(iy(p),ix(p),:) = reshape(localRotmLog(field.Rcontinuous(:,:,p)),1,1,3);
end
% SO(3) 的连续代表元仍可能跨主值 log 的 pi 切口，不能靠平滑掩盖。
assert(isscalar(maxLogVectorJumpRad) && isfinite(maxLogVectorJumpRad) ...
    && maxLogVectorJumpRad>0 && maxLogVectorJumpRad<=pi, 'Invalid log-cut threshold.');
logJumpX = sqrt(sum(diff(omegaGrid,1,2).^2,3));
logJumpY = sqrt(sum(diff(omegaGrid,1,1).^2,3));
logCutX = field.edgeX & logJumpX>maxLogVectorJumpRad;
logCutY = field.edgeY & logJumpY>maxLogVectorJumpRad;
edgeX = field.edgeX & ~logCutX; edgeY = field.edgeY & ~logCutY;
qualityRejectedGrid = false(ny,nx); qualityRejectedGrid(linIdx)=~qualityGood;
suspectGrid = false(ny,nx); suspectGrid(linIdx)=field.suspectLocal;
fprintf('Quality rejected: %d; invalid SO(3): %d; local suspects: %d (excluded=%d).\n', ...
    nnz(~qualityGood),nnz(field.invalidRotation),nnz(field.suspectLocal),excludeLocalOutliers);
fprintf('Independent domains: %d; loop-conflict edges: %d; log-cut edges: %d.\n', ...
    max(field.componentGrid(:)),nnz(field.loopConflictX)+nnz(field.loopConflictY), ...
    nnz(logCutX)+nnz(logCutY));
if any(field.loopConflictX(:)) || any(field.loopConflictY(:)) ...
        || any(logCutX(:)) || any(logCutY(:))
    warning('Some edges remain discontinuous; both derivative stages exclude them. Inspect diagnostics.');
end

%% 用同一公式、同一质量掩膜和同一边掩膜计算处理前后结果
% 仅移除对称代表元跳变；不把晶界/孪晶界两侧取向拉到一起。
[thetaRaw,~,~] = densityFromOmega(omegaRaw,grainGrid,validGrid,edgeX,edgeY,dx,dy);
[theta,kappaX,kappaY] = densityFromOmega(omegaGrid,grainGrid,validGrid,edgeX,edgeY,dx,dy);

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
% 分相运行不得覆盖另一相的结果。
baseName = [baseName '_phase' num2str(selectedPhaseIds)];
% 设定 CSV 表格的保存路径和文件名。
csvFile = fullfile(outDir,[baseName '_LieDisclinationDensity.csv']);
% 设定 MAT 数据文件的保存路径和文件名。
matFile = fullfile(outDir,[baseName '_LieDisclinationDensity.mat']);
% 设定四宫格 PNG 图的保存路径和文件名。
figFile = fullfile(outDir,[baseName '_LieDisclinationDensity.png']);

% 按网格原点及步长生成每个像素的二维 x、y 坐标。
[Xgrid,Ygrid] = meshgrid(min(x)+(0:nx-1)*dx, min(y)+(0:ny-1)*dy);
% 建立结果表；以下三行依次放入坐标、晶粒号、掩膜、三个分量及范数，并命名表列。
indexedGrid = false(ny,nx); indexedGrid(linIdx)=true;
invalidRotationGrid = false(ny,nx); invalidRotationGrid(linIdx)=field.invalidRotation;
raw13=thetaRaw(:,:,1); raw23=thetaRaw(:,:,2); raw33=thetaRaw(:,:,3);
T = table(Xgrid(:),Ygrid(:),grainGrid(:),indexedGrid(:), ...
    theta13(:),theta23(:),theta33(:),rhoTheta(:),validGrid(:), ...
    qualityRejectedGrid(:),invalidRotationGrid(:),suspectGrid(:), ...
    raw13(:),raw23(:),raw33(:), ...
    'VariableNames',{'x','y','grainId','indexed','theta13','theta23','theta33','rhoTheta', ...
    'usable','qualityRejected','invalidRotation','localSuspect','rawTheta13','rawTheta23','rawTheta33'});
% 把结果表写入 CSV 文件。
writetable(T,csvFile);
% 把主要结果和中间变量写入 MAT 文件；这条 save 语句跨三行。
save(matFile,'theta13','theta23','theta33','rhoTheta','omegaGrid', ...
    'kappaX','kappaY','grainGrid','validGrid','validResult','dx','dy', ...
    'grainBoundaryAngleDeg','inputFile','xVals','yVals','omegaRaw','thetaRaw', ...
    'field','qualityRules','qualityReport','qualityGood','qualityRejectedGrid','suspectGrid', ...
    'edgeX','edgeY','logCutX','logCutY','options','maxLogVectorJumpRad', ...
    'protectedBoundaryFile','symOps','indexedGrid','invalidRotationGrid', ...
    'inputPhaseIds','selectedPhaseIds','excludedPhaseIds','runtimeInfo');
writeContinuityDiagnostics(field,validResult,grainGrid,xPlot,yPlot,outDir,baseName,degree);

%% 保存处理前后分量比较及质量诊断（同分量共用色标）
compareFile = fullfile(outDir,[baseName '_LieComponents_before_after.png']);
compareComponents(thetaRaw,theta,grainGrid,xPlot,yPlot,dx,dy, ...
    {'\theta_{13}','\theta_{23}','\theta_{33}'},compareFile);
compareComponents(omegaRaw,omegaGrid,grainGrid,xPlot,yPlot,dx,dy, ...
    {'\omega_1 (rad)','\omega_2 (rad)','\omega_3 (rad)'}, ...
    fullfile(outDir,[baseName '_RotationVector_before_after.png']));
% 诊断图保留晶界；红/黄点仅为质量/局部异常标记，不是孪晶判据。
figQuality = figure('Color','w','Name','质量及连续性诊断');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile; imagesc(xPlot,yPlot,double(validGrid)); axis image; set(gca,'YDir','normal');
hold on; plotRasterBoundaries(grainGrid,xPlot,yPlot,dx,dy);
scatter(x(~qualityGood),y(~qualityGood),12,'r','filled');
scatter(x(field.suspectLocal),y(field.suspectLocal),18,'y'); hold off;
title('可用点；红：质量剔除，黄：局部可疑'); xlabel('x'); ylabel('y'); colorbar;
nexttile; imagesc(xPlot,yPlot,field.componentGrid); axis image; set(gca,'YDir','normal');
hold on; plotRasterBoundaries(grainGrid,xPlot,yPlot,dx,dy);
plotBlockedEdges(grainGrid,edgeX,edgeY,xPlot,yPlot,dx,dy); hold off;
title('独立连续域；红线：额外禁用差分边'); xlabel('x'); ylabel('y'); colorbar;
saveas(figQuality,fullfile(outDir,[baseName '_Quality_continuity.png']));

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
fprintf('Before/after comparison: %s\n',compareFile);
fprintf(['Check the plots for symmetry-representation jumps and boundary artifacts. ' ...
    'Do not interpret the scalar norm alone as the Frank-vector direction.\n']);

%% 脚本局部函数
% 定义局部函数：输入 3×3 旋转矩阵 R，输出主值旋转向量 w。
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
function [good,score,report] = checkIndexQuality(prop,n,rules)
% 按原字段逐项检查：仅记录字段名与阈值，不猜测缺失指标或插补坏点。
    good = true(n,1); score = zeros(n,1);
    names = fieldnames(prop); report = struct('field',{},'available',{},'rejected',{});
    for k = 1:numel(rules)
        match = find(strcmpi(names,rules(k).field));
        report(k).field = rules(k).field;
        report(k).available = ~isempty(match); report(k).rejected = 0;
        if isempty(match)
            warning('Quality field "%s" absent; its threshold was NOT checked.',rules(k).field);
            continue;
        end
        if numel(match)~=1, error('Ambiguous quality field: %s',rules(k).field); end
        values = prop.(names{match});
        if ~isnumeric(values) || ~isreal(values) || numel(values)~=n
            error('Quality field %s must be a real numeric value per point.',names{match});
        end
        values = double(values(:)); ok = isfinite(values);
        low = rules(k).minimum; high = rules(k).maximum;
        if ~isempty(low)
            assert(isscalar(low) && isfinite(low),'Invalid quality minimum.');
            ok = ok & values>=low;
        end
        if ~isempty(high)
            assert(isscalar(high) && isfinite(high),'Invalid quality maximum.');
            ok = ok & values<=high;
        end
        if ~isempty(low) && ~isempty(high), assert(low<=high,'Invalid quality interval.'); end
        good = good & ok; report(k).rejected = nnz(~ok);
        fprintf('Quality %s: min=%s max=%s, rejected %d/%d.\n', ...
            names{match},mat2str(low),mat2str(high),nnz(~ok),n);
        % 种子排序只在有单向质量定义时使用，分位次序不改变任何取向。
        if xor(isempty(low),isempty(high))
            direction = 1;
            if ~isempty(high), direction = -1; end
            eligible = find(ok);
            if ~isempty(eligible)
                v = direction*values(eligible); [~,~,rank] = unique(v);
                score(eligible) = score(eligible)+double(rank)/max(rank);
            end
        end
    end
    if ~any([report.available])
        warning('No configured quality fields found; quality remains unverified. Inspect indexing diagnostics.');
    end
end

function compareComponents(before,after,g,x,y,dx,dy,labels,fileName)
% 每列同一分量共用对称色限；不裁剪极值或各自归一化来制造改善。
    fig = figure('Color','w','Name','处理前后分量比较');
    tiledlayout(2,3,'Padding','compact','TileSpacing','compact');
    for component = 1:3
        a = before(:,:,component); b = after(:,:,component);
        shared = isfinite(a) & isfinite(b);
        values = [a(shared);b(shared)];
        limit = max(abs(values));
        if isempty(limit) || limit==0, limit = 1; end
        % 比较图只显示共同有限的像素；MAT 中保留两者各自的完整结果。
        a(~shared)=NaN; b(~shared)=NaN;
        for stage = 1:2
            nexttile((stage-1)*3+component);
            imageData = a; stageName = '原始代表元';
            if stage==2, imageData=b; stageName='连续等价代表元'; end
            h = imagesc(x,y,imageData); set(h,'AlphaData',isfinite(imageData));
            axis image; set(gca,'YDir','normal'); caxis([-limit,limit]); colorbar;
            hold on; plotRasterBoundaries(g,x,y,dx,dy); hold off;
            title(sprintf('%s | %s',stageName,labels{component})); xlabel('x'); ylabel('y');
        end
    end
    sgtitle('相同质量/差分掩膜；同分量相同色标');
    saveas(fig,fileName);
end

function plotBlockedEdges(g,edgeX,edgeY,x,y,dx,dy)
% 只标记同晶粒内禁用的边；跨晶粒边仍由黑色晶界线表示。
    blockedX = ~edgeX & g(:,1:end-1)>0 & g(:,1:end-1)==g(:,2:end);
    blockedY = ~edgeY & g(1:end-1,:)>0 & g(1:end-1,:)==g(2:end,:);
    [i,j] = find(blockedX); xm=x(j)+dx/2; ym=y(i);
    xv=[xm(:),xm(:),nan(numel(i),1)]';
    yv=[ym(:)-dy/2,ym(:)+dy/2,nan(numel(i),1)]';
    [i,j] = find(blockedY); xm=x(j); ym=y(i)+dy/2;
    xh=[xm(:)-dx/2,xm(:)+dx/2,nan(numel(i),1)]';
    yh=[ym(:),ym(:),nan(numel(i),1)]';
    line([xv(:);xh(:)],[yv(:);yh(:)],'Color','r','LineWidth',0.5);
end
