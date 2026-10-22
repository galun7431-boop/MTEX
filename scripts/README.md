# MTEX 向错密度脚本

实际计算入口为 `scripts/mtex_disclination_density_Lie_3.m`。
根目录同名脚本是兼容入口，统一执行该文件。项目中不存在
`mtex_disclination_density_Lie.m`。参数在实际入口顶部设置。

## 代表元选择与真实取向差

仅按选定晶相处理，保留默认 `phaseName='Magnesium'` 和 3° 晶粒重构阈值。
输入所有已索引相的 ID、名称和数量会打印，未处理相会警告并保存。
如需 Al2Y3，请把 `phaseName` 改成导入数据中的准确名称后单独运行。
输出文件带 `_phase<ID>` 后缀，避免分相运行互相覆盖。
多相不允许直接传入单群代表元模块；其他相在空间模板中是空洞。

只接受原坐标验证通过的轴对齐矩形网格；六角/畸变/散点报错，不重采样。
`linIdx` 根据实测 x、y 构造，不能使用导入顺序代替空间邻接。
同晶粒、同相、质量合格且没有显式保护边的四邻接图，按连通域独立处理：

1. 高质量种子开始传播，候选始终是该测点 `Rraw*S`。
2. 固定各连通域种子，通过离散坐标下降最小化所有图边的旋转角平方和。
3. 保存操作编号、原始/选定矩阵、能量历史、迭代次数与诊断标记。
4. 在原图上重新计算对称约化错配角，误差超过 `1e-7 rad` 即停止。

不做 SO(3) 投影、取向平均、平滑、插值或局部异常点剔除。
原始 EBSD 对象的取向不被修改。`maxContinuityAngleDeg=3` 和
`localOutlierAngleDeg=2` 仅为诊断，不能说明取向有错，也不会据此删边/删点。
真实局部 5° 甚至更大的晶内旋转不因角度大小自动屏蔽。
`excludeLocalOutliers` 必须为 false；独立的索引质量模块仍使用原有
MAD/CI 示例规则，缺失字段会警告，仪器阈值需自行核实。

50 次迭代是资源上限，不能证明全局最优。下降容差 `1e-12 rad²` 防止舍入引起
反复切换；候选能量差不超过 `1e-10 rad²` 作为近简并待检查标记。
未收敛连通域、近简并点、代表元残余增量大于物理约化角的回路边和原有
主值 log 切口会限制差分；这些不是“物理错配过大”的判据。
每个晶粒的规范独立，不计算跨晶界/相界的取向梯度，也不进行孪晶分类。
晶界错配分析应另外使用 MTEX 的对称约化方法，不用选定矩阵直接相减代替。

## API 与坐标约定

项目说明标注 MTEX 7.2.0，但此云端实际源码
`/workspace/dependencies/mtex/VERSION` 为 **mtex-6.1.1**。
没有 MATLAB/Octave，未启动这些 MTEX 类。已检查的真实源码：

- `geometry/@orientation/mtimes.m`：取向将晶体向量映射到试样向量。
- `geometry/@quaternion/symmetrise.m`：晶体操作 `q*CS.rot` 在右，试样操作在左。
- `geometry/@symmetry/properGroup.m`、`geometry/@rotation/matrix.m`：提取晶体正旋转。
- `EBSDAnalysis/@EBSD/calcGrains.m`：第二返回值为 grainId 数组，修复原脚本
  `[grains,ebsd]` 为 `[grains,ebsd.grainId]`。
- `geometry/@orientation/byMatrix.m` 和 `angle.m`：验证程序使用的官方接口。

没有施加试样对称性；矩阵实现不涉及四元数 q/-q 的符号问题。
在用户实际 MTEX 版本中还须运行下述对象测试。

## 差分与结果

保留原主值 log 和李代数截断公式；`localRotmLog`、`maskedDerivative`、
`crossField`、`densityFromOmega` 从脚本原样提取为独立函数，以便测试生产代码。
两个求导阶段都使用同一受保护边掩膜和有限值检查；中心模板不完整时用
同侧一阶差分，无同侧邻点保留 NaN。二阶计算不能通过中间 kappa 跳过边界。

保留已有向错密度、晶界叠加、前后对比及质量图；新增：

- `_Symmetry_diagnostics.png`：约化角/旋转增量分布、最终有效像素、
  屏蔽像素分类、原始相 ID、晶界与相界位置。
- `_Grain_diagnostics.csv`：各晶粒有效/屏蔽像素数和回路异常边数。
- MAT 中 `field` 保存约化角前后值、原始/选定增量、操作编号、
  近简并/未收敛/大物理增量诊断、相网格及边界掩膜。

约化角应保持不变；代表元增量变小本身不能证明向错密度正确。
绝对主值 log 加截断公式不对任意整体代表元选择保持严格不变，且有 pi 分支问题；
本次没有替换该定义。定量解释仍需原论文、参考场及真实 EBSD 验证。
未选相/未索引/无测点的密度均不计算，不作界面分布型密度解释。

## 验证与运行

已执行：11 个 MATLAB 文件的 MISS_HIT 静态检查（包括测试文件），
以及 `git diff --check`，均通过。**未执行任何 MATLAB 数值测试或真实 EBSD 运行**，
原因是环境没有 MATLAB、Octave 或许可证；不能把静态检查当作数值通过。

在 MATLAB 中从项目根目录执行：

```matlab
addpath('tests');
run_symmetry_validation; % 无需 MTEX；直接调用生产数值函数
results = runtests('tests/test_symmetry_continuity.m');
assertSuccess(results);
% 先启动用户安装的 MTEX，再执行：
validate_mtex_symmetry_api;
```

`run_symmetry_validation` 包含：A 均匀场随机等价表示；B 平缓场及真实 5° 总变化；
局部真实 5° 跃迁不删边；C 两个具有真实 20° 错配的晶粒及边界同侧密度；
D 两种不同正旋转群按相分别处理。A/B 的梯度及密度容差为 `1e-10`
（步长为 1 的合成场），物理约化角容差为 `1e-7 rad`。
对象验证另用 MTEX 官方 angle 对照约化角并检验矩阵右乘复合。
所有这些数值/对象测试目前均为**已编写、待运行**。
