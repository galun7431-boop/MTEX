# MTEX Scripts

MATLAB scripts for EBSD analysis,
twinning and disclination density.

## 向错密度脚本的晶体对称连续化

运行 `mtex_disclination_density_Lie_3.m` 前，先在 MATLAB 中启动 MTEX 7.2.0，
并修改脚本顶部的输入路径、相名、晶粒阈值及质量阈值。辅助函数
`prepareSymmetryContinuousField.m` 必须与脚本放在同一目录。

处理顺序为：索引相选择与晶粒重建 → 质量检查 → 同晶粒四邻接图 →
选择 `R*S` 晶体正旋转等价代表元 → 主值矩阵对数 → 受边界掩膜约束的两阶差分。
连续化不做取向平均、不施加试样对称性、不合并母晶与孪晶。晶界两侧和
被质量空洞分隔的连通域独立处理，不从远处跳过坏点建立邻接。

MAD 上限 1 度、CI 下限 0.1 是可配置的示例，须核实导出字段定义与仪器条件。
缺失字段会警告并记录为未检查；BC/IQ 没有通用阈值。默认标记但不自动剔除
局部约化取向差异常点，以免误删真实梯度；启用剔除须检查敏感性。
拒绝非有限、反射和明显偏离 SO(3) 的矩阵，只修正微小舍入误差。

晶粒重建默认阈值为 3 度，不执行孪晶合并。阈值本身不是孪晶识别算法，
若重建漏掉已知晶界/孪晶界，应提供 `protectedBoundaryFile`：MAT 文件包含
逻辑数组 `blockedX(ny,nx-1)`、`blockedY(ny-1,nx)` 和网格中心坐标 `xPlot`、
`yPlot`。true 表示禁止跨越该边；脚本会核对尺寸与实际扫描坐标。
四邻接中超过连续性阈值的边、回路代表元不一致的边和处理后 log 分支切口
也会禁止传播或差分。无法连续化的区域不应靠平滑强行连接。
晶界处仅求同侧导数，输出不包含跨晶界跳跃的分布型界面密度。

新增输出（在 EBSD 文件旁）：

- `_LieComponents_before_after.png`：原始/连续代表元的三个 theta 分量。
- `_RotationVector_before_after.png`：三个旋转向量分量比较。
- `_Quality_continuity.png`：质量剔除、局部可疑点、连续域及禁用边。
- MAT 文件额外保存原始/修正场、对称操作编号、质量报告、回路/分支诊断和参数。

前后结果使用相同质量点掩膜、相同差分边；同一分量共用色标，只显示双方
共同有限的像素，不用不同归一化制造“改善”。原有李代数截断公式未更改；
绝对主值 log 及其截断公式并非对任意代表元选择都不变，图像变化不能单独
证明模型正确。需用已知场、原论文公式及真实数据进一步验证。

静态检查（不执行 MATLAB）：

```bash
/workspace/tools/matlab-static/bin/mh_lint --ignore-config --input-encoding utf-8 \
  scripts/mtex_disclination_density_Lie_3.m scripts/prepareSymmetryContinuousField.m \
  tests/test_symmetry_continuity.m
```

以后获准运行 MATLAB 时，可从项目根目录执行：

```matlab
results = runtests('tests/test_symmetry_continuity.m');
assertSuccess(results);
```

这组八项合成回归测试覆盖代表元跳变、真实梯度、晶粒/孪晶量级边界、
错误晶粒标签、质量空洞、显式保护边界、无效旋转和局部可疑点。
本次修改未执行 `.m`，未生成真实 EBSD 前后图，也未验证数值结果。
