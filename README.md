# SCUDEM 2025 — Space Dust

**SCUDEM X 2025 · Outstanding Award**

A team modelling project on the atmospheric entry of space dust, coupling drag and ablation in a height–velocity–mass ODE. This release preserves the historical implementation alongside a guarded RK4 workflow with explicit termination states.

[Project case study / 项目案例](https://fury0831.github.io/projects/space-dust.html)

宇宙尘埃垂直进入大气时的减速与烧蚀模型，来自 SCUDEM 2025 Team 1075 Problem A 团队项目。代码耦合高度、向下速度与质量，使用手写 RK4 积分，并提供单例轨迹和小参数扫描。

团队成员：**Fengrui Zhang、Yiye Huang、Zening Wang**。Fengrui Zhang 负责建模、代码、可视化与最终视频录制。团队已同意按 [MIT License](LICENSE) 开源。

[项目演示视频](https://youtu.be/rn7Ku_a6Opg)。视频用于展示原团队项目，不代表当前整理版的每项数值修正都已在原演示中使用。

## 模型与数据

状态为 `[h, v, m]`，分别表示高度（m）、向下速度（m/s）和质量（kg）。模型假设球形恒密度粒子、垂直入射、恒定重力、静止指数大气，以及常数阻力和热传递系数。具体方程、假设及版本差异见 [MODEL_NOTES.md](MODEL_NOTES.md)。

**无需下载数据集。** 初始条件与物理常数确定数值实验，轨迹由求解器生成；这里没有观测数据拟合或实证验证。

## 安装与运行

推荐 Python 3.13.9，依赖 NumPy 2.3.5、Matplotlib 3.10.6，不需要 SciPy。在项目根目录执行：

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

# 单例轨迹。
python run.py demo --output results/demo

# 四组初始条件的小扫描。
python run.py scan --output results/scan
```

默认 `demo` 使用初始半径 **0.1 mm**、初始速度 **12 km/s**、初始高度 **100 km**，最长 **100 s**，基础时间步长 **0.01 s**。

默认 `scan` 使用半径 **0.05、0.1 mm** 与速度 **12、20 km/s** 的四个组合，同样从 100 km 开始，最长 100 s，基础步长 0.01 s。该小扫描不是原 notebook 的 250 点网格复现。

输出目录保存数值结果、摘要及图表。扫描图分别展示**终止时剩余质量**与**终止状态**，不将所有终点称作地面幸存。运行结果位于 `results/`，不作为源代码的一部分提交。

## 如何解释结果

请首先查看终止状态：

| 状态 | 含义 |
| --- | --- |
| `ground_reached` | 有效接受步跨过地面后，在该步内线性插值到 `h = 0` |
| `time_limit` | 达到最长积分时间；不表示粒子已经到达地面 |
| `mass_cutoff` | 质量达到声明的相对阈值；不表示物理上的完全烧蚀 |
| 数值或资源保护状态 | 阶段越界、非有限值或步数等限制阻止继续求解；不能解释为物理事件 |

摘要区分三个量：

- `end_mass_percent`：终止时的 `100 * m / m0`。
- `end_radius_percent`：终止时的 `100 * r / r0`；恒密度下等于 `100 * (m / m0)^(1/3)`。
- `ground_mass_percent`：仅在 `ground_reached` 时给出地面质量比例；其他状态为 `null`。

默认相对质量阈值为 **1e-12**。达到阈值只会标记 `mass_cutoff`，不会宣称“完全蒸发”。质量百分比和半径百分比也不是同一个量。

## 代码结构与整理改动

| 文件 | 用途 |
| --- | --- |
| `space_dust/historical.py` | 原核心函数及历史求解器、网格函数；保留历史行为，不自动运行 |
| `space_dust/model.py` | 半径/质量转换和有效域包装，保持原正域物理方程 |
| `space_dust/integrate.py` | 带阶段有效域检查、有限重试和明确终止状态的 RK4 |
| `space_dust/plotting.py` | 轨迹、终止质量及状态图 |
| `run.py` | `demo` / `scan` 命令行入口 |
| `tests/test_model.py` | 模型与数值保护的检查 |

本次模块化由 Codex 协助，包括命令行、阶段保护、事件与单位处理、文档和验证流程。原团队模型与后续数值整理明确区分；原 notebook 的内嵌图像、原题、模板和证书不纳入本仓库。

新的积分入口要求每个 RK4 阶段满足 `m > 0`、`v >= 0` 且状态有限，否则缩小为半步重试，并设最小步长、最大减半次数与最大步数以避免无限循环。地面跨越在有效接受步内插值处理。这些是**有效域保护，不是误差估计驱动的自适应精度控制**，不构成收敛或物理准确性保证。

## 验证状态

2026-10-25，在已有的 Python 3.13.9、NumPy 2.3.5、Matplotlib 3.10.6 环境中完成以下检查，没有安装新依赖：

```bash
python -m unittest discover -s tests -v
python run.py demo --t-max 2 --dt .001 --output results/smoke-demo
python run.py scan --t-max .2 --dt .001 --output results/smoke-scan
```

- 单元测试：**9 tests，OK**。
- 2 s 单例：退出码 0，终止状态为 `time_limit`。
- 0.2 s 四组合扫描：退出码 0，四个结果均为 `time_limit`。
- 这些结果的 `ground_mass_percent` 均为 `null`，没有被归类为地面到达。

这是当前模块化版本的**短程流程检查**，使用的时间上限和步长不同于默认示例。本轮没有运行默认 100 s 示例或历史 250 点网格，不能据此声称长期积分、全部参数空间或报告图已经复现。独立 DOP853 参考计算与这些短程命令属于不同验证范围，不能混为同一次运行。

原 notebook 保存的两个 0.1 mm 示例（12 和 20 km/s）均在 100 s 到达时间上限，文本明确写着未撞地、未完全烧毁。该历史记录不是本次整理版的运行结果。原 250 点大网格本次未重跑。

## 已知限制

- 这是简化模型，缺少入射角、分裂、辐射冷却、温度演化、变系数与真实大气结构；未进行观测校准。
- 不同初始半径在同一终止时间一般有不同高度与速度，不能声称剩余比例普遍与半径无关。
- 地面分布的逆推需要额外初始条件和观测信息；本项目没有实现逆求解器或证明所有反演都不可能。
- 可用原源码是质量状态实现，不包含报告所述独立半径状态求解器或纳米尺度图的生成代码。不能把当前运行称为这些报告图的完整复现。

Copyright (c) 2025 Fengrui Zhang, Yiye Huang, Zening Wang。使用与再分发请保留 [MIT 许可证](LICENSE) 要求的声明。

## Archived interactive figures / 原始交互图

[Open the figures on demand / 按需打开交互图](https://fury0831.github.io/scudem-2025-space-dust/)

| File | Metric | Original size |
| --- | --- | ---: |
| [meteor_simulation1.html](docs/meteor_simulation1.html) | Altitude / 高度 | 19.37 MB |
| [meteor_simulation2.html](docs/meteor_simulation2.html) | Radius remaining / 剩余半径 | 19.55 MB |
| [meteor_simulation3.html](docs/meteor_simulation3.html) | Velocity / 速度 | 20.59 MB |

These original files animate precomputed data, not a live ODE solver. Their exact original bytes are preserved; see [hash manifest](docs/ARCHIVE_MANIFEST.json). Opening a file starts its archived animation. The portfolio does not preload or embed the large files. Use a desktop browser and Reset axes if a previous frame leaves an unsuitable axis range.

Each figure contains 32 velocity frames (6–68 km/s, step 2) × 6 radii (0.001–100 μm). Of 192 trajectories, 54 record through 100 s, 125 stop at 0.01 s, and 13 end between 0.02 and 0.13 s. Stop reasons and the frame-generating solver were not supplied. A shortened curve or zero radius does not prove complete ablation; all recorded altitude endpoints remain above ground. The radius plot is not a mass-fraction plot.

The HTML velocity data are nonnegative, unlike one older report figure. Version equivalence has not been established. During local review all three plots rendered, parameter controls worked, and one `Uncaught (in promise)` console error was observed during velocity-slider/reset interactions; the plot remained usable. The historical files are preserved rather than silently repaired or regenerated.

三份文件保留原始字节和数据，属于预计算交互图，不是实时模拟。历史动画加载后会自动播放，建议桌面浏览。部分曲线非常短且无终止原因标签，不能称为地面存活结果；HTML 速度与报告旧图有差异，尚不能证明同一版本。第三张图在滑块/重置操作期间出现过一条 Promise 错误，仍可显示；没有为消除错误而更改历史文件。

Plotly.js and bundled components keep their original licenses. See [third-party notices](docs/licenses/THIRD_PARTY_NOTICES.md). No original report, presentation template, certificate or signatures are included.
