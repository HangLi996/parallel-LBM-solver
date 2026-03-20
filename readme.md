# 并行 LBM 流体求解器（Parallel LBM Solver）


## 项目简介（Overview）
- 基于 D2Q9 Lattice Boltzmann Method（LBM）的二维流体数值求解器，面向分布式 CPU 集群中的流场模拟与并行性能研究。
- 我主要完成了求解器的 MPI + OpenMP 分层并行化与 SIMD 向量化优化，并建立了 benchmark/profiling 流程用于性能分析。

---

## 项目特点（Features）
- 核心方法：LBM（D2Q9 速度格）。
- 并行策略：MPI（进程级）+ OpenMP（线程级）。
- 优化技术：编译器自动 SIMD、数据布局优化（SoA 思路）、分阶段计时分析。
- 性能分析：提供 benchmark 与 LIKWID profiling 脚本及结果数据。

---

## 方法与实现（Methods）
- 数值方法：采用 LBM 的 streaming / collision 框架，在离散速度格上推进分布函数并恢复宏观量。
- 算法结构：单步流程为 `stream -> communicate -> post-stream -> collide`，并记录初始化与计算阶段耗时。
- 并行设计：
  - **Domain decomposition**：支持 vertical、horizontal、2D block 三种划分；
  - **Halo exchange**：通过 `Messenger` 抽象与 MPI 非阻塞通信进行边界数据交换；
  - **Hybrid 并行**：计算密集循环用 OpenMP，通信阶段保持 MPI 主导。

---

## 性能分析（Performance）
- benchmark 内容：
  - `profiling/benchmark.sh`：不同网格规模下的运行时间测试；
  - `profiling/benchmark_results_simdVSnosimd.csv`：SIMD 与非 SIMD 对比（100 到 800 规模）。
- 关键结论：
  - 随问题规模增大，SIMD 版本优势更明显；
  - 在 `800x800` 时，SIMD 约 `23.67s`，非 SIMD 约 `143.86s`，约 `6x` 加速；
  - OpenMP 能明显降低计算核耗时，但高并行度下通信占比上升，整体趋于通信瓶颈。
- 可选图表：`profiling/timing_plot.png`、`profiling/speedup_plot.png`。

---

## 项目结构（Project Structure）
- `src/`：核心源码目录。
- `src/main.cpp`：参数解析、MPI 初始化、主流程与阶段计时。
- `src/LBM/`：LBM 核心计算与通信逻辑（含 `Simulation`、`Messenger`）。
- `src/Domains/`：域建模与划分、邻接关系初始化。
- `src/BoundaryConditions/`：边界条件实现。
- `src/VelocitySets/`：速度格定义（如 `d2q9`）。
- `src/Reporting/`：结果输出与统计。
- `profiling/`：基准测试脚本、性能分析脚本与结果文件。

---

## 编译与运行（Build & Run）
- 编译命令：

```bash
make clean
make lbm
```

- 可选分解版本：

```bash
make lbm_horizontal
make lbm_2d
```

- 运行命令：

```bash
mpirun -np <np> ./bin/main <P> <dx> <dy> <iterations>
```

- 示例：

```bash
mpirun -np 4 ./bin/main 4 80 80 1000
```

---

## 参数说明（Usage）
- `np`：`mpirun` 启动的 MPI 进程数。
- `P`：程序内部使用的进程数（建议与 `np` 一致）。
- `dx`：x 方向网格大小（默认 `80`）。
- `dy`：y 方向网格大小（默认 `80`）。
- `iterations`：迭代步数（默认 `11`）。
- 可调参数：
  - 分解策略：通过 `make lbm / lbm_horizontal / lbm_2d` 选择；
  - 编译优化：通过 `Makefile` 中 `CFLAGS` 切换 SIMD 与非 SIMD 配置；
  - OpenMP 线程数：通过环境变量 `OMP_NUM_THREADS` 设置。

---

## 运行环境（Environment）
### 硬件
- CPU：`Intel Xeon Platinum 8480+`（`Sapphire Rapids`）。
- 集群：CoolMUC-4（`106` 节点 / 每节点 `112` 物理核 / 总计 `11872` cores / 每节点 `512 GB` 内存）。

### 软件
- 编译器：`gcc/15.2.0`
- 构建工具：`ninja/1.12.1`
- MPI：`intel-mpi/2021.12.0`
- 其他：`OpenMP`（由编译参数 `-fopenmp` 启用）、`LIKWID`（可选，用于 profiling）。

---

## 报告（Report）
- 技术文档：`Project_report.pdf`
- 仓库链接（报告中引用）：[https://github.com/HangLi996/parallel-LBM-solver](https://github.com/HangLi996/parallel-LBM-solver)

---

## 注意事项（Notes）
- `2D block` 分解通常要求进程数可映射为规则二维进程网格（实践中建议平方数，如 4/9/16）。
- 网格规模建议与进程划分匹配，避免局部子域过小导致通信占比过高。
- 若要做 SIMD 对比，请确认 `Makefile` 中 `CFLAGS` 使用一致配置后再重新编译。
- 若运行失败，请先检查 `mpicxx/mpirun` 是否可用，以及模块环境是否正确加载。

---

## 参考（References）
- Succi, S. *The Lattice Boltzmann Equation for Fluid Dynamics and Beyond*. Oxford University Press, 2001.
- Kruger, T. et al. *The Lattice Boltzmann Method: Principles and Practice*. Springer, 2017.
- CoolMUC-4 文档：[https://doku.lrz.de/coolmuc-4-1082337877.html](https://doku.lrz.de/coolmuc-4-1082337877.html)
