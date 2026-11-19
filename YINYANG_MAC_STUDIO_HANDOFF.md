# AMRVAC 阴阳网格：全新 Mac Studio 开发交接

更新日期：2026-11-19。对象：没有旧电脑工作目录、编译缓存或聊天记录的新 Mac Studio，以及在该机器上接手的开发者／AI。

**先读结论：第一版固定网格框架进入收尾，不再把消除所有拼接误差作为交付条件。保留中点电场交换，以 HLL＋MC 为已测试的推荐配置；新的面平均 ghost 重构默认关闭。下一阶段优先讨论一个范围明确的太阳日冕／太阳风应用，而不是继续无止境优化局部重构。**

本文是交接时的状态快照。旧文档中的“CT 尚未启用”“下一步必须研究重构”等文字有历史上下文；先按本文区分当前状态，再查原始报告。本文没有宣布已经完成太阳风模型或已经确定下一阶段的全部物理方案。

## 1. 接手位置与代码版本

| 项目 | 内容 |
|---|---|
| 用户个人仓库 | <https://github.com/nanami0721/amrvac_ionization> |
| SSH 地址 | `git@github.com:nanami0721/amrvac_ionization.git` |
| 当前开发分支 | `codex/yinyang-grid-prep` |
| 已验证数值代码基线 | `d84eb55417f373bfbdce55acd6ede44a98735120` |
| 基线提交标题 | `Complete fixed-grid Yin-Yang CT integration and validation` |
| 本文提交 | 在上述基线之后，属于交接文档更新；后续应记录实际 checkout 的 SHA |
| 原始上游 | `git@github.com:amrvac/amrvac.git`；不是本次 push 的目标 |

旧机器工作区是 `/Users/nanami/workspace/2608_uawsom/yinyang-recovery-20261116/source`。**新机器不必创建这个路径，也不要把旧机器的其他 `/Users/nanami/codes/amrvac` checkout 当作当前开发版本。** 使用下面的新路径即可。

只下载这一个 Markdown 文件，不等于下载了源码或测试数据。先用有权限的 GitHub 身份 clone 正确仓库。新机器自行配置 SSH；不要把旧机器的私钥提交进仓库。

```sh
mkdir -p "$HOME/codes"
git clone --origin personal --branch codex/yinyang-grid-prep \
  git@github.com:nanami0721/amrvac_ionization.git "$HOME/codes/amrvac-yinyang"
cd "$HOME/codes/amrvac-yinyang"
git status --short
git branch --show-current
git remote -v
git rev-parse HEAD
git merge-base --is-ancestor d84eb55417f373bfbdce55acd6ede44a98735120 HEAD
```

最后一个命令应成功退出；否则先确认分支，不要在错误版本上开发。本文之后若分支继续演进，应比较后续提交，而不是盲目回退到旧 SHA。下载 ZIP 的用户也建议改用 clone，以便保留版本和差异证据。

## 2. 已完成、实验性和未完成部分

### 2.1 已实现且有相应验证的部分

| 模块 | 当前能力与证据边界 |
|---|---|
| 阴阳坐标与双 patch | 独立 `yinyang` 坐标标识；各 patch 复用球坐标度量。Yang 使用正交旋转 `(x,y,z) -> (-x,z,y)`。固定三维网格、全局时间步。 |
| 单元变量交换 | 原始量交换、笛卡尔矢量变换、角向 ghost 与角点填充、预计算供体／MPI 通信表。不是严格守恒 overset 耦合。 |
| GLM 路径 | 理想总能量 MHD；已有均匀场、输运、散度诊断和 GLM 重启／MPI 比较。 |
| 中心引力 | 已实现并完成 GLM 下静水平衡／潜势偶极场验证，包括十个声穿越时间和径向加密。**不能直接外推成 CT＋太阳风已验证。** |
| CT 演化集成 | 显式启用的固定网格 UCT-HLL、SSPRK3 staged 演化；边电场生成、跨 patch 交换、离散旋度更新已接通。 |
| 面／边同步 | 几何 owner、同 patch 重复边同步、跨 patch 边交换、面 ghost 填充、混合角点完成、径向 ghost 散度闭合。 |
| 变量一致性 | 面磁场与中心磁场同步；ghost 区磁能变化补偿以保持热压。不是新的低 beta 总能量算法。 |
| 电场积分约定 | donor 和 receiver 统一为原生 `dt × 边长 × 边中点电场`。真实输运测试显示收益，保留此实现。 |
| 回归与资源控制 | 仓库外构建、完整六配置库编译、MPI 对照、debug／优化对照、重复调用 stderr 检查、进程树与系统内存监视。六配置编译不代表支持一／二维阴阳网格。 |

### 2.2 已实现但默认关闭的实验选项

`yy_ct_face_moments=T`：用 24 个供体面的面积平均磁场拟合局部线性、无散的笛卡尔磁场，QR 与几何系数在通信表初始化时计算。

- 修正了把面平均磁场当作点值处理的问题：常量、线性无散场的 ghost 误差降至约 `1e-15`。
- 修正了正权拟合在网格节点附近因舍入而选取不同邻居的可复现性问题。
- **没有改善所测非线性磁环的完整输运；正反向 L2 略增约 0.009%／0.026%。**
- 独立拟合的各面不自动满足共同的离散散度约束；后续闭合仍可能放大误差。
- 因此默认 `F`，不能把它宣传成精度升级。不要在迁移时顺手开启。

### 2.3 尚未实现、尚未验证或明确延期

| 项目 | 状态／不要混淆的地方 |
|---|---|
| 完整 Yin–Yang CT 重启验收 | **待补。** GLM 重启与普通球坐标 staged CT 重启已有证据，不能替代真实阴阳 CT 的双 patch 面场／中心场连续运行对照。 |
| 太阳风物理边界 | 待设计和实现。当前实验 CT 要求径向所有变量为 `cont`，不是线 tied／光球驱动／特征太阳风流出边界。不能只改参数绕过拒绝检查。 |
| 1–20 个太阳半径应用 | 是拟讨论方向，尚无完成的开放日冕／太阳风模型与验收。已有小球壳静平衡不等于这个应用。 |
| CT＋引力的完整应用资格 | 组合仍需针对性验证；引力的历史长期平衡证据属于 GLM。 |
| 日冕热传导、辐射损失、加热／波模型 | 阴阳集成尚未完成；一般 AMRVAC 有某模块，不代表阴阳路径已支持。 |
| 低 beta 精度改进 | 研究暂缓。已有小网格 GLM 测试揭示压力／虚假流速误差，普通球坐标也有相近问题；没有解决这些误差，也没有建立 CT 低 beta 应用资格。 |
| 背景场分离、Boris 等 | 不属于当前已交付阴阳能力；现有限制仍生效。 |
| 联合面磁通重构＋离散散度约束 | 只形成了后续研究建议，没有实现；不是 V1 交付前置项。 |
| 严格全局守恒 overset／AMR／自动负载迁移／局部时间步 | 当前不支持或明确延期。patch 内 CT 散度控制不等于跨 patch 全局磁通守恒。 |
| 大规模扩展性与生产太阳爆发模拟 | 未完成资格验证；不能由单机小算例推断。 |

## 3. 固定采用的数值配置与限制

下列是配置片段，**不是可以独立运行的完整参数文件**。完整网格、初值、输出和边界定义从仓库测试中复制。

```fortran
&methodlist
 time_stepper='threestep'
 ct_stage_split=T
 yy_ct_experimental=T
 yy_ct_face_moments=F
 flux_scheme=20*'hll'
 limiter=20*'woodward'
/
&mhd_list
 typedivbfix='ct'
 type_ct='uct_hll'
 boundary_divbfix=6*.false.
/
&meshlist
 refine_max_level=1
 max_blocks=128
/
```

`woodward` 是当前配置中的 MC limiter。`usr_init` 先调用 `set_coordinate_system('yinyang')`，再激活物理模块。当前只支持三维、正半径球壳、匹配的径向网格和均匀角向网格；MPI 世界进程数为不小于 2 的偶数。patch 内 rank／网格编号与世界 rank 不同。禁止把两者混用。

源代码中的显式拒绝检查是能力边界，不应为跑通新模型直接删除。CT 当前只支持经典、不分裂的总能量 MHD；径向边界必须是 `cont`。开启 `yy_ct_face_moments` 而不开启实验 CT 会明确报错。

## 4. 已得到的关键结果：知道数值代表什么

细磁环输运：每 patch `16×32×96`，弱磁场峰值 `0.001`，磁环半径 `0.65`，球壳 `0.75–2.25`，固定 `dt=0.002`，终止时间 `0.6`。数值是终止时刻、体积加权的相对磁场 L2，不是所有日冕算例的误差上限。

| 配置 | 正向 L2 | 反向 L2 |
|---|---:|---:|
| 旧电场交换，TVDLF＋van Leer | 0.236885335 | 0.240640492 |
| 中点一致交换，TVDLF＋van Leer | 0.231467245 | 0.234869429 |
| 旧电场交换，HLL＋MC | 0.158538900 | 0.166044455 |
| **中点一致交换，HLL＋MC，当前推荐** | **0.149844810** | **0.155309843** |
| 再开启实验面平均 ghost 重构 | 0.149857689 | 0.155350753 |

中点一致交换的 HLL＋MC 正反向 L2 改善约 5.5%／6.5%。所测原生 active 散度保持舍入误差水平；中点版本的两／四 rank 单元和面输出逐位一致。实验面重构的 MPI 证据来自解析／组件测试，**没有声称它的细网格四 rank 完整输运已经跑过**。

此前“约 1.51 阶”的争论已经做过普通球坐标与分辨率对照，不能把它当作已证明的算法阶数缺陷重新启动无限排查。阅读 [球坐标收敛记录](doc/yinyang_spherical_convergence.md) 和 [恢复／内存记录](doc/yinyang_convergence_memory_recovery.md)。所有阶数都须说明具体场型、加密方向、分辨率区间和范数。

最近两个验证集合分别记录了 77 和 74 个正常算例；它们属于不同开发步骤，不应相加当作互不重叠的独立证据。详细计数、负测试和范围见 [进度与证据](doc/yinyang_fixed_grid_progress.md)。

## 5. 新机器环境：先建立本机工具链

旧机实测为 Apple Silicon `arm64`、GNU Fortran 14.2.0、Open MPI 5.0.5、Python 3.12。旧机 MPI 位于 `/usr/local/openmpi/bin`，Python 位于用户的 Conda 环境；**这些目录不必、也不应照搬到新机。** 新机可使用原生 Homebrew。版本变化要记录；不要求跨编译器逐位复现旧二进制，但同一工具链下既定 MPI／debug 回归不得自行放宽容差。

在正常 Terminal 中安装 Xcode Command Line Tools 和 Homebrew（若尚未安装），再执行：

```sh
xcode-select -p
# 若上一条失败，先运行 xcode-select --install，并完成系统安装界面。
# Homebrew 安装入口：https://brew.sh/
eval "$(/opt/homebrew/bin/brew shellenv)"
brew install gcc open-mpi python@3.12
"$(brew --prefix python@3.12)/bin/python3.12" -m venv "$HOME/.venvs/amrvac-yinyang"
source "$HOME/.venvs/amrvac-yinyang/bin/activate"
python -m pip install numpy matplotlib
uname -m
command -v gfortran mpif90 mpirun make perl rsync
mpif90 --showme:command
mpirun --version
gfortran --version
python --version
```

应为原生 arm64；MPI wrapper 应调用可用的 GNU Fortran。不要混合 Rosetta/x86_64 MPI 与 arm64 库。基础编译需要 GNU Fortran、MPI、make、Perl；NumPy 用于分析，Matplotlib 只用于可选绘图。若已有不同版本工具链，先记录、验证，不要无依据覆盖系统安装。

部分旧脚本会在 PATH 前加 `/usr/local/openmpi/bin`；新机这个目录不存在时会继续查找正常 PATH。若该目录存在，确认其工具链与当前选择一致。另有脚本直接执行绝对路径，见第 8 节；仅修改 PATH 不能修复这种情况。

## 6. 资源与工作区规则：接手后继续有效

1. 源码 checkout 只放源码、正式测试定义、文档和 Git 开发。**所有构建、模拟、分析、图片和日志在仓库外。** `AMRVAC_DIR` 必须指向外部源码副本，不能指回 checkout，也不能用可写符号链接绕过。
2. 系统实际占用不超过用户要求的总内存 70%；现有 watchdog 更保守：**65% 拒绝启动、68% 停止**。计算口径为 anonymous＋wired＋physical compressor，另记文件缓存与 swap；不要把文件缓存再加进去误判占用。
3. 编译进程树上限 **4096 MiB**；模拟／分析 **1024 MiB**；默认每 0.2 秒采样。这不是内核硬限额，必须预留余量。新 Studio 内存更多也不自动提高任务上限。
4. 编译使用 `make -j1`；先两 rank、`OMP_NUM_THREADS=1`、`OPENBLAS_NUM_THREADS=1`。编译并行度、MPI rank 数、OpenMP 线程数是三件事。不要按 CPU 核心数直接放大作业。
5. 测试显式使用 `max_blocks=128`，先小算例测量，再预测加密内存。若已有 runner 默认参数没有这个设置，先在外部测试定义中补齐；不要直接采用大块池默认值。
6. watchdog 返回非零、达到内存／时间阈值、预测超过安全门槛，都不算测试通过。检查 `summary.json` 的 `stop_reason`，保留日志；不能删除保护或提高阈值以掩盖失败。
7. macOS guard 需要 `/bin/ps`、`sysctl`、`vm_stat` 和 `libproc`。AI 沙箱若拒绝这些查询，应在允许的执行环境运行／申请相应执行权限，不能关闭监视。不要把权限错误当成系统内存不足。
8. 不重新执行此前接近 1 TiB 的大网格方案。详细日志只保存文件，控制台输出摘要和有界片段。正常 debug 模拟 stderr 应为空；不得重定向到 `/dev/null` 隐藏临时数组警告。
9. 共享 MHD 宏修改须按项目约定覆盖全部六配置，含 1D2V，并编译完整模块依赖。遵守数组生命周期、真实循环依赖和运行时诊断约束。提交前确认没有把构建产物带入仓库。

旧工作区的用户级 `AGENTS.md` 约束不一定随 clone 自动出现；本节保留了交接所需的核心规则。新机开发者仍须读取其实际环境中的有效 AGENTS.md，不能仅靠过去聊天记忆。

## 7. 不依赖旧机缓存的首次启动路线

下面是基于当前脚本接口核对的迁移步骤，**不是已在新 Mac Studio 上执行通过的报告**。不需要旧机的 Conda 路径、二进制或 `/private/tmp` 数据。按顺序执行；任一步失败就定位该步，不继续加大规模。

### 7.1 定义路径并验证监视器

激活上节虚拟环境，保持这些变量在当前 shell 中：

```sh
export YY_REPO="$HOME/codes/amrvac-yinyang"
export YY_RUN_ROOT="$HOME/amrvac-runs/yinyang-migration-$(date +%Y%m%d-%H%M%S)"
export YY_PY="$(command -v python)"
export YY_GUARD="$YY_REPO/tests/yinyang/interface_transport/convergence/memory_guard.py"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$YY_RUN_ROOT"
cd "$YY_REPO"
"$YY_PY" "$YY_GUARD" --out "$YY_RUN_ROOT/guard-probe" --limit-mib 1024 -- \
  "$YY_PY" -c 'import time; time.sleep(1)'
```

检查 `guard-probe/summary.json`：`returncode=0`、`stop_reason=null`，总内存应与 Studio 实际硬件一致。目录已存在时不要覆盖，改用新的任务根目录。

### 7.2 冷构建一个小型原生球坐标 CT 检查，生成本机缓存

这一步验证工具链和 CT 基础，**还不是阴阳 CT 验收**。仅从干净的已提交 checkout 导出；不要用 `git archive` 假装包含未提交修改。

```sh
set -e
set -o pipefail
export YY_BOOT="$YY_RUN_ROOT/01-native-ct"
test -z "$(git -C "$YY_REPO" status --porcelain)"
mkdir -p "$YY_BOOT/source" "$YY_BOOT/case"
git -C "$YY_REPO" rev-parse HEAD > "$YY_BOOT/source-commit.txt"
git -C "$YY_REPO" archive HEAD | tar -x -C "$YY_BOOT/source"
cp "$YY_BOOT/source/tests/yinyang/ct_preflight/mod_usr.t" "$YY_BOOT/case/"
cp "$YY_BOOT/source/tests/yinyang/ct_preflight/amrvac.par" "$YY_BOOT/case/"
"$YY_PY" - <<'PY'
import os
from pathlib import Path
p = Path(os.environ['YY_BOOT']) / 'case/amrvac.par'
s = p.read_text()
s = s.replace('&meshlist', '&meshlist\n max_blocks=128', 1)
s = s.replace('&methodlist', '&methodlist\n ct_stage_split=T', 1)
p.write_text(s)
PY
export AMRVAC_DIR="$YY_BOOT/source"
(
  cd "$YY_BOOT/case"
  "$YY_PY" "$YY_GUARD" --out "$YY_BOOT/setup-memory" --limit-mib 1024 -- \
    "$AMRVAC_DIR/setup.pl" -d=3 -arch=debug
  "$YY_PY" "$YY_GUARD" --out "$YY_BOOT/build-memory" --limit-mib 4096 -- \
    make 'F90=mpif90 --param ggc-min-expand=5 --param ggc-min-heapsize=4096' -j1
  "$YY_PY" "$YY_GUARD" --out "$YY_BOOT/run-memory" --limit-mib 1024 -- \
    mpirun --oversubscribe -np 2 ./amrvac
)
test ! -s "$YY_BOOT/run-memory/command.err"
```

GNU `--param` 设置降低编译器缓存开销，不是放宽浮点精度。冷构建可能较慢；guard 默认每项超时 3600 秒。超时后先看日志，不无限等待，也不因此抬高内存上限。

### 7.3 真正的阴阳 CT 小规模输运

```sh
cd "$YY_REPO"
export YY_PILOT="$YY_RUN_ROOT/02-yy-ct-pilot"
"$YY_PY" tests/yinyang/ct_transport/run.py "$YY_PILOT" \
  --cache-source "$YY_BOOT/source" --pilot \
  --radius 0.65 --radial-min 0.75 --radial-max 2.25
"$YY_PY" "$YY_GUARD" --out "$YY_RUN_ROOT/02-analysis-memory" --limit-mib 1024 -- \
  "$YY_PY" tests/yinyang/ct_transport/analyze.py "$YY_PILOT"
```

`--cache-source` 是刚在本机生成的完整外部源码目录，不是裸 `lib/`，也不是旧机绝对路径。runner 会复制源码、记录 commit／diff／哈希、校验缓存差异，并把 `AMRVAC_DIR` 指向自己的副本。

该 pilot 是已有 TVDLF＋van Leer 的两 rank 正向粗网格测试，目的在于迁移冒烟，**不应期待得到第 4 节 HLL＋MC 细网格的 0.149844810**。检查 `status.json`、`results.json`、`debug-d1-p1-s1-np2-a1/transport-history.txt` 和相应 `*-memory/`。应完成 300 步到 `t=0.6`、正常结束、无非法状态、运行 stderr 为空。

### 7.4 一次完整的边界／编译回归，然后停止扩大迁移验收

```sh
"$YY_PY" tests/yinyang/ct_boundary_sync/run.py "$YY_RUN_ROOT/03-boundary" \
  --cache-source "$YY_PILOT/source"
```

这一个 runner 已包含 debug／优化、uniform／stretched 径向网格、不同块布局、2／4／6 ranks、预期拒绝测试和六种完整库编译。负测试的非零退出是预期结果，必须由 runner 核对对应错误；正常模拟的 stderr 仍须为空。**新环境一旦通过这一组迁移门槛，不要再无条件重跑所有历史精度研究。**

随后记录 Studio 的硬件、编译器／MPI／Python 版本、源码 SHA、参数、结果摘要和内存峰值。跨机器差异先定位工具链和配置；没有理由时不要改算法或放宽断言。

## 8. 哪些测试不能从 clone 后直接照抄运行

| 入口／材料 | 新机处理方式 |
|---|---|
| `ct_preflight/build_stage_suite.py` | 支持冷构建，但其 `stage_suite.py` 写死 `/usr/local/openmpi/bin/mpirun`。Homebrew 新机先做路径可配置化并验证，或先用第 7 节路线。不要为了匹配旧路径创建指回旧机器／错误架构的链接。 |
| `run_validation.py` | 可指定 `--mpiexec`，用于 GLM；默认测试参数未显式写 `max_blocks=128`，迁移前应先补齐外部算例的资源设置。不要把它的 GLM 重启结果标为 CT。 |
| `ct_static/run.py`、`ct_accuracy/run.py`、`ct_ghost_fill/run.py` 等 | 需要已有外部 `source/lib` 缓存；第 7 节可生成。不要在所有 runner 之间循环等待“别人先提供缓存”。 |
| `ct_transport/run_controls.py` | 除源码缓存外，还读取旧 wide-transport 结果中的指定内存记录。任意新建 source 目录并不够；需重建对应前置 campaign，或针对新机增加独立的 preflight 入口。 |
| `ct_transport/run_representation.py` | 成对历史比较依赖 legacy/candidate 源码快照、reference-results 中的参数和内存记录。仓库里有 runner，不代表这些外部输入也已上传。 |
| `analyze_representation.py` 等历史比较 | 需匹配相应结果目录和标签；不能拿另一套 run 的 JSON 冒充。 |
| 旧 `/Users/nanami/...`、`/private/tmp/amrvac-tests/...` 引用 | 是证据位置，不是新机必须存在的环境。若没有迁移原数据，就从定义重建或明确标为旧机报告。 |

新机若要适配脚本路径／冷启动接口，这是允许的工程收尾；保持数值参数、阈值和验收语义，记录修改，不趁机扩展数值算法。

## 9. 旧机数据、文献和环境：哪些没有随 Git 迁移

Git 已包含数值源码、测试定义、分析脚本和正式结果报告。**没有包含大规模 `.dat`、可执行文件、编译库、完整 stdout/stderr、外部 source 快照、Zotero PDF 或参考代码库。** 仅靠 push 本文不能声称这些都已搬到 Studio。

如果需要保留原始证据，在旧机器清理前另做外部备份；新机开发不以取得这些备份为启动前提：

| 旧机材料 | 原位置／用途 |
|---|---|
| 最近面重构审计 | `/private/tmp/amrvac-tests/ct-face-20261119-baseline/final-audit.json` 与 `final-definitions/` |
| 面重构完整结果 | 同前缀的 `final2`、`ghost-final`、`transport`、`switch`、`boundary`、`patch` |
| 中点电场审计与输运 | `/private/tmp/amrvac-tests/ct-midpoint-20261119-baseline/`、同前缀 `transport`、`analytic-final` 等 |
| 原生网格对照 | `/private/tmp/amrvac-tests/ct-controls-20261119-suite/` |
| 引力长期平衡证据 | `/Users/nanami/workspace/2608_uawsom/yinyang-gravity-20261117/` |
| Zotero 文献 | 用户的本地 Zotero 库／附件；通过 Zotero 同步或独立备份迁移，不依赖旧条目本地路径 |
| 外部参考实现 | `/Users/nanami/workspace/2608_MAP/`，并未作为本仓库一部分上传 |

旧临时目录可能因清理／重启丢失。不要把路径当成持久归档保证。备份应保留 `amrvac.par`、`mod_usr.t`、source commit/diff/manifest、分析结果及需要复核的原始输出；记录目录清单、大小和 SHA-256 校验。不要把旧 `lib/` 当作新机可直接复用的工具链。备份后再校验，不要无差别 `git clean -fdx` 或删除唯一结果。

文献起点：Kageyama & Sato 2004（10.1029/2004GC000734）；Luo & Fan 2025 Yin–Yang MFE（10.3847/1538-4365/adf8d1）；Wongwathanarat 等 2010（10.1051/0004-6361/200913435）；Peng 等 2006（10.1256/qj.05.18）。用户已经在旧 Zotero 下载 Wongwathanarat 和 Peng，但新机应确认附件确已迁移。早期 [调研记录](doc/yinyang_preparation.md) 中的“尚未下载／尚未读完”是当时状态，不是当前附件清单。不要声称 MFE 的全部物理模块已在 AMRVAC 复现。

## 10. 下一阶段计划：先定范围，按关口结束

### 已与用户确定的边界

V1 终点是可复现、稳定、误差和限制明确的固定网格版本。**不追求拼接误差消失，不要求所有算例严格二阶，不把 ghost 研究继续滚动加入必做清单。** 只有实际稳定性／散度／MPI 缺陷，或一个明确科学任务在事先约定容限下失败且证据指向拼接，才重新打开算法研究。

### 建议给 Studio 的顺序

| 阶段 | 工作 | 交付／停止条件 |
|---|---|---|
| M0：迁移验收 | 按第 7 节建立本机环境、完成小型 CT 与边界回归，记录工具链和资源。 | 能独立冷启动、正常运行、关键回归通过。完成后停止扩展历史测试。 |
| M1：与用户确定应用契约 | 讨论简化偶极场日冕／太阳风；确认是否采用 1–20 Rsun、能量／加热近似、重力和边界物理、研究关心的输出。 | 写一页范围与验收表：具体配置、终止时间、对照、内存预算、量化容限、最多的测试轮次。**这些选择尚未获逐项确认。** |
| M2：必要可靠性缺口 | 在上长时间应用前补真实 Yin–Yang CT 重启对照：连续与分段续跑、两 patch 面场和中心场、元数据及同 MPI 规模。 | 对照达标或明确定位修复；不顺带开发可变 MPI 重启／AMR。 |
| M3：最小物理实现 | 按 M1 契约补 CT 物理径向边界和引力组合验证。若选择太阳风，可先用廉价一维／普通球坐标 Parker 类控制，再进入三维偶极场。 | 初值、单位、边界／入射特征数量、能量假设明确；不靠删除 startup 限制跑通。 |
| M4：单个三维示范 | 使用已冻结的中点 EMF、HLL＋MC、默认 ghost 方案。按约定做两级分辨率、一个有限时长稳定性窗口和需要的原生对照。 | 验收通过则交付参数、脚本、结果、资源与限制；失败则只处理已定位瓶颈，不无限加算例。 |

M1 建议明确测量：径向质量通量、速度／密度／温度剖面、开放磁通或所关心的拓扑量、原生散度、正定性、边界与拼接分区误差，以及实际耗时／内存。物理模型、目标量不同，不能用一个统一的“15% 磁场误差”或“2%”替代全部验收。

**默认不纳入下一阶段：**低 beta 方法专项、联合守恒重构、AMR、Boris／背景分离、完整波加热与辐射传导、生产级爆发、跨机器性能竞赛。只有用户明确选定应用确实需要时，再单独排期。

## 11. 源码与文档导航

| 位置 | 接手时看什么 |
|---|---|
| `src/mod_yinyang.t` | mesh／plan 适配、实际 runtime CT 接线、变量与边界同步、配置限制、实验开关 |
| `src/mod_yinyang_staggered.t` | 面／边供体表、矢量重构、实验面平均矩、QR、节点 tie 规则 |
| `src/mod_yinyang_ct_coupling.t` | 边积分近似、面交换、ghost 散度闭合 |
| `src/mod_yinyang_patch_sync.t` | 同 patch 面／边 owner 与重复量同步 |
| `src/mod_advance.t`、`src/mod_ghostcells_update.t` | RK-stage 和 ghost 生命周期中的实际调用位置 |
| `src/mhd/mod_mhd_phys.t`、`src/physics/mod_physics.t` | CT API、中心磁场／能量同步及支持模式 |
| `src/io/mod_input_output.t` | 参数解析、默认关闭实验选项、输出／重启相关接口 |
| [CT 接口说明](doc/yinyang_ct_stage_api.md) | 当前集成、默认值、物理边界限制与实验功能 |
| [开发进度和证据](doc/yinyang_fixed_grid_progress.md) | 按开发阶段记录的测试、正／负结果、范围冻结；历史叙述不能覆盖后续结论 |
| [引力验证](doc/yinyang_gravity_validation.md)、[低 beta 验证](doc/yinyang_low_beta_validation.md) | GLM 证据及不能外推到 CT／开放日冕的限制 |
| `tests/yinyang/ct_static/`、`ct_transport/` | 实际 CT 演化、输运、原生对照与分析 |
| `tests/yinyang/ct_accuracy/`、`ct_ghost_fill/`、`ct_boundary_sync/`、`ct_patch_sync/` | 算子诊断、面平均身份检查、边界／变量同步、MPI 一致性 |

## 12. 给新机器 AI 的第一条任务（可直接复制）

> 请先阅读仓库根目录的 `YINYANG_MAC_STUDIO_HANDOFF.md`，核对当前分支和 commit。先只执行 M0 移机验收：使用本机原生工具链，全部构建和测试放到仓库外，维持系统 65% 启动／68% 停止、编译 4 GiB、模拟和分析 1 GiB 的保护。不依赖旧机器的绝对路径或缓存，不开启 `yy_ct_face_moments` 作为应用默认，不重新启动拼接精度专项。记录已验证、待验证及失败原因；M0 完成后，向我汇报并讨论 M1 的物理范围与验收标准，再开始应用开发。未经我要求不 push。
