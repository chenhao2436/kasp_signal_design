# KASP Signal Design Workflow

[English](README_EN.md) | **简体中文**

**KASP Signal Design Workflow** 用于从 BSA 分析结果中筛选适合后续 KASP 引物设计的候选 SNP。流程将 SNP-index 信号排序、VCF 基因型过滤、候选位点整理和参考序列提取连接起来，并通过 Snakemake 管理多区间任务。

> **项目状态：研究工作流 / Alpha。** 当前仓库包含可运行的 Snakemake 规则和脚本，但配置文件仍含有特定服务器路径；独立服务器部署前需要替换路径并完成实际数据验证。仓库中的旧版 Bash 总控脚本依赖未纳入仓库的 INI 配置文件，暂不建议作为默认入口。

## 功能

- 从指定 BSA 峰区中提取左右两侧候选位点，并按 `DELTA_SNP_INDEX` 排序。
- 支持 flank 和 bin 两种候选区域选择模式。
- 从 VCF 检查候选位点的亲本基因型、池基因型及邻近变异。
- 过滤非 SNP、多等位位点、亲本非纯合或亲本基因型相同的位点。
- 整理每侧推荐位点和 `extract_fa.py` 所需输入。
- 从参考基因组提取 SNP 上下游序列，并生成带有 `[REF/ALT]` 标记的序列表。
- 使用 Snakemake 对多个目标区间执行统一流程。

## 工作流程

```text
SNP-index / DELTA_SNP_INDEX table
              |
              v
    01_rank_snpindex_candidates
              |
              v
       ranked candidate SNPs
              |
              +---- Joint VCF + sample genotypes
              |                 |
              v                 v
       02_filter_isolated_parental_snps
              |
              v
         PASS / FAIL sites
              |
              v
     03_prepare_extract_fa_input
              |
              v
       extract_fa.py + reference FASTA
              |
              v
     KASP flanking sequence table
```

## 仓库结构

```text
kasp_signal_design/
├── Snakefile
├── environment.yml             # Snakemake 驱动环境
├── data/
│   ├── *_config.yaml           # 分析参数配置示例/工作配置
│   ├── *_regions.tsv           # 待分析区间
│   └── kasp_rule_env.yaml      # Snakemake rule 运行环境
├── scripts/
│   ├── 20260418_kasp_01_rank_snpindex_candidates.py
│   ├── 20260418_kasp_02_filter_isolated_parental_snps.py
│   ├── 20260418_kasp_03_prepare_extract_fa_input.py
│   ├── extract_fa.py
│   └── 20260418_kasp_00_run_pipeline.sh  # 旧版 Bash 总控脚本
├── doc/                        # 分步骤说明和历史运行记录
├── results/                    # 分析结果
└── temp/                       # 日志及临时文件
```

## 软件要求

### Snakemake 驱动环境

- Linux 环境（推荐在 Linux 服务器运行）。
- Conda / Mamba。
- Snakemake。

使用仓库中的环境文件创建驱动环境：

```bash
conda env create -f environment.yml
conda activate kasp-snakemake
snakemake --version
```

### Snakemake rule 环境

规则环境定义在 `data/kasp_rule_env.yaml`，包括：

- Python 3.12
- bcftools
- htslib
- pyfaidx

运行时使用 `--use-conda`，Snakemake 会根据规则环境文件创建独立环境。当前 Python 候选排序、过滤和整理脚本主要使用 Python 标准库；VCF 基因型检查通过外部 `bcftools` 调用，序列提取使用 `pyfaidx`。

不需要为本仓库单独维护 pip 版 `requirements.txt`：当前依赖已经按 Snakemake 驱动环境和 rule 环境拆分。

## 输入数据

每个分析区间需要在配置文件中提供以下信息：

| 配置键 | 说明 |
|---|---|
| `snpindex_path` | BSA SNP-index / plot 结果表 |
| `vcf_path` | bgzip 压缩并建立索引的联合 VCF |
| `reference_fasta` | 与 VCF 对应的参考基因组 FASTA |
| `parent1_sample` / `parent2_sample` | 两个亲本在 VCF 中的样本名 |
| `high_pool_sample` / `low_pool_sample` | 高池和低池在 VCF 中的样本名 |
| `flank_no_variant_bp` | 目标 SNP 邻域过滤距离，当前配置通常为 50 bp |
| `target_markers_per_side` | 每侧目标推荐标记数量 |
| `initial_rank_limit` | 初始候选检查规模 |
| `expansion_step_bp` | 候选不足时向峰区中心扩展的步长 |
| `regions_file` | 区间 ID 与染色体区间对应表（可选，默认路径见 Snakefile） |

SNP-index 输入表应包含 `CHROM`、`POS`、`DELTA_SNP_INDEX` 等列；无表头输入时，脚本按默认列位置解析。VCF 必须包含配置中填写的样本名，并有可用的 `.tbi` 或 `.csi` 索引。参考 FASTA 的染色体命名应与 SNP-index 和 VCF 保持一致。

**重要：** 仓库中的配置文件保留了原始分析环境中的绝对路径和样本名。运行前必须检查并修改这些值；不要直接把示例路径当作通用路径。

## 快速开始

### 1. 检查配置

当前 `Snakefile` 顶部默认加载：

```text
data/70_Chr11_255900000_260500000_config.yaml
```

请先打开该文件，修改 `snpindex_path`、`vcf_path`、`reference_fasta` 和样本名等参数。区间列表由配置中的 `regions_file` 指定。

如需切换另一套区间配置，请同步修改 `Snakefile` 的 `configfile` 指向目标 YAML 文件，并确认相应 regions TSV 存在。

### 2. 先做 dry-run

```bash
snakemake -n
```

检查 Snakemake 展开的规则和目标输出是否符合预期。dry-run 不会执行分析命令。

### 3. 正式运行

```bash
snakemake --use-conda -j 1
```

建议第一次先使用单任务运行，确认输入路径、索引和样本名正确后再增加并行度。

### 4. 输出位置

每个 region 在 `results/<region_id>/` 下生成：

| 文件 | 内容 |
|---|---|
| `01_ranked_candidates.tsv` | 左右侧候选位点排序表 |
| `01_rank_summary.tsv` | 候选数、窗口范围和扩展层级汇总 |
| `02_checked_sites.tsv` | 所有已检查位点及过滤原因 |
| `02_pass_sites.tsv` | 通过硬过滤的候选位点 |
| `03_extract_fa_input.tsv` | 给序列提取脚本的四列输入，无表头 |
| `03_recommended_markers.tsv` | 推荐标记及候选不足原因 |
| `04_kasp_sequences.tsv` | SNP 上下游序列及 REF/ALT 标记 |

运行日志位于 `temp/<region_id>/`。

## 过滤逻辑说明

当前候选过滤脚本主要检查：

1. 目标位点是否能在 VCF 中找到。
2. 位点是否为双等位 SNP，而非 InDel 或多等位位点。
3. 两个亲本是否分别为纯合基因型。
4. 两个亲本的基因型是否不同。
5. 高池和低池的样本基因型是否可读取。
6. 位点上下游指定 flank 范围内是否存在其他变异。

默认配置通常以 50 bp 作为邻近变异检查范围，每侧目标推荐 2 个标记。具体值以实际 YAML 配置为准。

**候选通过不代表 KASP 引物已经通过实验验证。** 输出序列是供后续引物设计和实验筛选使用的候选材料。

## 已知限制与注意事项

- 当前配置含有服务器绝对路径，需要用户按自己的数据目录修改。
- 当前 Snakefile 默认加载特定的 Chr11 配置，不是通用的所有区间默认配置。
- 旧版 `scripts/20260418_kasp_00_run_pipeline.sh` 默认引用仓库中未提供的 `20260418_kasp_pipeline.ini`；本 README 推荐使用 Snakemake 入口。
- 当前没有提供可公开分发的最小测试 VCF、SNP-index 数据和参考基因组，因此尚不能宣称跨服务器端到端验证通过。
- 参考序列提取脚本会尝试不同染色体名称格式，但正式分析前仍需确认参考序列与 VCF 的组装版本一致。
- KASP 引物的特异性、扩增效率和分型聚类仍需通过专门引物设计工具及实验验证。

## 验证建议

每次分析至少检查：

- `snakemake -n` 是否显示预期的 region 和输出文件。
- VCF 是否可读取、是否包含配置的四个样本、索引是否存在。
- `02_checked_sites.tsv` 的失败原因是否合理。
- `03_recommended_markers.tsv` 是否达到每侧目标数量。
- 输出序列中参考碱基是否与参考基因组一致。

## 路线图

- [ ] 将服务器绝对路径替换为通用配置模板。
- [ ] 为 Snakemake 和旧版 Bash 入口统一参数定义。
- [ ] 增加小型模拟数据和自动化测试。
- [ ] 增加配置校验命令和更清晰的运行报告。
- [ ] 增加 KASP 引物设计结果的可选下游检查。

## 引用与问题反馈

本仓库用于研究分析流程。引用时请注明仓库地址、使用的 Git commit，以及实际使用的软件和版本。

问题反馈请提供配置文件（去除私人路径和敏感信息）、命令、日志片段和输入文件格式说明。
