#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""nee_recon —— 论文《Lineage structure and molecular diversity of microbial
strategies across terrestrial Mars analogues》v6 图件复现 CLI。

设计原则
- 全部使用相对路径；仓库根由 --input-root 指定（默认=本脚本上一级目录）。
- 产物写入 --output-root（默认 <input-root>/figures_out）。
- **绘图不重新发明**：十张图由 `scripts/figures/manuscript_figures.py` 中
  从作者端原始绘图代码**逐段移植**的函数产出（见该文件头部来源索引）。
- 输入仅为本仓库 data/derived 与 data/upstream 下已打包的小体积衍生表 / 源表 / 原始矢量。
- 不声称 raw-to-paper 全流程复现：重型上游（拼接/注释/系统发育/贝叶斯/AI 模型）
  未随仓库打包，逐条点名见 reproducibility_report.md。

子命令
  list                    列出十张图及状态
  verify                  校验 manifest 中列出的已打包输入表是否存在
  plot <label>            从已打包衍生表/源表/原始矢量重绘图件
                          （fig1/fig2/fig3/fig4/fig5/fig6/figS1/figS2/figS3/figS4）
  plot all                批量重绘全部十张图
  joint                   运行联合构型支持诊断（可移植 CLI，见 scripts/joint/）

退出码：0 成功；2 需要外部输入；1 错误。
"""
import argparse
import csv
import importlib.util
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
MANIFEST = "manifest/figure_to_analysis.tsv"
FIGMOD = REPO / "scripts" / "figures" / "manuscript_figures.py"

# 十张图所需的小体积衍生表 / 原始矢量 / 源树均已打包进 data/upstream
# （逐项见 manifest/figure_to_analysis.tsv 与 data/upstream/SOURCES.tsv）。
# Fig6 采用「原始矢量 + 52 轮已批准标签修正」链，几何与数据不变。
UNAVAILABLE = {}

_MOD = None


def _figures():
    """载入 scripts/figures/manuscript_figures.py（原始绘图代码移植版）。"""
    global _MOD
    if _MOD is None:
        spec = importlib.util.spec_from_file_location("manuscript_figures", FIGMOD)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        _MOD = mod
    return _MOD


def load_manifest(root):
    p = root / MANIFEST
    with p.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def cmd_list(root):
    print(f"{'label':<7}{'status':<26}{'copied_input_tables'}")
    for r in load_manifest(root):
        print(f"{r['figure_label']:<7}{r['status']:<26}{r['copied_figure_input_tables']}")
    return 0


def cmd_verify(root):
    bad = 0
    for r in load_manifest(root):
        cells = [c for c in r["copied_figure_input_tables"].split(";")
                 if c and not c.startswith("data/derived (none") and not c.startswith("none")]
        miss = [c for c in cells if not (root / c).exists()]
        state = "OK" if not miss else "MISSING:" + ",".join(miss)
        if miss:
            bad += 1
        print(f"[{state}] {r['figure_label']}  ({len(cells)} inputs)")
    print(f"\nverify: {bad} figure(s) with missing bundled inputs")
    return 0 if bad == 0 else 1


def cmd_joint(root, out):
    """运行 scripts/joint/joint_configuration_support.py（可移植联合构型诊断）。"""
    import runpy
    data_root = root / "data/upstream"
    if not data_root.is_dir():
        print(f"[error] 缺少已打包输入目录：{data_root}")
        return 1
    try:
        script = root / "scripts/joint/joint_configuration_support.py"
        mod = runpy.run_path(str(script))
        outdir = out / "joint"
        mod["run"](data_root, outdir)
        made = sorted(p.name for p in outdir.glob("*"))
        print(f"[ok] joint -> {os.path.relpath(outdir, root)}/ ({', '.join(made)})")
        return 0
    except Exception as e:  # noqa: BLE001
        print(f"[error] joint: {type(e).__name__}: {e}")
        return 1


def cmd_plot(root, out, label):
    out.mkdir(parents=True, exist_ok=True)
    plotters = _figures().PLOTTERS
    labels = list(plotters) if label == "all" else [label]
    rc = 0
    for lb in labels:
        if lb in UNAVAILABLE:
            print(f"[external-input-required] {lb}: {UNAVAILABLE[lb]}")
            rc = max(rc, 2)
            continue
        fn = plotters.get(lb)
        if not fn:
            print(f"[error] unknown figure label: {lb}")
            rc = max(rc, 1)
            continue
        try:
            made = fn(root, out)
            print(f"[ok] {lb} -> " + ", ".join(os.path.relpath(m, root) for m in made))
        except Exception as e:  # noqa: BLE001
            print(f"[error] {lb}: {type(e).__name__}: {e}")
            rc = max(rc, 1)
    return rc


def main(argv=None):
    ap = argparse.ArgumentParser(description="Mars-analogue manuscript v6 figure reproduction CLI")
    ap.add_argument("--input-root", default=str(REPO), help="repository root (default: repo containing this bin/)")
    ap.add_argument("--output-root", default=None, help="output dir (default: <input-root>/figures_out)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list")
    sub.add_parser("verify")
    sub.add_parser("joint")
    pp = sub.add_parser("plot")
    pp.add_argument("label")
    a = ap.parse_args(argv)
    root = Path(a.input_root).resolve()
    out = Path(a.output_root).resolve() if a.output_root else root / "figures_out"
    if a.cmd == "list":
        return cmd_list(root)
    if a.cmd == "verify":
        return cmd_verify(root)
    if a.cmd == "joint":
        return cmd_joint(root, out)
    if a.cmd == "plot":
        return cmd_plot(root, out, a.label)
    ap.error("unknown command")
    return 1


if __name__ == "__main__":
    sys.exit(main())
