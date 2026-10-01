# analysis/ — 作者端联合构型诊断原始文件（原样保留）

本目录保留 root host 新增的联合构型支持审计**原始**文件，供溯源：

| 文件 | 说明 |
|---|---|
| `analysis_contract.md` | 冻结的科学判据（4 对 Cox/Cut 共检、门槛 2/3/5、group≥3；描述性，无 p 值/无环境模型） |
| `joint_configuration_support.py` | 作者端脚本；以 `--revision-root` 相对定位各轮 `host/` 目录 |

**公开复现入口**为去私有化、可移植的 `../scripts/joint/joint_configuration_support.py`
（以 `--data-root`，默认 `data/upstream`，定位已打包输入），其输出与
`../data/upstream/joint/reference_outputs/` **逐字节一致**。原始脚本中的相对路径
（`49_.../host/...`、`48_.../host/...`）在本仓库内**不可直接运行**，仅作溯源。
