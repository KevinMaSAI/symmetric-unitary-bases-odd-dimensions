# Orthogonal bases of symmetric unitary matrices in odd dimensions

Author / 作者: **Qiyuan Ma**, Shanghai Jiao Tong University (上海交通大学), ORCID [0009-0000-9532-0821](https://orcid.org/0009-0000-9532-0821).

This repository contains the public preprint and the verification/certification code for a **partial result**. It is a preprint and has not undergone academic peer review.

本仓库收录该预印本及其验证 / 认证代码。这是**部分结果**，尚未经过学术同行评审。

## Scope (important) / 范围（重要）

This paper does **not** prove the existence of a symmetric unitary orthogonal basis for all odd dimensions.

It proves, by an exact n = 5 construction and a computer-assisted Krawczyk / Arb argument, that such a basis exists for every odd integer **n with 3 ≤ n ≤ 67**. Consequently the conjecture of Girard, Leung, Levick, Li, Paulsen, Poon and Watrous holds for all n ≤ 68. The case of general odd n remains open.

本文**没有**解决全部奇数维情形，只证明了每个奇数 **3 ≤ n ≤ 67** 存在由对称酉矩阵组成的正交基底，因而上述猜想对所有 n ≤ 68 成立。一般奇数维情形仍然开放。

## Identifiers and files / 标识与文件

- Hexagon: [hexagon:2611.00505v1](https://hexagonmath.org/2611.00505v1)
- AiraXiv: [2611.0027](https://airaxiv.com/papers/view/2611.0027/)
- Zenodo Version 1.1: [10.5281/zenodo.24246318](https://doi.org/10.5281/zenodo.24246318)

The PDF and code here are byte-identical to the Zenodo Version 1.1 files `paper.pdf` and `code.zip` (MD5 `a9e6292ab58605d54ffc0ef7c29c6760` and `b33e1f70c6bbee1c9b25332b2cb2e697`). The manuscript header still reads “Research preprint, version 1.0” (15 November 2026); that is the text on the PDF, while the Zenodo record of these same files is Version 1.1. The TeX source is the `source.zip` that accompanies this deposit (`paper/main.tex` only; the bibliography is already inlined).

此处的 PDF 与代码与 Zenodo Version 1.1 的 `paper.pdf`、`code.zip` 逐字节一致。文稿题头仍写 “version 1.0”，这是 PDF 上的文字；同一批文件在 Zenodo 记录中的版本号是 1.1。

## License / 许可

- Paper (`paper/paper.pdf`, `paper/main.tex`): [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) (see `LICENSE-PAPER`).
- Code (`code/`): MIT License, already present as `code/LICENSE` in the original code archive. This repository does not relicense the code.

论文为 CC BY 4.0；代码沿用原归档中已有的 MIT 许可证，本仓库不另行再授权。

## AI disclosure / AI 使用说明

The author used AI assistance in the literature search, in designing and running the computational experiments, in writing the accompanying verification and certification code, and in drafting the manuscript. Tools used: Codex (OpenAI) and Claude (Anthropic), as research and writing assistants. Specific version numbers were not recorded. AI is not listed as an author. The human author is solely responsible for verifying all results, for the correctness of the proofs and computations, and for the final text.

作者在文献检索、计算实验的设计与运行、验证 / 认证代码的编写以及文稿起草中使用了 AI 辅助（Codex（OpenAI）与 Claude（Anthropic））。未记录具体版本号。AI 不是作者。人类作者对核实全部结果、证明与计算的正确性以及最终文本负责。

## Repository layout / 仓库结构

- `paper/paper.pdf` — compiled preprint (8 pages).
- `paper/main.tex` — LaTeX source corresponding to that PDF.
- `code/` — verification scripts, archived Krawczyk centres for odd n = 3, 5, …, 67, reference outputs and certification logs, plus `code/README.md` and `code/requirements.txt`.
- `CITATION.cff` — citation metadata.

Do not treat `code/logs/` as a complete proof by themselves; they are transcripts of the certification runs described in the paper and in `code/README.md`.

## How to reproduce the short checks / 如何复现短检查

See `code/README.md` for the full command list, including recertifying every odd n ≤ 67. A minimal check is:

```bash
cd code
python3 -m pip install -r requirements.txt
python3 verify_exact.py
python3 certify.py 5 --from-center
```

`verify_exact.py` should end with `ALL EXACT CHECKS PASSED`. The certification line should read `CERTIFIED: unique exact zero in box`. Recertifying all 33 archived centres, or running `python3 certify.py n` without `--from-center` (a new search), is optional and longer; those experiments were not re-run when this repository snapshot was prepared.

`verify_exact.py` 只需要 Python 标准库和 numpy；`certify.py` 需要 numpy、scipy 和 python-flint（含 Arb）。测试环境见 `code/README.md`。

## Status / 状态

This is an unrefereed preprint of a **finite-range** computer-assisted result. Please do not describe it as a proof for all odd n.

这是未经同行评审的预印本，结果只覆盖有限的奇数范围，请不要将其说成对全部奇数维的证明。
