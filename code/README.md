# Reproducibility code

Supporting files for the preprint

> Orthogonal bases of symmetric unitary matrices in odd dimensions and the mixed-unitary rank of the symmetric Werner–Holevo channel

This archive is intended for a Zenodo deposit together with `paper.pdf`. Author: Qiyuan Ma (Shanghai Jiao Tong University, ORCID 0009-0000-9532-0821).

## Contents

| Path | Role |
| --- | --- |
| `verify_exact.py` | Exact-arithmetic verification of the explicit $n=5$ basis (Python `fractions` only; `numpy` is used only for a numeric rank sanity check at the end). |
| `verify_sympy.py` | Independent Gröbner-basis verification of the same $n=5$ basis. |
| `reduce_eqs.py` | Expansions of $E_{12}$ and $E_{23}$ used in Step 3 of the $n=5$ proof. |
| `certify.py` | Krawczyk / Arb certification of a zero of the square system $F_n$ for odd $n$. |
| `centers/cert_sol_n*.npy` | Archived box centres for every odd $n$ with $3\le n\le 67$ (the centres used for Table 1 of the paper). |
| `logs/` | Original certification logs, the table recertification log `recertify_from_centers.out`, `package_rerun.out` from the packaging check described below, and `review_recertify_all.out` (full re-run of all 33 cases from this tree during the final review). |
| `verify_exact.out`, `verify_sympy.out` | Reference outputs from the original runs. |
| `LICENSE` | MIT license for this code. |
| `requirements.txt` | Python dependencies. |

## Requirements

Tested on 15 November 2026 (Asia/Shanghai) with:

- Python 3.13.5
- python-flint 0.9.0 (FLINT 3.6.0, including Arb ball arithmetic)
- numpy 2.5.4
- scipy 1.18.1
- sympy 1.14.0 (needed only for `verify_sympy.py` and `reduce_eqs.py`)

Install, for example:

```bash
python3 -m pip install -r requirements.txt
```

`verify_exact.py` needs only the Python standard library plus numpy (the numpy import is used solely for a floating-point rank check and is not part of the exact proof). `certify.py` needs numpy, scipy, and python-flint.

## Commands

Unpack and enter the directory:

```bash
unzip code.zip
cd sym-unitary-code
```

### 1. Exact $n=5$ verification (about 1 second)

```bash
python3 verify_exact.py
```

Expected terminal output ends with:

```
ALL EXACT CHECKS PASSED
```

All symmetry, unitarity, Gram-matrix and lemma checks should print `True`, and the Gram check should report `failures = 0`. A reference transcript is in `verify_exact.out`.

### 2. Recertify archived centres (the proof of Table 1)

```bash
python3 certify.py n --from-center
```

where `n` is an odd integer from 3 to 67. This loads `centers/cert_sol_n{n}.npy` and evaluates the Krawczyk operator in 128-bit Arb arithmetic on the box of radius $10^{-9}$. It does **not** search for a new centre.

Expected line (numbers must match Table 1 of the paper):

```
n={n}: Krawczyk box radius 1e-09, max|K-x~|=...  ->  CERTIFIED: unique exact zero in box
```

To recertify every odd $n\le 67$:

```bash
for n in $(seq 3 2 67); do python3 certify.py $n --from-center; done
```

A full original transcript is `logs/recertify_from_centers.out`.

### 3. Optional independent $n=5$ check (about 50 seconds)

```bash
python3 verify_sympy.py
```

Expected output:

```
symmetric & unitary: True
Gram == 5I failures: 0
ALL PASSED
```

### 4. Optional search + certify (does not use the archived centre)

```bash
python3 certify.py n
```

This searches from a random starting point and, if it finds a floating-point approximate zero, certifies it. A new `cert_sol_n{n}.npy` is written in the current working directory. For reproducing Table 1, use `--from-center` instead.

### 5. Optional expansions for the $n=5$ proof, Step 3

```bash
python3 reduce_eqs.py
```

## Packaging check (ran from this tree, not uploaded)

The following commands were run from a copy of this directory on 15 November 2026, 18:55:50 +0800 (Asia/Shanghai), with the versions listed above. Full transcript: `logs/package_rerun.out`.

```bash
python3 verify_exact.py
python3 certify.py 5 --from-center
python3 certify.py 9 --from-center
python3 certify.py 67 --from-center
```

Results:

- `verify_exact.py`: `ALL EXACT CHECKS PASSED`.
- $n=5$: `max|K-x~|=5.619e-15` → `CERTIFIED: unique exact zero in box` (matches Table 1).
- $n=9$: `max|K-x~|=5.107e-15` → `CERTIFIED: unique exact zero in box` (matches Table 1).
- $n=67$: `max|K-x~|=1.888e-12` → `CERTIFIED: unique exact zero in box` (matches Table 1; wall-clock time was 6.6 s on this machine; the table lists 5.7 s, and timings vary).

The residual `max|F|` and `cond(J)` printed on the first line of `certify.py` are floating-point diagnostics and may vary slightly; the certificate is the `CERTIFIED` line and the rigorous `max|K-x~|` bound.

The validity of the computer-assisted proof depends on the correctness of Arb's ball arithmetic and of Krawczyk's theorem. Floating-point arithmetic is used only to produce the centre and the preconditioner $Y$, neither of which needs to be accurate for the certificate to be valid.

## License

This code is released under the MIT License (see `LICENSE`). The paper PDF is released under CC-BY-4.0.

## Final review re-run

On 15 November 2026 (Asia/Shanghai) all 33 archived centres (odd n from 3 to 67) were re-certified from this tree with `certify.py n --from-center`; all printed `CERTIFIED` and every `max|K-x~|` value matches Table 1. Transcript: `logs/review_recertify_all.out`.
