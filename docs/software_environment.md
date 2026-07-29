# Software environment

The validation run used:

| Component | Version |
|---|---:|
| R | 4.5.3 |
| readxl | 1.4.5 |
| lme4 | 2.0.1 |
| lmerTest | 3.2.1 |
| splines | 4.5.3 |
| zoo | 1.8.15 |
| sandwich | 3.1.1 |
| mice | 3.19.0 |
| mediation | 4.5.1 |
| MASS | 7.3-65 |
| ordinal | 2025.12-29 |
| geepack | 1.3-13 |
| ggplot2 | 4.0.2 |

Install dependencies explicitly before execution. The scripts intentionally do
not call `install.packages()` because automatic installation makes a clinical
analysis difficult to audit and may silently change package versions.

Required packages by analysis are checked at runtime. Each completed analysis
writes `session_info.txt`, input MD5 hashes, seed settings, and key sample
dimensions to its generated output directory.
