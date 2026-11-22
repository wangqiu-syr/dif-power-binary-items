TABLE_SIM1
| Method | Power (3 DIF items) | Type I error (7 clean items) |
|---|---|---|
| MIMIC, baseline specification (z ≥ 1.96, one-sided) | .348 [.329, .366] | .026 [.022, .030] |
| IRT-LRT (difLRT, purified) | .323 [.305, .340] | .122 [.114, .130] |
| Rasch tree + anchored Wald test | .121 [.106, .137] | .041 [.035, .048] |

TABLE_EXACT
| | θ known | Item's own proportions (= baseline MIMIC) | Clean anchor | All other items as anchor |
|---|---|---|---|---|
| Power at 50 + 50, mean of i1, i5, i10 | .39 | .34 | .36 | .25 |
| Type I error, mean of the 7 clean items | .05 | .05 | .05 | .08 |
| n per group for .80 (range over i1, i5, i10) | 125–165 | 151–193 | 142–181 | 221–297 |

TABLE_TYPE1
| Method | No DIF, μ = 0 | No DIF, μ = 0.5 | One-way DIF 0.8, μ = 0 | One-way DIF 0.8, μ = 0.5 | Mixed DIF 0.8, μ = 0.5 |
|---|---|---|---|---|---|
| Benchmark: θ known | .049 / .052 / .057 / .050 | .048 / .049 / .054 / .053 | .050 / .054 / .057 / .053 | .044 / .055 / .046 / .045 | .047 / .051 / .052 / .048 |
| **LRT family** | | | | | |
| difLRT, purified (baseline; n = 50, μ = 0 only) | .089 / — / — / — | — | .121 / — / — / — | — | — |
| CML LRT, purified | .055 / .046 / .054 / .048 | .050 / .050 / .051 / .049 | .069 / .086 / .082 / .058 | .066 / .081 / .070 / .061 | .057 / .056 / .055 / .046 |
| Mantel–Haenszel, purified | .048 / .044 / .055 / .047 | .051 / .053 / .047 / .048 | .060 / .079 / .081 / .056 | .058 / .086 / .070 / .060 | .052 / .053 / .058 / .045 |
| CML LRT, robust anchor (k = 5) | .046 / .041 / .052 / .043 | .039 / .041 / .044 / .035 | .061 / .080 / .067 / .060 | .065 / .073 / .066 / .054 | .047 / .047 / .048 / .038 |
| CML LRT, true clean anchor (benchmark) | .059 / .051 / .058 / .053 | .054 / .053 / .054 / .050 | .050 / .058 / .055 / .050 | .056 / .057 / .051 / .050 | .051 / .049 / .048 / .049 |
| **MIMIC family** | | | | | |
| MIMIC, baseline specification | .046 / .052 / .054 / .055 | .169 / .291 / .416 / .603 | .055 / .061 / .050 / .050 | .172 / .282 / .429 / .631 | .171 / .297 / .407 / .596 |
| MIMIC with scale purification (M-SP) | .056 / .050 / .060 / .052 | .057 / .056 / .053 / .054 | .080 / .086 / .087 / .060 | .070 / .086 / .077 / .058 | .057 / .055 / .054 / .047 |
| **Rasch tree family** | | | | | |
| Rasch tree + anchortest (baseline) | .012 / .010 / .015 / .011 | .009 / .014 / .009 / .011 | .043 / .090 / .093 / .085 | .042 / .085 / .088 / .073 | .045 / .069 / .075 / .071 |
| anchortest without the tree's gate | .074 / .063 / .073 / .067 | .064 / .064 / .066 / .063 | .095 / .117 / .101 / .085 | .086 / .111 / .096 / .073 | .070 / .079 / .077 / .071 |

TABLE_POWER
| Method | One-way DIF 0.8, μ = 0 | One-way DIF 0.8, μ = 0.5 | Mixed DIF 0.8, μ = 0.5 | One-way DIF 0.5, μ = 0 | Mixed DIF 0.5, μ = 0.5 |
|---|---|---|---|---|---|
| Benchmark: θ known | .382 / .707 / .832 / .963 | .394 / .678 / .847 / .970 | .384 / .658 / .820 / .957 | .162 / .362 / .474 / .662 | .202 / .327 / .443 / .682 |
| **LRT family** | | | | | |
| difLRT, purified (baseline; n = 50, μ = 0 only) | .322 / — / — / — | — | — | — | — |
| CML LRT, purified | .238 / .484 / .683 / .901 | .268 / .469 / .721 / .920 | .362 / .604 / .767 / .923 | .110 / .232 / .297 / .502 | .197 / .311 / .418 / .613 |
| Mantel–Haenszel, purified | .213 / .496 / .670 / .902 | .233 / .466 / .683 / .908 | .333 / .568 / .742 / .912 | .084 / .219 / .291 / .491 | .170 / .297 / .397 / .600 |
| CML LRT, robust anchor (k = 5) | .216 / .446 / .632 / .879 | .236 / .433 / .669 / .897 | .330 / .549 / .744 / .904 | .103 / .216 / .276 / .454 | .163 / .273 / .374 / .577 |
| CML LRT, true clean anchor (benchmark) | .343 / .633 / .797 / .934 | .357 / .614 / .804 / .958 | .353 / .587 / .788 / .932 | .143 / .333 / .432 / .622 | .187 / .288 / .401 / .613 |
| **MIMIC family** | | | | | |
| MIMIC, baseline specification | .321 / .624 / .772 / .928 | .086 / .149 / .168 / .264 | .302 / .399 / .450 / .496 | .162 / .311 / .401 / .610 | .223 / .312 / .344 / .362 |
| MIMIC with scale purification (M-SP) | .246 / .476 / .676 / .898 | .238 / .478 / .667 / .901 | .352 / .553 / .729 / .894 | .107 / .238 / .298 / .521 | .193 / .303 / .398 / .591 |
| **Rasch tree family** | | | | | |
| Rasch tree + anchortest (baseline) | .125 / .432 / .661 / .911 | .140 / .408 / .670 / .930 | .242 / .510 / .766 / .924 | .032 / .135 / .217 / .428 | .067 / .181 / .340 / .575 |
| anchortest without the tree's gate | .267 / .524 / .713 / .916 | .303 / .523 / .728 / .939 | .404 / .614 / .787 / .924 | .134 / .273 / .343 / .541 | .217 / .352 / .459 / .637 |

MCSE_NOTE
largest MCSE: Type I 0.014 (no-DIF cells 0.012); power 0.024

TABLE_MATRIX
| Method | A1 | A2 | A3 | A4 | A5 | A6 | B1 | B2 | B3 | B4 | B5 | B6 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **LRT family** | | | | | | | | | | | | |
| difLRT, purified (baseline)ᵃ | — | Semi | Yes | Semi† | No | Yes | — | — | No | Semi | No | Semi |
| CML LRT, purified | Yes | Semi | Yes | Yes | Yes | Yes | Semi | Yes | Semi | Yes | No | Semi |
| Mantel–Haenszel, purified | Yes | Semi | Yes | Semi† | Yes | Yes | Semi | Yes | Semi | Yes | No | Semi |
| CML LRT, robust anchor (k = 5) | Yes | Semi | Yes | Yes | Yes | Yes | Semi | Yes | Semi | Yes | No | Semi |
| CML LRT, true clean anchor (benchmark) | Yes | Yes | Yes | Yes | Yes | Yes | Yes | Yes | Yes | Yes | No | Yes |
| **MIMIC family** | | | | | | | | | | | | |
| MIMIC, baseline specification | Semi | No | Yes | Semi† | Yes | Yes | No | No | Yes | Yes | No | No |
| MIMIC with scale purification (M-SP) | Yes | Semi | Yes | Semi† | Yes | Yes | Semi | Yes | Semi | Yes | No | Semi |
| **Rasch tree family** | | | | | | | | | | | | |
| Rasch tree + anchortest (baseline)ᵇ | Yes | Semi | Yes | Semi† | Semi | No | Semi | Yes | No | Semi | Yes | Semi |
| anchortest without the tree's gate | Yes | Semi | Yes | Semi† | Yes | Yes | No | Yes | No | Semi | No | Semi |
