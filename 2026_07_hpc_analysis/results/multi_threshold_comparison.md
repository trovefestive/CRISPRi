## Gene-Level Coverage Analysis: Multi-Threshold Comparison

### Set A Coverage by Threshold

| Metric | ≥300× (Count) | ≥300× (%) | ≥200× (Count) | ≥200× (%) | ≥100× (Count) | ≥100× (%) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Total Target Genes | 18,710 | 100% | 18,710 | 100% | 18,710 | 100% |
| All 3 gRNAs above threshold | 14,557 | 77.80% | 16,338 | 87.32% | 17,843 | 95.37% |
| Exactly 2 gRNAs above threshold | 3,573 | 19.10% | 2,136 | 11.42% | 789 | 4.22% |
| Exactly 1 gRNA above threshold | 530 | 2.83% | 221 | 1.18% | 74 | 0.40% |
| 0 gRNAs above threshold | 48 | 0.26% | 13 | 0.07% | 2 | 0.01% |
| Cumulative: ≥2 gRNAs above threshold | 18,130 | 96.90% | 18,474 | 98.74% | 18,632 | 99.58% |

### Set B Coverage by Threshold

| Metric | ≥300× (Count) | ≥300× (%) | ≥200× (Count) | ≥200× (%) | ≥100× (Count) | ≥100× (%) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Total Target Genes | 18,708 | 100% | 18,708 | 100% | 18,708 | 100% |
| All 3 gRNAs above threshold | 13,428 | 71.78% | 15,424 | 82.45% | 17,295 | 92.45% |
| Exactly 2 gRNAs above threshold | 4,294 | 22.95% | 2,835 | 15.15% | 1,270 | 6.79% |
| Exactly 1 gRNA above threshold | 882 | 4.71% | 418 | 2.23% | 142 | 0.76% |
| 0 gRNAs above threshold | 103 | 0.55% | 30 | 0.16% | 0 | 0.00% |
| Cumulative: ≥2 gRNAs above threshold | 17,722 | 94.73% | 18,259 | 97.60% | 18,565 | 99.24% |

### Key Metrics Comparison

| Set | Metric | ≥300× | ≥200× | ≥100× | Improvement (300→100) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| Set A | All 3 gRNAs above threshold | 77.80% | 87.32% | 95.37% | +17.56% |
| Set A | ≥2 gRNAs above threshold | 96.90% | 98.74% | 99.58% | +2.68% |
| Set A | 0 gRNAs above threshold (dropout) | 0.26% | 0.07% | 0.01% | -0.25% |
| Set B | All 3 gRNAs above threshold | 71.78% | 82.45% | 92.45% | +20.67% |
| Set B | ≥2 gRNAs above threshold | 94.73% | 97.60% | 99.24% | +4.51% |
| Set B | 0 gRNAs above threshold (dropout) | 0.55% | 0.16% | 0.00% | -0.55% |
