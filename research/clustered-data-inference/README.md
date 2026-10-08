# Confidence intervals with clustered data

I examined how treating clustered observations as independent affects confidence intervals and false-positive rates. I generated balanced Gaussian data with treatment assigned to whole clusters. I compared individual-level OLS with an analysis of cluster means across sixteen scenarios, with 2,000 replications each. With twenty clusters of thirty observations and ICC = 0.30 under the null, nominal 95% interval coverage was 46.5% for individual OLS and 95.2% for cluster means.

## Results

| Method | 95% interval coverage | False-positive rate |
| --- | ---: | ---: |
| Individual-level OLS | 46.5% | 53.5% |
| Cluster means | 95.2% | 4.8% |

This table shows the null scenario with twenty clusters, thirty observations per cluster, and ICC = 0.30. The values come from [summary.csv](results/summary.csv). The Monte Carlo SE for cluster-mean coverage is about 0.48 percentage points.

![Interval coverage across null scenarios](results/coverage.png)

I varied cluster count (20 or 50), cluster size (10 or 30), ICC (0.05 or 0.30), and treatment effect (0 or 0.20 marginal SD). I assigned half the clusters to treatment and used seed 20261008. [design.csv](results/design.csv) gives all sixteen scenarios. The summary reports bias, RMSE, empirical SD, estimated SE, coverage, rejection rates, and Monte Carlo SEs.

Both methods estimate the same difference in means in this balanced design. Their standard errors and degrees of freedom differ. I use a pooled-variance t interval for cluster means, with degrees of freedom equal to the number of clusters minus two. Plot bars show ±1.96 Monte Carlo SE.

<details>
<summary>Code and files</summary>

I created and ran this simulated-data project on 8 October 2026. It uses generated data, not thesis or participant records. It requires base R only.

```bash
Rscript analysis.R
```

For a shorter run, I can set the number of replications:

```bash
Rscript analysis.R 500
```

The script checks its coefficient and SE calculations against `lm()` on the first replication.

```text
clustered-data-inference/
├── .gitignore
├── README.md
├── analysis.R
└── results/
    ├── coverage.png
    ├── design.csv
    ├── run_metadata.txt
    ├── session_info.txt
    └── summary.csv
```

</details>

## Limitations

I studied equal-sized clusters with Gaussian outcomes and independent random intercepts. The results do not cover unequal cluster sizes, observational exposures, missing data, or treatment assigned within clusters. Two thousand replications leave Monte Carlo uncertainty, which I report alongside each estimate.

## Credits

Project author: **Laura Maria Fetz**. This is a new methods project, separate from my earlier thesis and published work.
