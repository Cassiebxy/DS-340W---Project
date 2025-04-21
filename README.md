# Forward Mortality Modeling using Regional Death Data  
### Author: Xinyi Bao  
### Course: DS 340 W

## Overview

This project explores a forward-looking mortality forecasting framework by integrating principal component analysis (PCA) with stochastic differential equations (SDEs), specifically the Ornstein–Uhlenbeck process. It aims to improve model interpretability and forecasting performance in regional mortality settings without relying on traditional age-grouped data.

The work builds upon concepts presented in the paper **"AI in Actuarial Science – Part 2"** (Richman, 2021), with added emphasis on:
- Using **city-level monthly mortality data** instead of age-specific mortality
- Implementing a **dynamic simulation approach** to capture mortality shocks such as COVID-19

---

## Repository Structure(File Description)
```
.
├── Final Implement.R: Full implementation of this project's PCA–SDE mortality model
├── local_mortality.csv: Project dataset, Monthly regional death data (2015–2022) from 24 global subregions
├── AI_in_Actuarial_Science-master.zip: Source code from the parent paper (Richman 2021)
├── Keras - Mortality.R: Original mortality model from the parent paper
├── Keras - Mortality(re-do).R: Re-run version of the parent paper's code using my environment
├── Mx_1x1.txt: Original dataset from the parent paper
└── README.md: You're here!
```

---

## Dataset Description

- **local_mortality.csv** contains monthly all-cause death counts for 24 subnational regions from 2015 to 2022.
  - Countries include Argentina, India, and China.
  - Data sourced from the **World Mortality Dataset** GitHub repository.
  - Columns include: `year`, `time`, `deaths`, `country_name`, `local_unit_name`, etc.

---

## How to Run

1. Clone the repo or download this folder.
2. Open `Final Implement.R` in RStudio.
3. Install required packages:
   ```r
   install.packages(c("data.table", "lubridate", "ggplot2", "reshape2", "dplyr"))
   ```
4. Make sure `local_mortality.csv` is in the same directory.
5. Run the script section-by-section to:
   - Preprocess and filter data
   - Apply PCA to extract temporal mortality factors
   - Simulate future values using the OU process
   - Reconstruct and visualize regional mortality forecasts

---

## Results

- The project compares 1PC vs 2PC models.
- 2PC model shows improved fit in high-volatility regions (e.g., Tamil Nadu), while 1PC provides a cleaner fit in stable regions.
- Visualization and evaluation are provided in the script for further analysis.

---

## Citation

If using or referencing this work, please cite:
- Richman, R. (2021). *AI in Actuarial Science – Part 2*, Annals of Actuarial Science, 15(2), 230–258. https://doi.org/10.1017/S174849952000024X
- World Mortality Dataset: https://github.com/akarlinsky/world_mortality

---

For questions or contributions, feel free to open an issue or contact xzb5043@psu.edu.

