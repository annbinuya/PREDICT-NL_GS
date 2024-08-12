# Prognostication and treatment predictions for early-stage breast cancer: incorporating the 70-gene signature into the PREDICT prognostication model

The 70-gene signature (70-GS) has been shown to identify women at low-risk of distant recurrence who can safely forgo adjuvant chemotherapy. Incorporating this GS into the well-validated and widely used PREDICT breast cancer model could improve the model’s ability to estimate breast cancer prognosis, and thereby further reduce overtreatment and its long-term impact on patients’ quality of life. We incorporated the 70-GS into PREDICT-v2.3 and assessed the new PREDICT-GS model’s ability to predict 5-year risk of breast cancer death.

The prognostic effect of 70-GS was obtained from the MINDACT study (for details, see: [https://github.com/annbinuya/PREDICT-NL_Methods])

# Syntax files
| File                   | Description             |
| :----                  | :----                   |
| 1_Cleaning.Rmd                     | Cleaning, preparation, and imputation of the Netherlands Cancer Registry (here, "NKR") dataset.
| 2_NKR_Descriptive_imp.Rmd          | Descriptive tables, imputed dataset.
| 2_NKR_Descriptive_raw.Rmd          | Descriptive tables, raw dataset.
| 3_PREDICTv23_validation.Rmd        | Validation of the PREDICTv2.3 model.
| 4_PREDICT-GS_validation.Rmd        | Validation of the PREDICT-GS model.
| 5_PREDICT_P-O.Rmd                  | Observed and predicted events at 5 years.
| pool_perf.R                        | Pooling function and other functions used.

# Contact
Mary Ann E. Binuya <br/>
Netherlands Cancer Institute <br/>
[m.binuya@nki.nl](m.binuya@nki.nl)

# Authors
| Author                 | Role   | Description             |
| :----                  | :----: | :----                   |
| Mary Ann Binuya   | Author | Development and support |
| Ellen Engelhardt  | Author | Development and review  |
| Paul Pharoah      | Author | Review  |
| Coralie Poncet    | Author | Review  |
| Emiel Rutgers     | Author | Review  |
| Martine Piccart   | Author | Review  |
| Fatima Cargoso    | Author | Review  |
| Laura van 't Veer | Author | Review  |
| Ewout Steyerberg  | Author | Development and review  |
| Sabine Linn       | Author | Development and review  |
| Marjanka Schmidt  | Author | Development and review   |
