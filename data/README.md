# Data dictionary

One row per patient. Column names are converted to lower case by `R/01_prepare_data.R`; decimal commas are accepted. Times: `_0h` at sepsis diagnosis (admission), `_24h` and `_72h` after it.

| Variable | Description | Unit |
|---|---|---|
| `id` | Patient identifier | |
| `idade` (or `age`) | Age | years |
| `sex_masc` | Male sex | 0/1 |
| `bmi` | Body mass index | kg/m² |
| `saps_3` | SAPS 3 score | |
| `sofa_0h`, `sofa_72h` | SOFA score | |
| `fc_max_0h` | Maximum heart rate | beats/min |
| `ench_cap_0h` | Capillary refill time | s |
| `bh_cum_0h` | Cumulative fluid balance | mL |
| `vasopr_max_0h` | Maximum norepinephrine dose | µg·kg⁻¹·min⁻¹ |
| `lactate_0h`, `lactate_24h`, `lactate_72h` | Arterial lactate | mmol/L |
| `svo2_0h` | Central venous oxygen saturation | % |
| `gapco2_0h` | Venous-to-arterial CO₂ difference | mmHg |
| `troponine_0h` | Troponin | ng/L |
| `albumin_0h` | Albumin | g/dL |
| `pf_0h` | PaO₂/FiO₂ ratio | |
| `crp_0h`, `crp_72h` | C-reactive protein | mg/L |
| `leuc_0h`, `leuc_72h` | Total leukocyte count | cells/mm³ |
| `vm` | Mechanical ventilation | 0/1 |
| `corticoide_0h` | Corticosteroid therapy | 0/1 |
| `mort_uti` | ICU mortality | 0/1 |
| `mort_hosp` | In-hospital mortality | 0/1 |
| `mort_28d` | 28-day mortality | 0/1 |
| `mort_6m` | 6-month mortality | 0/1 |
| `il1_0h`, `il1_72h` | IL-1β | pg/mL |
| `il6_0h`, `il6_72h` | IL-6 | pg/mL |
| `il10_0h`, `il10_72h` | IL-10 | pg/mL |
| `tnfa_0h`, `tnfa_72h` | TNF-α | pg/mL |
| `gdf15_0h`, `gdf15_72h` | GDF15 | pg/mL |
| `fgf21_0h`, `fgf21_72h` | FGF21 | pg/mL |
| `leptin_0h`, `leptin_72h` | Leptin | pg/mL |
| `ghrelin_0h`, `ghrelin_72h` | Ghrelin | pg/mL |
| `r_0h`, `r_72h` | LEAK respiration in the absence of adenylates (Ln) | pmol O₂·s⁻¹·mg⁻¹ |
| `ci_0h`, `ci_72h` | Complex I-linked respiration | pmol O₂·s⁻¹·mg⁻¹ |
| `cii_0h`, `cii_72h` | Complex II-linked respiration | pmol O₂·s⁻¹·mg⁻¹ |
| `bce_0h`, `bce_72h` | Bioenergetic coupling efficiency, (P − Ln)/P | |

Renal replacement therapy and source of infection are summarised in Table 1; set their column names at the top of `R/02_table1.R`.
