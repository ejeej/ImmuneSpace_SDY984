---
title: "Project_ImmuneSpace"
author: "Shakir Suleimanov"
date: "2022-11-27"
output: 
  html_document:
    keep_md: TRUE
---




```r
library(dplyr)
library(purrr)
library(gtsummary)
library(pheatmap)
```

Загрузим данные


```r
all_noNorm_withResponse_eset <- readRDS("all_noNorm_withResponse_eset.rds")

all_data <- all_noNorm_withResponse_eset@phenoData@data
```

Отфильтруем данные по VZ


```r
sdy984 <- all_data %>% 
  filter(study_accession == "SDY984")
```

Посмотрим на структуру данных


```r
str(sdy984)
```

```
## 'data.frame':	140 obs. of  55 variables:
##  $ participant_id                : chr  "SUB174087.984" "SUB174135.984" "SUB174114.984" "SUB174093.984" ...
##  $ age_imputed                   : num  60 60 60 60 60 60 60 60 60 60 ...
##  $ gender                        : chr  "Female" "Female" "Female" "Male" ...
##  $ race                          : chr  "White" "White" "White" "White" ...
##  $ ethnicity                     : chr  "Not Hispanic or Latino" "Not Hispanic or Latino" "Not Hispanic or Latino" "Not Hispanic or Latino" ...
##  $ study_accession               : chr  "SDY984" "SDY984" "SDY984" "SDY984" ...
##  $ vaccine                       : chr  "Zostavax" "Zostavax" "Zostavax" "Zostavax" ...
##  $ vaccine_type                  : chr  "Live virus" "Live virus" "Live virus" "Live virus" ...
##  $ adjuvant                      : chr  "VZV" "VZV" "VZV" "VZV" ...
##  $ pathogen                      : chr  "Varicella Zoster" "Varicella Zoster" "Varicella Zoster" "Varicella Zoster" ...
##  $ arm_accession                 : chr  "ARM3537" "ARM3537" "ARM3537" "ARM3537" ...
##  $ uid                           : chr  "SUB174087.984_0_Days_BS934468" "SUB174135.984_0_Days_BS934516" "SUB174114.984_0_Days_BS934495" "SUB174093.984_0_Days_BS934474" ...
##  $ biosample_accession           : chr  "BS934468" "BS934516" "BS934495" "BS934474" ...
##  $ study_time_collected          : num  0 0 0 0 0 0 0 0 0 0 ...
##  $ study_time_collected_unit     : chr  "Days" "Days" "Days" "Days" ...
##  $ time_post_last_vax            : num  0 0 0 0 0 0 0 0 0 0 ...
##  $ unit_post_last_vax            : chr  "Days" "Days" "Days" "Days" ...
##  $ age_reported                  : num  60 60 60 60 60 60 60 60 60 60 ...
##  $ exposure_material_reported    : chr  "Zostavax" "Zostavax" "Zostavax" "Zostavax" ...
##  $ exposure_process_preferred    : chr  "vaccination" "vaccination" "vaccination" "vaccination" ...
##  $ matrix                        : chr  "SDY984_PBMC_Elderly_Geo" "SDY984_PBMC_Elderly_Geo" "SDY984_PBMC_Elderly_Geo" "SDY984_PBMC_Elderly_Geo" ...
##  $ gsm                           : chr  "GSM2094441" "GSM2094443" "GSM2094445" "GSM2094447" ...
##  $ Hispanic                      : num  0 0 0 0 0 0 0 0 0 1 ...
##  $ White                         : num  1 1 1 1 1 1 1 1 1 1 ...
##  $ Asian                         : num  0 0 0 0 0 0 0 0 0 0 ...
##  $ Black                         : num  0 0 0 0 0 0 0 0 0 0 ...
##  $ cell_type                     : chr  "PBMC" "PBMC" "PBMC" "PBMC" ...
##  $ cohort                        : chr  "elderly" "elderly" "elderly" "elderly" ...
##  $ featureSetName                : chr  "HGU133_plus_PM" "HGU133_plus_PM" "HGU133_plus_PM" "HGU133_plus_PM" ...
##  $ featureSetName2               : chr  "HGU133_plus_PM" "HGU133_plus_PM" "HGU133_plus_PM" "HGU133_plus_PM" ...
##  $ featureSetVendor              : chr  "Affymetrix" "Affymetrix" "Affymetrix" "Affymetrix" ...
##  $ geBatchName                   : chr  "SDY984" "SDY984" "SDY984" "SDY984" ...
##  $ irpBatchName                  : chr  "SDY984_Zostavax" "SDY984_Zostavax" "SDY984_Zostavax" "SDY984_Zostavax" ...
##  $ y_chrom_present_timepoint     : logi  FALSE TRUE FALSE FALSE TRUE FALSE ...
##  $ y_chrom_present               : logi  FALSE FALSE FALSE TRUE FALSE FALSE ...
##  $ failedYchromQC                : logi  FALSE TRUE FALSE TRUE TRUE FALSE ...
##  $ assay                         : chr  "elisa" "elisa" "elisa" "elisa" ...
##  $ MFC                           : num  -0.07932 1.40438 0.3132 0.00766 3.14552 ...
##  $ maxRBA                        : num  NA NA NA NA NA NA NA NA NA NA ...
##  $ maxStrain_MFC                 : chr  "IgG" "IgG" "IgG" "IgG" ...
##  $ maxStrain_RBA                 : chr  NA NA NA NA ...
##  $ ImmResp_baseline_value_MFC    : num  13.8 13.2 13.3 12.8 11.9 ...
##  $ ImmResp_baseline_timepoint_MFC: num  0 0 0 0 0 0 0 0 0 0 ...
##  $ ImmResp_postVax_value_MFC     : num  13.7 14.6 13.7 12.8 15.1 ...
##  $ ImmResp_postVax_timepoint_MFC : num  30 30 30 30 30 30 30 30 30 30 ...
##  $ ImmResp_baseline_value_RBA    : num  NA NA NA NA NA NA NA NA NA NA ...
##  $ ImmResp_baseline_timepoint_RBA: num  NA NA NA NA NA NA NA NA NA NA ...
##  $ ImmResp_postVax_value_RBA     : num  NA NA NA NA NA NA NA NA NA NA ...
##  $ ImmResp_postVax_timepoint_RBA : num  NA NA NA NA NA NA NA NA NA NA ...
##  $ MFC_p30                       : Factor w/ 3 levels "lowResponder",..: 1 3 1 1 3 1 2 2 2 1 ...
##  $ MFC_p40                       : Factor w/ 3 levels "lowResponder",..: 1 3 1 1 3 1 1 2 2 1 ...
##  $ maxRBA_p30                    : Factor w/ 3 levels "lowResponder",..: NA NA NA NA NA NA NA NA NA NA ...
##  $ maxRBA_p40                    : Factor w/ 3 levels "lowResponder",..: NA NA NA NA NA NA NA NA NA NA ...
##  $ maxRBA_p50                    : Factor w/ 3 levels "lowResponder",..: NA NA NA NA NA NA NA NA NA NA ...
##  $ numAssays                     : int  1 1 1 1 1 1 1 1 1 1 ...
```

```r
#summary(yellow_fever_data)
```

Сделаем Table 1 по полу


```r
sdy984_description_table_gender <- sdy984 %>%
  mutate(across(c(gender, race, ethnicity, study_accession, vaccine, vaccine_type, adjuvant,     pathogen, cell_type, cohort, featureSetVendor, assay, MFC_p30), ~ as.factor(.x))) %>%
  select(gender, race, ethnicity, cohort, MFC, ImmResp_baseline_value_MFC, MFC_p30) %>%
  tbl_summary(by = gender)

sdy984_description_table_gender
```

```{=html}
<div id="kfyfeilmsl" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#kfyfeilmsl .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#kfyfeilmsl .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#kfyfeilmsl .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#kfyfeilmsl .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#kfyfeilmsl .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#kfyfeilmsl .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#kfyfeilmsl .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#kfyfeilmsl .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#kfyfeilmsl .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#kfyfeilmsl .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#kfyfeilmsl .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#kfyfeilmsl .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#kfyfeilmsl .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#kfyfeilmsl .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#kfyfeilmsl .gt_from_md > :first-child {
  margin-top: 0;
}

#kfyfeilmsl .gt_from_md > :last-child {
  margin-bottom: 0;
}

#kfyfeilmsl .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#kfyfeilmsl .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#kfyfeilmsl .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#kfyfeilmsl .gt_row_group_first td {
  border-top-width: 2px;
}

#kfyfeilmsl .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#kfyfeilmsl .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#kfyfeilmsl .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#kfyfeilmsl .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#kfyfeilmsl .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#kfyfeilmsl .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#kfyfeilmsl .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#kfyfeilmsl .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#kfyfeilmsl .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#kfyfeilmsl .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-left: 4px;
  padding-right: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#kfyfeilmsl .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#kfyfeilmsl .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#kfyfeilmsl .gt_left {
  text-align: left;
}

#kfyfeilmsl .gt_center {
  text-align: center;
}

#kfyfeilmsl .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#kfyfeilmsl .gt_font_normal {
  font-weight: normal;
}

#kfyfeilmsl .gt_font_bold {
  font-weight: bold;
}

#kfyfeilmsl .gt_font_italic {
  font-style: italic;
}

#kfyfeilmsl .gt_super {
  font-size: 65%;
}

#kfyfeilmsl .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 75%;
  vertical-align: 0.4em;
}

#kfyfeilmsl .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#kfyfeilmsl .gt_indent_1 {
  text-indent: 5px;
}

#kfyfeilmsl .gt_indent_2 {
  text-indent: 10px;
}

#kfyfeilmsl .gt_indent_3 {
  text-indent: 15px;
}

#kfyfeilmsl .gt_indent_4 {
  text-indent: 20px;
}

#kfyfeilmsl .gt_indent_5 {
  text-indent: 25px;
}
</style>
<table class="gt_table">
  
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="&lt;strong&gt;Characteristic&lt;/strong&gt;"><strong>Characteristic</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;strong&gt;Female&lt;/strong&gt;, N = 84&lt;sup class=&quot;gt_footnote_marks&quot;&gt;1&lt;/sup&gt;"><strong>Female</strong>, N = 84<sup class="gt_footnote_marks">1</sup></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;strong&gt;Male&lt;/strong&gt;, N = 56&lt;sup class=&quot;gt_footnote_marks&quot;&gt;1&lt;/sup&gt;"><strong>Male</strong>, N = 56<sup class="gt_footnote_marks">1</sup></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="label" class="gt_row gt_left">race</td>
<td headers="stat_1" class="gt_row gt_center"></td>
<td headers="stat_2" class="gt_row gt_center"></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Black or African American</td>
<td headers="stat_1" class="gt_row gt_center">24 (29%)</td>
<td headers="stat_2" class="gt_row gt_center">12 (21%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">0 (0%)</td>
<td headers="stat_2" class="gt_row gt_center">8 (14%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    White</td>
<td headers="stat_1" class="gt_row gt_center">60 (71%)</td>
<td headers="stat_2" class="gt_row gt_center">36 (64%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">ethnicity</td>
<td headers="stat_1" class="gt_row gt_center"></td>
<td headers="stat_2" class="gt_row gt_center"></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Hispanic or Latino</td>
<td headers="stat_1" class="gt_row gt_center">4 (4.8%)</td>
<td headers="stat_2" class="gt_row gt_center">8 (14%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Not Hispanic or Latino</td>
<td headers="stat_1" class="gt_row gt_center">80 (95%)</td>
<td headers="stat_2" class="gt_row gt_center">48 (86%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">cohort</td>
<td headers="stat_1" class="gt_row gt_center"></td>
<td headers="stat_2" class="gt_row gt_center"></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    elderly</td>
<td headers="stat_1" class="gt_row gt_center">48 (57%)</td>
<td headers="stat_2" class="gt_row gt_center">28 (50%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    young</td>
<td headers="stat_1" class="gt_row gt_center">36 (43%)</td>
<td headers="stat_2" class="gt_row gt_center">28 (50%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">MFC</td>
<td headers="stat_1" class="gt_row gt_center">0.67 (0.28, 1.36)</td>
<td headers="stat_2" class="gt_row gt_center">0.45 (0.31, 0.92)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">ImmResp_baseline_value_MFC</td>
<td headers="stat_1" class="gt_row gt_center">13.58 (13.06, 14.13)</td>
<td headers="stat_2" class="gt_row gt_center">13.22 (12.85, 14.13)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">MFC_p30</td>
<td headers="stat_1" class="gt_row gt_center"></td>
<td headers="stat_2" class="gt_row gt_center"></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    lowResponder</td>
<td headers="stat_1" class="gt_row gt_center">28 (33%)</td>
<td headers="stat_2" class="gt_row gt_center">16 (29%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    moderateResponder</td>
<td headers="stat_1" class="gt_row gt_center">24 (29%)</td>
<td headers="stat_2" class="gt_row gt_center">28 (50%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    highResponder</td>
<td headers="stat_1" class="gt_row gt_center">32 (38%)</td>
<td headers="stat_2" class="gt_row gt_center">12 (21%)</td></tr>
  </tbody>
  
  <tfoot class="gt_footnotes">
    <tr>
      <td class="gt_footnote" colspan="3"><sup class="gt_footnote_marks">1</sup> n (%); Median (IQR)</td>
    </tr>
  </tfoot>
</table>
</div>
```

Сделаем Table 1 по возрасту


```r
sdy984_description_table_age <- sdy984 %>%
  mutate(across(c(gender, race, ethnicity, study_accession, vaccine, vaccine_type, adjuvant, pathogen, cell_type, cohort, featureSetVendor, assay, MFC_p30), ~ as.factor(.x))) %>%
  select(!where(is.character)) %>%
  select(c(gender, race, ethnicity, cohort, MFC, ImmResp_baseline_value_MFC, MFC_p30)) %>%
  tbl_summary(by = cohort)

sdy984_description_table_age
```

```{=html}
<div id="geywegafeq" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#geywegafeq .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#geywegafeq .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#geywegafeq .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#geywegafeq .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#geywegafeq .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#geywegafeq .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#geywegafeq .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#geywegafeq .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#geywegafeq .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#geywegafeq .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#geywegafeq .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#geywegafeq .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#geywegafeq .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#geywegafeq .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#geywegafeq .gt_from_md > :first-child {
  margin-top: 0;
}

#geywegafeq .gt_from_md > :last-child {
  margin-bottom: 0;
}

#geywegafeq .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#geywegafeq .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#geywegafeq .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#geywegafeq .gt_row_group_first td {
  border-top-width: 2px;
}

#geywegafeq .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#geywegafeq .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#geywegafeq .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#geywegafeq .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#geywegafeq .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#geywegafeq .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#geywegafeq .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#geywegafeq .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#geywegafeq .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#geywegafeq .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-left: 4px;
  padding-right: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#geywegafeq .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#geywegafeq .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#geywegafeq .gt_left {
  text-align: left;
}

#geywegafeq .gt_center {
  text-align: center;
}

#geywegafeq .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#geywegafeq .gt_font_normal {
  font-weight: normal;
}

#geywegafeq .gt_font_bold {
  font-weight: bold;
}

#geywegafeq .gt_font_italic {
  font-style: italic;
}

#geywegafeq .gt_super {
  font-size: 65%;
}

#geywegafeq .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 75%;
  vertical-align: 0.4em;
}

#geywegafeq .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#geywegafeq .gt_indent_1 {
  text-indent: 5px;
}

#geywegafeq .gt_indent_2 {
  text-indent: 10px;
}

#geywegafeq .gt_indent_3 {
  text-indent: 15px;
}

#geywegafeq .gt_indent_4 {
  text-indent: 20px;
}

#geywegafeq .gt_indent_5 {
  text-indent: 25px;
}
</style>
<table class="gt_table">
  
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="&lt;strong&gt;Characteristic&lt;/strong&gt;"><strong>Characteristic</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;strong&gt;elderly&lt;/strong&gt;, N = 76&lt;sup class=&quot;gt_footnote_marks&quot;&gt;1&lt;/sup&gt;"><strong>elderly</strong>, N = 76<sup class="gt_footnote_marks">1</sup></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1" scope="col" id="&lt;strong&gt;young&lt;/strong&gt;, N = 64&lt;sup class=&quot;gt_footnote_marks&quot;&gt;1&lt;/sup&gt;"><strong>young</strong>, N = 64<sup class="gt_footnote_marks">1</sup></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="label" class="gt_row gt_left">gender</td>
<td headers="stat_1" class="gt_row gt_center"></td>
<td headers="stat_2" class="gt_row gt_center"></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Female</td>
<td headers="stat_1" class="gt_row gt_center">48 (63%)</td>
<td headers="stat_2" class="gt_row gt_center">36 (56%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Male</td>
<td headers="stat_1" class="gt_row gt_center">28 (37%)</td>
<td headers="stat_2" class="gt_row gt_center">28 (44%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">race</td>
<td headers="stat_1" class="gt_row gt_center"></td>
<td headers="stat_2" class="gt_row gt_center"></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Black or African American</td>
<td headers="stat_1" class="gt_row gt_center">4 (5.3%)</td>
<td headers="stat_2" class="gt_row gt_center">32 (50%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Unknown</td>
<td headers="stat_1" class="gt_row gt_center">0 (0%)</td>
<td headers="stat_2" class="gt_row gt_center">8 (12%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    White</td>
<td headers="stat_1" class="gt_row gt_center">72 (95%)</td>
<td headers="stat_2" class="gt_row gt_center">24 (38%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">ethnicity</td>
<td headers="stat_1" class="gt_row gt_center"></td>
<td headers="stat_2" class="gt_row gt_center"></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Hispanic or Latino</td>
<td headers="stat_1" class="gt_row gt_center">8 (11%)</td>
<td headers="stat_2" class="gt_row gt_center">4 (6.2%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    Not Hispanic or Latino</td>
<td headers="stat_1" class="gt_row gt_center">68 (89%)</td>
<td headers="stat_2" class="gt_row gt_center">60 (94%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">MFC</td>
<td headers="stat_1" class="gt_row gt_center">0.40 (0.13, 0.67)</td>
<td headers="stat_2" class="gt_row gt_center">0.93 (0.59, 1.40)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">ImmResp_baseline_value_MFC</td>
<td headers="stat_1" class="gt_row gt_center">13.19 (12.82, 13.77)</td>
<td headers="stat_2" class="gt_row gt_center">13.65 (12.98, 14.35)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">MFC_p30</td>
<td headers="stat_1" class="gt_row gt_center"></td>
<td headers="stat_2" class="gt_row gt_center"></td></tr>
    <tr><td headers="label" class="gt_row gt_left">    lowResponder</td>
<td headers="stat_1" class="gt_row gt_center">32 (42%)</td>
<td headers="stat_2" class="gt_row gt_center">12 (19%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    moderateResponder</td>
<td headers="stat_1" class="gt_row gt_center">28 (37%)</td>
<td headers="stat_2" class="gt_row gt_center">24 (38%)</td></tr>
    <tr><td headers="label" class="gt_row gt_left">    highResponder</td>
<td headers="stat_1" class="gt_row gt_center">16 (21%)</td>
<td headers="stat_2" class="gt_row gt_center">28 (44%)</td></tr>
  </tbody>
  
  <tfoot class="gt_footnotes">
    <tr>
      <td class="gt_footnote" colspan="3"><sup class="gt_footnote_marks">1</sup> n (%); Median (IQR)</td>
    </tr>
  </tfoot>
</table>
</div>
```

Создание таблицы по данным экспрессии


```r
expr_data <- as.data.frame(all_noNorm_withResponse_eset@assayData[["exprs"]])

expr_data_sdy984 <- expr_data[,sdy984$uid[sdy984$study_accession == "SDY984"]]
```

Работа с таблицей по данным экспрессии. Убираем пропущенные значения и пришиваем данные по ответу на вакцинацию.


```r
expr_data_sdy984_trans <- expr_data_sdy984 %>%
  na.omit() %>%
  t()

expression_data <- expr_data_sdy984_trans %>%
  as.data.frame() %>%
  mutate(id = rownames(expr_data_sdy984_trans)) %>%
  inner_join(., sdy984[,c("uid", "MFC_p30", "study_time_collected")], by = c("id" = "uid")) %>%
  select(id, MFC_p30, study_time_collected, everything()) %>%
  filter(!MFC_p30 == "moderateResponder")
```


Проводим T-тест для данных на 0 день. Сравниваем между собой high и low респондеров, чтобы обнаружить гены, которые статистически значимо различаются


```r
expression_data_0_day <- expression_data %>%
  filter(study_time_collected == 0)

t_test_results_0day <- map2_dbl(expression_data_0_day[,-c(1:3)], expression_data_0_day["MFC_p30"], ~wilcox.test(.x ~ .y)$p.value)

adj_ttest_0_day <- as.data.frame(p.adjust(t_test_results_0day, "BH"))

genes_0_day_stat <- adj_ttest_0_day %>%
  rename(p_value = `p.adjust(t_test_results_0day, "BH")`) %>%
  filter(p_value < 0.05)

genes_0_day_stat
```

```
## [1] p_value
## <0 rows> (or 0-length row.names)
```

T-test для всех генов 1 дня


```r
expression_data_1_day <- expression_data %>%
  filter(study_time_collected == 1)

t_test_results_1day <- map2_dbl(expression_data_1_day[,-c(1:3)], expression_data_1_day["MFC_p30"], ~wilcox.test(.x ~ .y)$p.value)

adj_ttest_1_day <- as.data.frame(p.adjust(t_test_results_1day, "BH"))

genes_1_day_stat <- adj_ttest_1_day %>%
  rename(p_value = `p.adjust(t_test_results_1day, "BH")`) %>%
  filter(p_value < 0.05)

genes_1_day_stat
```

```
## [1] p_value
## <0 rows> (or 0-length row.names)
```


T-test для всех генов 3 дня


```r
expression_data_3_day <- expression_data %>%
  filter(study_time_collected == 3)

t_test_results_3day <- map2_dbl(expression_data_3_day[,-c(1:3)], expression_data_3_day["MFC_p30"], ~wilcox.test(.x ~ .y)$p.value)

adj_ttest_3_day <- as.data.frame(p.adjust(t_test_results_3day, "BH"))

genes_3_day_stat <- adj_ttest_3_day %>%
  rename(p_value = `p.adjust(t_test_results_3day, "BH")`) %>%
  filter(p_value < 0.05)

genes_3_day_stat
```

```
## [1] p_value
## <0 rows> (or 0-length row.names)
```


T-test для всех генов 7 дня


```r
expression_data_7_day <- expression_data %>%
  filter(study_time_collected == 7)

t_test_results_7day <- map2_dbl(expression_data_7_day[,-c(1:3)], expression_data_7_day["MFC_p30"], ~wilcox.test(.x ~ .y)$p.value)

adj_ttest_7_day <- as.data.frame(p.adjust(t_test_results_7day, "BH"))

genes_7_day_stat <- adj_ttest_7_day %>%
  rename(p_value = `p.adjust(t_test_results_7day, "BH")`) %>%
  filter(p_value < 0.05)

genes_7_day_stat
```

```
## [1] p_value
## <0 rows> (or 0-length row.names)
```


