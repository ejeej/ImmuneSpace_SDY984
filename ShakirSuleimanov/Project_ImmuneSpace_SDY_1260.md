---
title: "Project_ImmuneSpace_SDY1260"
author: "Shakir Suleimanov"
date: "2022-11-27"
output:
  html_document:
    code_folding: hide
    toc: yes
    toc_float: true
    toc_depth: 3
    theme: cerulean
    keep_md: true
#editor_options:
#  chunk_output_type: console
#bibliography: biblio.bib
#link-citations: yes
#runtime: shiny
---

<style type="text/css">

.math {
font-size: small;
}

</style>
---

  

```r
library(tidyverse)
library(labelled)
library(gtsummary)
library(ggbeeswarm)
library(ggrepel)
library(lemon)
library(ComplexUpset)
library(colorspace)
library(ggpubr)
library(lme4)
library(lmerTest)
library(marginaleffects)
library(broom.mixed)
library(qvalue)
library(knitr)
library(shiny)


options(knitr.kable.NA = '')
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, error = FALSE)
list("style_number-arg:big.mark" = "") %>% set_gtsummary_theme()
```


```r
# Function for upset plots
upset_plot <- function(data, sets_names, rowcolors, title, subtitle, showh = TRUE, ab_hjust = -0.95, nletters = 1, show_num = TRUE) {
  
  p1 <- upset(
    data, sets_names, name = "gene", sort_sets = FALSE,
    sort_intersections = "ascending", sort_intersections_by = c("degree", "cardinality"),
    wrap = TRUE, set_sizes = FALSE,
    base_annotations = list("Inclusive intersection size" = (
      intersection_size(
        mode = "intersect",
        # mapping = aes(fill = exclusive_intersection),
        text_mapping = aes(y = !!get_size_mode("inclusive_intersection"),
                           colour = "on_background"),
        text = list(size = 3.5)
      ) +
        scale_y_continuous(expand = c(0.1,0)) +
        labs(title = sprintf("%s. Incl. intersection size", LETTERS[nletters*2 - 1])) +
        theme(legend.position = "none",
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.ticks.y = element_line(color = "grey50", linewidth = 0.2),
              axis.text.y = element_text(size = 8),
              axis.title.y = element_blank(),
              plot.title = element_text(hjust = ab_hjust, size = 11, face = "bold")))),
    matrix = intersection_matrix(
      outline_color = list(active = "#4F6980", inactive = "white"),
      segment = geom_segment(color = "#4F6980")
    ) +
      scale_color_manual(values = c("TRUE" = "#4F6980", "FALSE" = "white"), guide = "none"),
    stripes = upset_stripes(colors = rowcolors),
    themes = upset_modify_themes(list(
      intersections_matrix = theme(panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   axis.text.y = element_text(size = 10),
                                   axis.title.x = element_blank()))))
  p2 <- upset(
    data, sets_names, name = "gene", sort_sets = FALSE, 
    sort_intersections = "ascending", sort_intersections_by = c("degree", "cardinality"),
    wrap = TRUE, set_sizes = FALSE,
    base_annotations = list("Intersection size" = (
      intersection_size(
        text_mapping = aes(y = !!get_size_mode("exclusive_intersection"),
                           colour = "on_background"),
        text = list(size = 3.5)
      ) + 
        annotate(geom = "text", x = Inf, y = Inf,
                 label = ifelse(show_num, sprintf("%d\nunique genes", nrow(data)), ""),
                 vjust = 1, hjust = 1) +
        scale_y_continuous(expand = c(0.1,0)) +
        labs(title = sprintf("%s. Excl. intersection size", LETTERS[nletters*2])) +
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.ticks.y = element_line(color = "grey50", linewidth = 0.2),
              axis.text.y = element_text(size = 8),
              axis.title.y = element_blank(),
              plot.title = element_text(hjust = ab_hjust, size = 11, face = "bold")))),
    matrix = intersection_matrix(
      outline_color = list(active = "#4F6980", inactive = "white"),
      segment = geom_segment(color = "#4F6980")
    ) +
      scale_color_manual(values = c("TRUE" = "#4F6980", "FALSE" = "white"), guide = "none"),
    stripes = upset_stripes(colors = rowcolors),
    themes = upset_modify_themes(list(
      intersections_matrix = theme(panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(),
                                   axis.text.y = element_text(size = 10),
                                   axis.title.x = element_blank()))))
  
  if (showh) {
    p1 +
      labs(title = title, subtitle = subtitle) +
      theme(plot.title = element_text(size = 13, face = "bold"),
            plot.subtitle = element_text(size = 10)) +
      p2
  } else {
    ((p1 +
        labs(title = title, subtitle = subtitle) +
        theme(plot.title = element_text(size = 13, face = "bold"),
              plot.subtitle = element_text(size = 10))) 
     /
       (p2))
  }
}
# Function for getting BPs (GO) for genes strored in "gene" column of df
bp_genesF <- function(df) {
  AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db, 
    df %>% pull(gene) %>% unique(), 
    "GO", "SYMBOL") %>% 
    filter(ONTOLOGY == "BP") %>%
    mutate(TERM = AnnotationDbi::Term(GO), DEFINITION = AnnotationDbi::Definition(GO)) %>%
    filter(!is.na(TERM) & TERM != "biological_process") %>%
    transmute(SYMBOL, GO, TERM, DEFINITION) %>%
    unique() %>%
    group_by(TERM) %>%
    summarise(GOBPID = unique(GO), def = unique(DEFINITION), n = n(), genes = paste(SYMBOL, collapse = ", ")) %>%
    arrange(-n)
}
# Function to add an image with link
image_link <- function(image, url,...){
  htmltools::a(
    href = url,
    htmltools::img(src = image,...)
  )
}
```


```r
# Data https://www.immunespace.org/project/home/Integrative_Public_Study/begin.view?SDY=IS2
# all_noNorm_withResponse_eset <- readRDS("data/all_noNorm_withResponse_eset.Rds")
all_noNorm_withResponse_eset <- readRDS("all_noNorm_withResponse_eset.Rds")
# Data for participants of the study SDY984 -------------------------------------
df_subj <- all_noNorm_withResponse_eset@phenoData@data %>%
  filter(study_accession == "SDY1260")
timepoints <- unique(df_subj$study_time_collected)
# Baseline data (0 days) for descriptive statistics
# excluding moderate responders (MFC_p40)
df_subj_baseline <- df_subj %>%
  filter(study_time_collected == 0 & MFC_p40 != "moderateResponder") %>% 
  transmute(participant_id,
            arm_accession = factor(arm_accession), 
            age = age_imputed,
            gender = factor(gender), 
            race = fct_infreq(factor(race)),
            ethnicity = factor(ethnicity),
            igg_baseline = ImmResp_baseline_value_MFC,
            response = factor(MFC_p40, c("lowResponder", "highResponder"), 
                              c("Low Responder", "High Responder")))
var_label(df_subj_baseline) <-
  list(arm_accession = "Study arm",
       age = "Age, yrs",
       gender = "Gender",
       race = "Race",
       ethnicity = "Ethnicity",
       igg_baseline = "log2(IgG), ELISA",
       response = "Vaccination response (MFC_p40)")
# Data on gene expression, SDY984 ----------------
expr_t <- all_noNorm_withResponse_eset@assayData$exprs
expr_t <- expr_t[,df_subj$uid]
# Expression matrix --> dataframe
# excluding moderate responders (MFC_p40)
# + participants' characteristics
df_expr <- expr_t %>%
  t() %>%
  as_tibble(rownames = "uid") %>%
  separate(uid, c("participant_id", "time", NA, NA), 
           sep = "_", remove = FALSE, convert = TRUE) %>%
  filter(participant_id %in% df_subj_baseline$participant_id) %>%
  left_join(df_subj %>% transmute(uid, igg_postvax = ImmResp_postVax_value_MFC), 
            by = "uid") %>%
  left_join(df_subj_baseline %>% select(participant_id, arm_accession, gender, 
                                        race, ethnicity, igg_baseline, response), 
            by = "participant_id")
# Exclude genes with missing expression
df_expr <- df_expr %>%
  select(function(x) sum(is.na(x)) != length(x)) %>%
  select(participant_id, uid, arm_accession, gender, race, ethnicity,
         igg_baseline, time, response, igg_postvax, everything())
# Long dataframe
df_expr_long <- df_expr %>%
  pivot_longer(A1BG:last_col(), names_to = "gene", values_to = "expr")
# All offsprings for GO BP = immune system process | defense response | 
# cytokine production | response to cytokine | MHC
gobp_immresp <- tibble(
  GOBPID = c("GO:0002376", GO.db::GOBPOFFSPRING[["GO:0002376"]], 
             "GO:0006952", GO.db::GOBPOFFSPRING[["GO:0006952"]],
             "GO:0001816", GO.db::GOBPOFFSPRING[["GO:0001816"]],
             "GO:0034097", GO.db::GOBPOFFSPRING[["GO:0034097"]],
             "GO:0002396", GO.db::GOBPOFFSPRING[["GO:0002396"]]),
  TERM = AnnotationDbi::Term(GOBPID)) %>%
  unique()
```


## **Descriptive statistics: participants' characteristics at the Baseline**

<br>
  
Descriptive statistics for participants' characteristics at the baseline (before vaccination) are presented in Table 1. As we are going to compare gene expression between participants with strong and weak response to MSPV4 and MCV4 (meningococcu) vaccine in our study, we showed statistics for these two groups in Table 1. All volunteers were present at age 30. Unfortunately, a lot of data was not specified by the researchers.


```r
tbl_summary(
  df_subj_baseline %>% select(-participant_id), 
  by = "response",
  type = all_continuous() ~ "continuous2",
  statistic = list(
    all_continuous() ~ c("{mean} ({sd})", 
                         "{median} ({p25}-{p75})", "{min}-{max}")),
  digits = list(igg_baseline ~ rep(1,7)),
  missing_text = "Н.Д.") %>%
  add_stat_label(label = list(
    all_continuous() ~ c("Mean (SD)",
                         "Median (Q1-Q3)", "Range"))) %>%
  add_p(pvalue_fun = function(x) style_pvalue(x, digits = 3)) %>%
  modify_header(all_stat_cols() ~ "**{level}**<br>N = {n}") %>%
  modify_footnote(p.value ~ "p-value: Study arm, Gender, Race - Pearson's Chi-squared test;<br> 
  Ethnicity - Fisher's exact test; IgG - Mann-Whitney test") %>%
  bold_labels() %>%
  as_kable_extra(caption = "<b>Table 1. Baseline characteristics of the study participants.</b>",
                 addtl_fmt = FALSE) %>% 
  kableExtra::row_spec(0, bold = TRUE) %>%
  kableExtra::kable_classic(full_width = FALSE, position = "left", font_size = 14,
                            html_font = "\"Source Sans Pro\", helvetica, sans-serif")
```

<table style='NAborder-bottom: 0; font-size: 14px; font-family: "Source Sans Pro", helvetica, sans-serif; width: auto !important; ' class=" lightable-classic">
<caption style="font-size: initial !important;"><b>Table 1. Baseline characteristics of the study participants.</b></caption>
 <thead>
  <tr>
   <th style="text-align:left;font-weight: bold;"> Characteristic </th>
   <th style="text-align:center;font-weight: bold;"> Low Responder<br>N = 12 </th>
   <th style="text-align:center;font-weight: bold;"> High Responder<br>N = 12 </th>
   <th style="text-align:center;font-weight: bold;"> p-value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">Study arm, n (%)</span> </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> &gt;0.999 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> ARM4355 </td>
   <td style="text-align:center;"> 7 (58%) </td>
   <td style="text-align:center;"> 7 (58%) </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> ARM4356 </td>
   <td style="text-align:center;"> 5 (42%) </td>
   <td style="text-align:center;"> 5 (42%) </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">Age, yrs, n (%)</span> </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> 30 </td>
   <td style="text-align:center;"> 12 (100%) </td>
   <td style="text-align:center;"> 12 (100%) </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">Gender, n (%)</span> </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Not Specified </td>
   <td style="text-align:center;"> 12 (100%) </td>
   <td style="text-align:center;"> 12 (100%) </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">Race, n (%)</span> </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Not Specified </td>
   <td style="text-align:center;"> 12 (100%) </td>
   <td style="text-align:center;"> 12 (100%) </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">Ethnicity, n (%)</span> </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Not Specified </td>
   <td style="text-align:center;"> 12 (100%) </td>
   <td style="text-align:center;"> 12 (100%) </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">log2(IgG), ELISA</span> </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 0.024 </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Mean (SD) </td>
   <td style="text-align:center;"> 2.5 (1.8) </td>
   <td style="text-align:center;"> 1.1 (1.4) </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Median (Q1-Q3) </td>
   <td style="text-align:center;"> 2.3 (1.8-3.6) </td>
   <td style="text-align:center;"> 1.0 (0.2-1.9) </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> Range </td>
   <td style="text-align:center;"> -0.6-5.8 </td>
   <td style="text-align:center;"> -1.3-4.2 </td>
   <td style="text-align:center;">  </td>
  </tr>
</tbody>
<tfoot><tr><td style="padding: 0; " colspan="100%">
<sup>1</sup> p-value: Study arm, Gender, Race - Pearson's Chi-squared test;<br> 
Ethnicity - Fisher's exact test; IgG - Mann-Whitney test</td></tr></tfoot>
</table>

<br>

## **Descriptive statistics: gene expression**

<br>

The initial dataset contained data in 26925 genes' expression. For the purposes of our research we excluded those genes which had missings on their expression for all participants - it resulted in the data for 16146 genes. We estimate median absolute deviation (MAD) of the expression for these genes (without differentiation by time points) and leave 5 thousand genes with the maximum MAD (genes with the largest variation in expression).


```r
genes_maxvar <- df_expr_long %>%
  group_by(gene) %>%
  summarise(mad_expr = mad(expr)) %>%
  arrange(-mad_expr) %>%
  slice_head(n = 5000)
df_expr_long_fin <- df_expr_long %>%
  filter(gene %in% genes_maxvar$gene) %>%
  mutate(timeF = factor(time, timepoints, ifelse(timepoints == 0, "Baseline", paste(timepoints, "d."))))
```

<br>
  
Below histograms for expression are given for 10 randomly selected genes. 


```r
df_expr_long_fin %>% 
  filter(gene %in% sample(genes_maxvar$gene, size = 10, replace = FALSE)) %>%
  mutate(gene = fct_reorder(gene, expr, max)) %>% 
  ggplot(aes(x = expr, fill = response)) +
  geom_histogram(alpha = 0.6, position = "identity") +
  scale_y_continuous(expand = c(0,0,0.1,0)) +
  scale_x_continuous(breaks = seq(0, max(df_expr_long_fin$expr), 2)) +
  scale_fill_manual(values = c("#E15759", "#4E79A7")) +
  facet_grid(gene ~ timeF, scales = "free", switch = "y") +
  labs(x = "Gene expression", y = element_blank(), fill = element_blank(),
       title = "Figure 1. Distribution of the random 10 genes expressions by time",
       caption = "Counts on the Y-axis (the same scale for all plots)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_text(size = 10, color = "black", face = "bold", angle = 0, hjust = 1),
        strip.text.x = element_text(size = 11, color = "black", face = "bold"),
        plot.title = element_text(size = 13, face = "bold"),
        plot.caption = element_text(size = 10, color = "black", face = "italic", hjust = 0))
```

![](Project_ImmuneSpace_SDY_1260_files/figure-html/fig1_hists-1.png)<!-- -->

<br>
  
## **DEGs: pairwise Mann-Whitney tests**
  
<br>
  
At every time point for every gene we compared its expression between low and high responders using the Mann-Whitney test, correct the obtained p-values with Benjamini-Hochberg method (for the point-wise control of the FDR) and count the number of genes with unadjusted and adjusted p-value less than 0.05. The results are given in Table 2.


```r
mwtest_p <- suppressWarnings({df_expr_long_fin %>%
    group_by(time, timeF, gene) %>%
    summarise(p = wilcox.test(expr ~ response)$p.value) %>%
    group_by(time) %>%
    mutate(p_adj = p.adjust(p, method = "BH")) %>%
    ungroup() %>%
    arrange(p_adj) %>%
    mutate(p_group = cut(p, c(0, 0.05, 1.1), c("< 0.05", "$\\geq$ 0.05"), right = FALSE),
           p_adj_group = cut(p_adj, c(0, 0.05, 1.1), c("< 0.05", "$\\geq$ 0.05"), right = FALSE))
})
genes_sig_mw <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, mwtest_p$gene[mwtest_p$p_adj < 0.05], 
                                      c("GENENAME", "GO"), "SYMBOL") %>% 
  filter(ONTOLOGY == "BP") %>%
  mutate(TERM = AnnotationDbi::Term(GO), DEFINITION = AnnotationDbi::Definition(GO)) %>%
  dplyr::select(SYMBOL, GENENAME, TERM, DEFINITION) %>%
  unique()
tbl_summary(
  mwtest_p %>% select(timeF, p_group, p_adj_group), 
  by = "timeF",
  label = list(p_group = "P-value",
               p_adj_group = "Adjusted p-value")) %>%
  modify_footnote(everything() ~ NA) %>%
  modify_header(label ~ "", 
                all_stat_cols() ~ "**{level}**") %>%
  bold_labels() %>%
  as_kable_extra(caption = "<b>Table 2. P-values of the Mann-Whitney tests before and after<br>
                 adjustment under Benjamini & Hochberg method</b>,<br>number of genes (%)") %>% 
  kableExtra::row_spec(0, bold = TRUE) %>%
  kableExtra::kable_classic(full_width = FALSE, position = "left", font_size = 14,
                            html_font = "\"Source Sans Pro\", helvetica, sans-serif")
```

<table class=" lightable-classic" style='font-size: 14px; font-family: "Source Sans Pro", helvetica, sans-serif; width: auto !important; '>
<caption style="font-size: initial !important;">
<b>Table 2. P-values of the Mann-Whitney tests before and after<br>
                 adjustment under Benjamini &amp; Hochberg method</b>,<br>number of genes (%)</caption>
 <thead>
  <tr>
   <th style="text-align:left;font-weight: bold;">  </th>
   <th style="text-align:center;font-weight: bold;"> Baseline </th>
   <th style="text-align:center;font-weight: bold;"> 3 d. </th>
   <th style="text-align:center;font-weight: bold;"> 7 d. </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">P-value</span> </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> &lt; 0.05 </td>
   <td style="text-align:center;"> 242 (4.8%) </td>
   <td style="text-align:center;"> 375 (7.5%) </td>
   <td style="text-align:center;"> 181 (3.6%) </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> $\geq$ 0.05 </td>
   <td style="text-align:center;"> 4758 (95%) </td>
   <td style="text-align:center;"> 4625 (92%) </td>
   <td style="text-align:center;"> 4819 (96%) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> <span style=" font-weight: bold;    ">Adjusted p-value</span> </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> &lt; 0.05 </td>
   <td style="text-align:center;"> 0 (0%) </td>
   <td style="text-align:center;"> 0 (0%) </td>
   <td style="text-align:center;"> 0 (0%) </td>
  </tr>
  <tr>
   <td style="text-align:left;padding-left: 2em;" indentlevel="1"> $\geq$ 0.05 </td>
   <td style="text-align:center;"> 5000 (100%) </td>
   <td style="text-align:center;"> 5000 (100%) </td>
   <td style="text-align:center;"> 5000 (100%) </td>
  </tr>
</tbody>
</table>
Unfortunately, mostly because of lack of the group size, we did not observed any statistical difference at gene expression between low and high responders. 

<br>
  
  ## **DEGs: linear mixed models**
  
  <br>
  

```r
# Function for quantiles of the empirical distribution function --> ranks
rank_quant <- function(x) {
  ecdf_x <- ecdf(x)
  ecdf_inv <- function(v) {quantile(x, v)}
  qq <- ecdf_inv(ecdf_x(x))
  qq <- names(qq) %>% gsub("%", "", .) %>% as.numeric(.)/100
  qq
}
# "Quantile" ranks of genes by expression for every participant at each time point
df_expr_long_fin <- df_expr_long_fin %>%
  group_by(participant_id, time) %>%
  mutate(expr_rank = rank_quant(expr)) %>%
  ungroup()
# rm(all_noNorm_withResponse_eset, df_expr_long, df_expr, expr_t)
```

### **Model description**

<br>
  
  For every gene out of 5 thousand genes with the largest MAD we will estimate **linear mixed effects model** of the following kind:
  
  $expr_{it} = \beta_0 + \beta_{0i} + \beta_1*time_{j} + \beta_2*response_{i} + \beta_3*time_j * response_i + \epsilon_{ij}$, where

- $expr_{ij}$ - _i_-th individual's gene expression at the _j_-th time point; we will estimate separate specifications for expression measured as:

  - the initially given data, 

  - ranks of genes by expression, calculated for each participant at each time point as the quantile of the empirical distribution function for all genes' expression for this participant at this time point, 

- $time_j$ - time from vaccination (days), 

- $response_i$ - _i_-th participant's response to Zostavax vaccine measured at 30 days from vaccination (1 - high responder, 0 - low responder), 

- $\beta_0$ - regression intercept (global average expression for this particular gene among all participants and time points),

- $\beta_{0i}$ - individual random effect (modeled for each participant as normally distributed variable with zero mean and variance which characterize the between-subject variation of the mean gene expression around the global average gene expression),

- $\beta_{1,2,3}$ - regression coefficients (fixed effects for time, response and their interaction), 

- $\epsilon_{ij}$ - residuals (it is assumed that they are normally distributed with zero mean and variance which characterize the within-subject variation of the gene expression (in time)).

This model allows to estimate whether the response effect changes significantly over time (or, equivalently, whether the change in expression differs significantly depending on the response to vaccination). $\beta_2$ estimate in this model can be interpreted as the ATE of the high response at the Baseline ($time=0$), i.e. the baseline difference in the average gene expression between high and low responders. For other time points ATE can be calculated as $\beta_2 + \beta_3*time$. Therefore, $\beta_3$ estimate can be interpreted as the average change in the ATE of the high response with every additional day after vaccination.

<br>

Linear mixed models will be estimated with the `lmer` function out of the `lme4` package [@bates_fitting_2015] for R (version 4.1.1), average marginal effects will be estimated with the means of the `marginaleffects` function out of the package with the same name [@arel-bundock_marginaleffects_2023]. It can be noted that if we give particular time points to the latter function after the mixed model with the interaction term (e.g., $time=$ 0,3,7), then near the $response$ variable we will get the desired values for the ATE of the high response at these points ($\beta_2 + \beta_3*time$) with the corresponding p-values. If to not pass the time points into this function, near the $response$ variable we will get the estimate for the ATE of the high response from the model without the interaction term.


```r
# Linear mixed models
lmer_q <- purrr::quietly(.f = lmer)
lmer_fits_SDY1260 <- df_expr_long_fin %>%
  nest(-gene) %>%
  mutate(
    fit_noint = map(data, ~ lmer_q(expr ~ time + response + (1|participant_id), .x)),
    fit_int = map(data, ~ lmer_q(expr ~ time*response + (1|participant_id), .x)),
    fit_nointrank = map(data, ~ lmer_q(expr_rank ~ time + response + (1|participant_id), .x)),
    fit_intrank = map(data, ~ lmer_q(expr_rank ~ time*response + (1|participant_id), .x)))
lmer_res_SDY1260 <- lmer_fits_SDY1260 %>%
  mutate(
    beta_p_noint = map(fit_noint, ~ .x$result %>% broom.mixed::tidy()),
    beta_p_nointrank = map(fit_nointrank, ~ .x$result %>% broom.mixed::tidy()),
    beta_p_int = map(fit_int, ~ .x$result %>% broom.mixed::tidy()),
    beta_p_intrank = map(fit_intrank, ~ .x$result %>% broom.mixed::tidy()),
    ame_int = map(
      fit_int, ~ .x$result %>% 
        marginaleffects(newdata = datagrid(time = timepoints, response = levels(df_subj_baseline$response)), eps = 0.001)),
    ame_intrank = map(
      fit_intrank, ~ .x$result %>% 
        marginaleffects(newdata = datagrid(time = timepoints, response = levels(df_subj_baseline$response)), eps = 0.001))) %>%
  select(-data, -contains("fit"))
# # Save results to rds (without sending to github)
  saveRDS(lmer_res_SDY1260, "lmer_res_SDY1260.rds")
# # lmer_res <- readRDS("OlgaMironenko/res/lmer_res.rds")
# # lmer_res <- readRDS(file.path("..", "OlgaMironenko", "res", "lmer_res.rds"))
# Betas and p-values for models with interaction term
lmer_betas_int_SDY1260 <- full_join(
  lmer_res_SDY1260 %>%
    select(gene, beta_p_int) %>%
    unnest(beta_p_int) %>%
    select(gene, term, estimate, p.value),
  lmer_res_SDY1260 %>%
    select(gene, beta_p_intrank) %>%
    unnest(beta_p_intrank) %>%
    select(gene, term, estimate, p.value),
  by = c("gene", "term"), suffix = c("_init", "_rank"))
# AMEs/ATEs and p-values for models with interaction term
lmer_ames_int_SDY1260 <- full_join(
  lmer_res_SDY1260 %>%
    select(gene, ame_int) %>%
    unnest(ame_int) %>%
    select(gene, term, time, response, dydx, p.value),
  lmer_res_SDY1260 %>%
    select(gene, ame_intrank) %>%
    unnest(ame_intrank) %>%
    select(gene, term, time, response, dydx, p.value),
  by = c("gene", "term", "time", "response"), suffix = c("_init", "_rank"))
# Betas and p-values for models without interaction term
lmer_betas_noint_SDY1260 <- full_join(
  lmer_res_SDY1260 %>%
    select(gene, beta_p_noint) %>%
    unnest(beta_p_noint) %>%
    select(gene, term, estimate, p.value),
  lmer_res_SDY1260 %>%
    select(gene, beta_p_nointrank) %>%
    unnest(beta_p_nointrank) %>%
    select(gene, term, estimate, p.value),
  by = c("gene", "term"), suffix = c("_init", "_rank"))
# # Saving results to rds
saveRDS(lmer_betas_int_SDY1260, "lmer_betas_int_SDY1260.rds")
saveRDS(lmer_betas_noint_SDY1260, "lmer_betas_noint_SDY1260.rds")
saveRDS(lmer_ames_int_SDY1260, "lmer_ames_int_SDY1260.rds")
```


```r
# lmer_betas_int <- readRDS("OlgaMironenko/res/lmer_betas_int.rds")
# lmer_betas_noint <- readRDS("OlgaMironenko/res/lmer_betas_noint.rds")
# lmer_ames_int <- readRDS("OlgaMironenko/res/lmer_ames_int.rds")
lmer_betas_int_SDY1260 <- readRDS("lmer_betas_int_SDY1260.rds")
lmer_betas_noint_SDY1260 <- readRDS("lmer_betas_noint_SDY1260.rds")
lmer_ames_int_SDY1260 <- readRDS("lmer_ames_int_SDY1260.rds")
critv <- -log10(0.05)
# Data frame for all betas and p-values for all LMMs
lmer_res_long <- bind_rows(
  lmer_betas_int_SDY1260 %>%
    filter(grepl("response", term)) %>%
    transmute(gene, model = "int", term = ifelse(grepl("time", term), "b3", "b2"), 
              estimate_init, estimate_rank, p.value_init, p.value_rank,
              time = NA) %>%
    pivot_longer(cols = -c(gene, model, term, time), 
                 names_pattern = "(estimate|p.value)_(.+)$",
                 names_to = c(".value", "expr")),
  lmer_betas_noint_SDY1260 %>%
    filter(grepl("response", term)) %>%
    pivot_longer(cols = -c(gene, term), names_pattern = "(estimate|p.value)_(.+)$",
                 names_to = c(".value", "expr")) %>%
    mutate(model = "noint", term = "b2", time = NA),
  lmer_ames_int_SDY1260 %>%
    filter(term == "response" & response == "Low Responder") %>%
    select(-term, -response) %>%
    pivot_longer(cols = -c(gene, time), names_pattern = "(dydx|p.value)_(.+)$",
                 names_to = c(".value", "expr")) %>%
    mutate(model = "int", term = sprintf("ame%d", time)) %>%
    rename(estimate = dydx)
) %>%
  mutate(logp = -log10(p.value),
         sig = logp > critv,
         model_var = sprintf("%s_%s_%s", expr, model, term),
         expr = factor(expr, c("init", "rank"), c("Expression as is", "Ranks by expression")),
         time = factor(time, timepoints, labels = ifelse(timepoints == 0, "Baseline", paste(timepoints, "d."))))
# Significant genes
genes_sig <- map(
  lmer_res_long %>%
    filter(sig) %>%
    split(f = as.factor(.$model_var)),
  ~ .x %>% pull(gene))
models <- names(genes_sig)
models_num <- as.numeric(str_extract(models, "\\d+"))
genes_sig_lbls <- sprintf("%s, %s (%s)",
                          ifelse(grepl("init",models), "Expr.", "Ranks"),
                          ifelse(grepl("ame",models), "sig. ATE", 
                                 ifelse(models_num == 2, "sig. b2", "sig. b3")),
                          ifelse(grepl("ame0|_int_b2", models), "Baseline",
                                 ifelse(grepl("b3", models), "Change",
                                        ifelse(grepl("noint", models), "Overall",
                                               sprintf("%d d.", models_num))))) %>%
  as.list() %>% setNames(models)
# Simultaneous significance for the main and interaction terms
lmer_sig_int <- lmer_res_long %>%
  filter(model == "int" & !grepl("ame", term)) %>%
  pivot_wider(id_cols = c(gene, expr), names_from = "term", values_from = c("logp", "sig"))
genes_sig$init_int_b2b3 <- lmer_sig_int %>% filter(sig_b2 & sig_b3 & expr == "Expression as is") %>% pull(gene)
genes_sig$rank_int_b2b3 <- lmer_sig_int %>% filter(sig_b2 & sig_b3 & expr != "Expression as is") %>% pull(gene)
genes_sig_lbls$init_int_b2b3 <- "Expr., sig. b2 and b3"
genes_sig_lbls$rank_int_b2b3 <- "Ranks, sig. b2 and b3"
# Whether each gene was significant in any model
genes_maxvar <- cbind(
  genes_maxvar,
  map_dfc(genes_sig, ~ genes_maxvar$gene %in% .x))
```

<br>

### **Results for linear mixed models**

<br>

#### **Differentially expressed gene sets**

<br>

Upon estimation of the linear mixed models we got 65 **genes with significant coefficients** for both the main response effect and its interaction with time in the specifications with the initial data for expression (71 in the specification with the ranks) - these genes were not only differentially expressed at the Baseline, but differed in the dynamics of the mean expression after vaccination in dependence of the response strength. In addition, 291 (284) genes had significant main effect for the response only, i.e. were differentially expressed at the Baseline, but got similar dynamics of the expression among low and high responders. For 138 (141) genes the main effect of the response was insignificant, while the interaction term was significant, i.e. these genes got similar average expression before vaccination, but diverged in its dynamics afterwards.

We can **compare p-values for the main effect of the response and response-time interaction** for each gene using the scatter diagram, showing p-value for the main effect ($\beta_2$ estimate) on the X axis and p-value for the interaction term ($\beta_3$ estimate) on the Y axis. For better clarity, we use $-log_{10}$-transformation for both of p-values. In each quadrant except the lower left (here are genes with both insignificant effects), we label 5 genes with the lowest p-values.


```r
lmer_p_int_plt <- lmer_sig_int %>%
  group_by(expr, sig_b2, sig_b3) %>%
  arrange(-pmax(logp_b2, logp_b3)) %>%
  mutate(sig_plot = row_number() %in% c(1:5)) %>%
  ungroup() %>%
  mutate(sig_plot = ifelse(!sig_b2 & !sig_b3, FALSE, sig_plot))
ggplot(lmer_p_int_plt, aes(x = logp_b2, y = logp_b3)) +
  geom_point(aes(color = sig_plot), alpha = 0.5, show.legend = FALSE) +
  geom_hline(aes(yintercept = 1), linewidth = 0.7, color = "#E15759", linetype = "dashed") +
  geom_vline(aes(xintercept = 1), linewidth = 0.7, color = "#E15759", linetype = "dashed") +
  geom_hline(aes(yintercept = critv), linewidth = 0.7, color = "#E15759") +
  geom_vline(aes(xintercept = critv), linewidth = 0.7, color = "#E15759") +
  geom_text_repel(aes(label = gene), lmer_p_int_plt %>% filter(sig_plot), size = 3) +
  scale_x_continuous(expand = c(0.02, 0)) +
  scale_y_continuous(expand = c(0.02, 0)) +
  scale_color_manual(values = c("#86BCB6", "#F28E2B")) +
  facet_rep_wrap(~ expr, nrow = 2, ncol = 2, repeat.tick.labels = TRUE) +
  labs(x = bquote(-log[10](p)~", response"), y = bquote(-log[10](p)~", response*time"), 
       title = "Figure 3. P-values for regression coefficients estimates",
       subtitle = "Linear mixed models with the interaction term",
       caption = "Solid red line is for p = 0.05, dashed red line is for p = 0.1. The higher the axis value, the lower is the p-value.") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        # axis.line = element_blank(),
        panel.grid.major = element_line(linewidth = .2, color = '#ebebebFF'),
panel.grid.minor = element_line(linewidth = .1, color = '#ebebebFF'),
strip.background = element_blank(),
strip.text.x = element_text(size = 12, color = "black", face = "bold"),
plot.title = element_text(size = 13, face = "bold"),
plot.subtitle = element_text(size = 10),
plot.caption = element_text(size = 10, color = "black", face = "italic", hjust = 0))
```

![](Project_ImmuneSpace_SDY_1260_files/figure-html/fig3_p12-1.png)<!-- -->

In addition, we can see **how much significant gene sets, revealed by the significance of the coefficients near response and its interaction with time, are intersected** (Figure 4).


```r
upset_plt_df <- cbind(
  gene = genes_maxvar$gene,
  map_dfc(genes_sig[grepl("_int_b\\d$", names(genes_sig))], ~ genes_maxvar$gene %in% .x)) %>%
  filter(if_any(-gene, ~.)) %>%
  setNames(c("gene", genes_sig_lbls[names(.)[-1]]))
upset_plot(upset_plt_df, rev(colnames(upset_plt_df)[-1]), 
           c(rep("#A0CBE8", 2), rep("#F1CE63", 2)),
           "Figure 4. Intersection of gene sets with significant coefficients",
           "Linear mixed models with the interaction term", FALSE, -0.3) +
  labs(caption = "sig. b2 (Baseline) = p < 0.05 for the coeff-t near the main effect of the response,\n sig. b3 (Change) = p < 0.05 for the coeff-t near the response * time interaction.") +
  theme(plot.caption = element_text(size = 10, color = "black", face = "italic", hjust = 0))
```

![](Project_ImmuneSpace_SDY_1260_files/figure-html/fig4_upset_int_noint-1.png)<!-- -->

<br>
  

```r
# Genes with insignificant coefficient at baseline, but significant at change
genes_sig_b3 <- genes_maxvar %>% 
  filter(!init_int_b2 & !rank_int_b2 & (init_int_b3 | rank_int_b3))
# Genes with significant coefficient at baseline, but insignificant at change
genes_sig_b2 <- genes_maxvar %>% 
  filter((init_int_b2 | rank_int_b2) & !init_int_b3 & !rank_int_b3)
```

Combining the results for the both specifications of gene expression we found **343 unique genes with the significant main effect for response (in any of specifications) and its insignificant interaction with time (in both specifications)**, and **160 unique genes with significant interaction (in any specification) and insignificant main term (in both specifications)**. We can say that the average expression for the former genes was significantly different between high and low responders before vaccination, but its dynamics was similar for these groups after vaccination. On the contrary, the latter genes got similar average expression at the Baseline for both groups, but differed in its change over time.

Using information from Gene Ontology for each of these two gene sets we extract those **signalling paths (biological processes)** into which the products of these genes are involved. Figure 5 (1) represents the revealed processes with the largest number of genes. Figure 5 (2) lists the most represented immune-related processes (those of them which were missing among ones with the significant timre-response interaction are marked in blue color).


```r
gobp_sig_b2 <- bp_genesF(genes_sig_b2)
gobp_sig_b3 <- bp_genesF(genes_sig_b3)
gobp_sig_b2b3_plot <- bind_rows(
  gobp_sig_b2 %>% 
    mutate(time = "Sig. b2 only\n(Baseline)"),
  gobp_sig_b3 %>% 
    mutate(time = "Sig. b3 only\n(Change)")) %>%
  mutate(immresp = factor(GOBPID %in% gobp_immresp$GOBPID, labels = c("Not related to immune response", "Related to immune response"))) %>%
  arrange(-n) 
g1 <- gobp_sig_b2b3_plot %>% filter(time == "Sig. b2 only\n(Baseline)") %>% head(20) %>% pull(GOBPID)
g2 <- gobp_sig_b2b3_plot %>% filter(time != "Sig. b2 only\n(Baseline)") %>% head(20) %>% pull(GOBPID)
plt_df <- gobp_sig_b2b3_plot %>% 
  filter(GOBPID %in% c(g1,g2)) %>%
  mutate(TERM = fct_reorder(factor(TERM), n))
ggplot() +
  geom_bar(aes(y = TERM, x = n, fill = immresp), plt_df, stat = "identity", width = 0.6) +
  geom_vline(xintercept = 0) +
  facet_grid(~ time, scales = "fixed") +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("#76B7B2", "#F28E2B")) +
  labs(y = element_blank(), fill = element_blank(),
       x = "Number of genes with sig. coefficient",
       title = "Figure 5 (1). Biological processes (GO), represented by the largest number of\ngenes with significant coefficients",
       subtitle = "Coefficients after linear mixed models") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        legend.justification = 1,
        legend.key.size = unit(0.8, "lines"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 11),
        panel.grid.major.x = element_line(linewidth = .2, color = '#ebebebFF'),
        panel.grid.minor.x = element_line(linewidth = .1, color = '#ebebebFF'),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 11, color = "black", face = "bold"),
        plot.title = element_text(size = 13, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(size = 10))
```

![](Project_ImmuneSpace_SDY_1260_files/figure-html/fig5_bps_b2b3-1.png)<!-- -->

```r
g1 <- gobp_sig_b2b3_plot %>% filter(GOBPID %in% gobp_immresp$GOBPID & time == "Sig. b2 only\n(Baseline)") %>% head(15) %>% pull(GOBPID)
g2 <- gobp_sig_b2b3_plot %>% filter(GOBPID %in% gobp_immresp$GOBPID & time != "Sig. b2 only\n(Baseline)") %>% head(15) %>% pull(GOBPID)
plt_df <- gobp_sig_b2b3_plot %>% 
  filter(GOBPID %in% c(g1,g2)) %>%
  mutate(TERM = fct_reorder(factor(TERM), n))
ggplot() +
  geom_bar(aes(y = TERM, x = n), plt_df, fill = "#FFBE7D", stat = "identity", width = 0.6) +
  geom_vline(xintercept = 0) +
  facet_grid(~ time, scales = "fixed") +
  scale_x_continuous(expand = c(0,0)) +
  labs(y = element_blank(), 
       x = "Number of genes with sig. coefficient",
       title = "Figure 5 (2). Immune-related biological processes (GO), represented by the largest number of\ngenes with significant coefficients",
       subtitle = "Coefficients after linear mixed models") +
  theme_classic(base_size = 12) +
  theme(legend.justification = 1,
        legend.key.size = unit(0.8, "lines"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 11),
        axis.text.y = element_text(
          face = ifelse(levels(plt_df$TERM) %in% plt_df$TERM[plt_df$time != "Sig. b2 only\n(Baseline)"],
                        "plain", "bold"),
          color = ifelse(levels(plt_df$TERM) %in% plt_df$TERM[plt_df$time != "Sig. b2 only\n(Baseline)"],
                         "black", "#4E79A7")),
                         panel.grid.major.x = element_line(linewidth = .2, color = '#ebebebFF'),
        panel.grid.minor.x = element_line(linewidth = .1, color = '#ebebebFF'),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 11, color = "black", face = "bold"),
        panel.spacing.x = unit(0.7, "lines"),
        plot.title = element_text(size = 13, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(size = 10))
```

![](Project_ImmuneSpace_SDY_1260_files/figure-html/fig5_bps_b2b3-2.png)<!-- -->

<br>
  
  Besides regression coefficients, we can also estimate **the average treatment effect (ATE) of the response at particular time points**. Figures 6(1) and 6(2) shows the intersection of the genes sets with the significant ATE between time points.


```r
upset_plt_df_1 <- genes_maxvar %>%
  select(contains("init_int_ame")) %>%
  filter(if_any(everything(), ~.)) %>%
  setNames(genes_sig_lbls[names(.)])
upset_plot(upset_plt_df_1, rev(colnames(upset_plt_df_1)), 
           rep("#F1CE63", 4),
           "Figure 6(1). Intersection of gene sets with significant ATEs,\nExpression as is",
           "Linear mixed models with the interaction term", TRUE, -0.8)
```

![](Project_ImmuneSpace_SDY_1260_files/figure-html/fig6_upset_ame-1.png)<!-- -->

```r
upset_plt_df_2 <- genes_maxvar %>%
  select(contains("rank_int_ame")) %>%
  filter(if_any(everything(), ~.)) %>%
  setNames(genes_sig_lbls[names(.)])
upset_plot(upset_plt_df_2, rev(colnames(upset_plt_df_2)), 
           rep("#A0CBE8", 4),
           "Figure 6(2). Intersection of gene sets with significant ATEs,\nRanks by expression",
           "Linear mixed models with the interaction term", TRUE, -0.8)
```

![](Project_ImmuneSpace_SDY_1260_files/figure-html/fig6_upset_ame-2.png)<!-- -->

The noticeably small number of genes have the significant ATE over the whole range of time points when expression was assessed (0-7 days). There were also quite few genes with insignificant baseline ATE of response and significant ATE from 1 to 7 days, as well genes with significant ATE at the Baseline and insignificant subsequent ATEs. On the contrast, there are two modes: **genes for which the ATE became significant only at 7 days after vaccination and genes for which the ATE turned from significant at 0-3 days to insignificant at 7 days**.


```r
# Genes with insignificant ATE at 7d., but significant at 0-3 d.
genes_sig_ame03 <- genes_maxvar %>% 
  filter((!init_int_ame7 & !rank_int_ame7 & init_int_ame0 & init_int_ame3) |
           (!init_int_ame7 & !rank_int_ame7 & rank_int_ame0 & rank_int_ame3))
# Genes with significant ATE at 7d., but insignificant at 0-3 d.
genes_sig_ame7 <- genes_maxvar %>% 
  filter((init_int_ame7 | rank_int_ame7) & 
           !init_int_ame0 & !init_int_ame3 &
           !rank_int_ame0 & !rank_int_ame3)
```

After combining the results for the specifications with the initial data on expression and ranks of genes by expression we obtain 165 genes with significant 7-day ATE (in any specification) and insignificant 0 day ATEs (in both specifications), as well as 253 genes with significant 0-3-day ATEs (in any specification) and insignificant 7-day ATE (in both specifications). The most represented by these two gene sets **biological processes** are shown in Figures 7 (1) and 7(2).


```r
gobp_sig_ame03 <- bp_genesF(genes_sig_ame03)
gobp_sig_ame7 <- bp_genesF(genes_sig_ame7)
gobp_sig_ame_plot <- bind_rows(
  gobp_sig_ame03 %>% 
    mutate(time = "Sig. at 0-3d. only"),
  gobp_sig_ame7 %>% 
    mutate(time = "Sig. at 7d. only")) %>%
  mutate(immresp = factor(GOBPID %in% gobp_immresp$GOBPID, labels = c("Not related to immune response", "Related to immune response"))) %>%
  arrange(-n)
g1 <- gobp_sig_ame_plot %>% filter(time == "Sig. at 0-3d. only") %>% head(17) %>% pull(GOBPID)
g2 <- gobp_sig_ame_plot %>% filter(time != "Sig. at 0-3d. only") %>% head(17) %>% pull(GOBPID)
plt_df <- gobp_sig_ame_plot %>% 
  filter(GOBPID %in% c(g1,g2)) %>%
  mutate(TERM = fct_reorder(factor(TERM), n))
ggplot() +
  geom_bar(aes(y = TERM, x = n, fill = immresp), plt_df, stat = "identity", width = 0.6) +
  geom_vline(xintercept = 0) +
  facet_grid(~ time, scales = "fixed") +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("#76B7B2", "#F28E2B")) +
  labs(y = element_blank(), fill = element_blank(),
       x = "Number of genes with sig. ATE",
       title = "Figure 7 (1). Biological processes (GO), represented by the largest number of\ngenes with significant ATEs",
       subtitle = "ATEs after linear mixed models") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        legend.justification = 1,
        legend.key.size = unit(0.8, "lines"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 10),
        panel.grid.major.x = element_line(linewidth = .2, color = '#ebebebFF'),
        panel.grid.minor.x = element_line(linewidth = .1, color = '#ebebebFF'),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 11, color = "black", face = "bold"),
        plot.title = element_text(size = 13, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(size = 10))
```

![](Project_ImmuneSpace_SDY_1260_files/figure-html/fig7_bps_ame-1.png)<!-- -->

```r
g1 <- gobp_sig_ame_plot %>% filter(GOBPID %in% gobp_immresp$GOBPID & time == "Sig. at 0-3d. only") %>% head(15) %>% pull(GOBPID)
g2 <- gobp_sig_ame_plot %>% filter(GOBPID %in% gobp_immresp$GOBPID & time != "Sig. at 0-3d. only") %>% head(15) %>% pull(GOBPID)
plt_df <- gobp_sig_ame_plot %>% 
  filter(GOBPID %in% c(g1,g2)) %>%
  mutate(TERM = fct_reorder(factor(TERM), n))
ggplot() +
  geom_bar(aes(y = TERM, x = n), plt_df, fill = "#FFBE7D", stat = "identity", width = 0.6) +
  geom_vline(xintercept = 0) +
  facet_grid(~ time, scales = "fixed") +
  scale_x_continuous(expand = c(0,0), breaks = seq(0,20,2)) +
  labs(y = element_blank(), fill = element_blank(),
       x = "Number of genes with sig. ATE",
       title = "Figure 7 (2). Immune-related biological processes (GO), represented by the largest number of\ngenes with significant ATEs",
       subtitle = "ATEs after linear mixed models") +
  theme_classic(base_size = 12) +
  theme(legend.justification = 1,
        legend.key.size = unit(0.8, "lines"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.text.y = element_text(
          face = ifelse(levels(plt_df$TERM) %in% plt_df$TERM[plt_df$time != "Sig. at 0-3d. only"] &
                          levels(plt_df$TERM) %in% plt_df$TERM[plt_df$time == "Sig. at 0-3d. only"],
                        "plain", "bold"),
          color = ifelse(levels(plt_df$TERM) %in% plt_df$TERM[plt_df$time != "Sig. at 0-3d. only"] &
                           levels(plt_df$TERM) %in% plt_df$TERM[plt_df$time == "Sig. at 0-3d. only"],
                         "black", 
                         ifelse(levels(plt_df$TERM) %in% plt_df$TERM[plt_df$time != "Sig. at 0-3d. only"],
                                "#E15759", "#4E79A7"))),
                                panel.grid.major.x = element_line(linewidth = .2, color = '#ebebebFF'),
        panel.grid.minor.x = element_line(linewidth = .1, color = '#ebebebFF'),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 11, color = "black", face = "bold"),
        plot.title = element_text(size = 13, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(size = 10))
```

![](Project_ImmuneSpace_SDY_1260_files/figure-html/fig7_bps_ame-2.png)<!-- -->

<br>
  
#### **Significance vs. size of the average response effect**
  
<br>
  
At the next step in the analysis of the results obtained from the linear mixed models we can **match the statistical significance of the average response effect to its size**. 

Figure 8 contains the volcano plots for the ATEs at particular time points (on the X axis) and $-log_{10}$-transformed p-values for them (on the Y axis). Points for 10 genes with the largest negative and largest positive significant ATEs are highlighted at every point.

The same plots are drawn for **the average response effect at particular time points** estimated after linear mixed models with the response-time interaction (Figure 9).


```r
lmer_v_int_plt <- lmer_res_long %>%
  filter(grepl("ame", term)) %>%
  group_by(expr, time, sig) %>%
  arrange(estimate) %>%
  mutate(sig_plot = row_number() %in% c(1:10, (n()-4):n())) %>%
  ungroup() %>%
  mutate(sig_plot = ifelse(!sig, FALSE, sig_plot))
ggplot(lmer_v_int_plt, aes(x = estimate, y = logp)) +
  geom_point(aes(color = sig_plot), alpha = 0.5, show.legend = FALSE) +
  geom_hline(aes(yintercept = 1), linewidth = 0.7, color = "#E15759", linetype = "dashed") +
  geom_hline(aes(yintercept = critv), linewidth = 0.7, color = "red") +
  geom_text_repel(aes(label = gene), lmer_v_int_plt %>% filter(sig_plot), size = 3) +
  scale_x_continuous(expand = c(0.02, 0)) +
  scale_y_continuous(expand = c(0.02, 0)) +
  scale_color_manual(values = c("#86BCB6", "#F28E2B")) +
  facet_rep_grid(time ~ expr, repeat.tick.labels = TRUE, scales = "free", switch = "y") +
  labs(x = "ATE for response at the time point", y = bquote(-log[10](p[ATE])), 
       title = "Figure 8. P-value vs. effect size for response",
       subtitle = "Average marginal effects after linear mixed models with the interaction term",
       caption = "Solid red line is for p = 0.05, dashed red line is for p = 0.1.\nThe higher the Y axis value, the lower is the p-value.") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.major = element_line(linewidth = .2, color = '#ebebebFF'),
        panel.grid.minor = element_line(linewidth = .1, color = '#ebebebFF'),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, color = "black", face = "bold"),
        strip.placement = "outside",
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10),
        plot.caption = element_text(size = 10, color = "black", face = "italic", hjust = 0))
```

![](Project_ImmuneSpace_SDY_1260_files/figure-html/fig8_volcano_time-1.png)<!-- -->


<br>
  
Let's call genes with ATE > 0.25 (in the intitial measurements units for expressiob) as **upregulated** and genes with ATE < -0.25 as **downregulated** at the corresponding time point. The average expression for every gene in the former gene set is significantly higher among high responders in comparison to the low ones, for the latter it is significantly lower. Using Gene Ontology we obtained annotations for the corresponding **biological processes** for these sets, chose 10 of the most represented of them for each point and showed them at Figure 9.


```r
gobp_time <- bind_rows(
  map_dfr(paste0("ame", timepoints),
          ~ bp_genesF(lmer_v_int_plt %>% filter(term == .x & sig & expr == "Expression as is" & estimate > 0.25)) %>%
            mutate(reg = "+", time = as.numeric(str_extract(.x, "\\d")))),
  map_dfr(paste0("ame", timepoints),
          ~ bp_genesF(lmer_v_int_plt %>% filter(term == .x & sig & expr == "Expression as is" & estimate < -0.25)) %>%
            mutate(reg = "-", time = as.numeric(str_extract(.x, "\\d"))))) %>%
  mutate(time = factor(time, timepoints, ifelse(timepoints == 0, "Baseline", paste(timepoints, "d."))))
gobp_time_plot <- gobp_time %>%
  group_by(time, GOBPID) %>%
  summarise(term = unique(TERM), def = unique(def),
            n_up = n[reg == "+"], n_d = n[reg == "-"], n_sum = sum(n)) %>%
  arrange(-n_sum)
g1 <- map_dfr(levels(gobp_time_plot$time), ~ gobp_time_plot %>% filter(time == .x) %>% head(10)) %>%
  pull(GOBPID) %>% unique()
gobp_time_plot <- gobp_time_plot %>%
  filter(GOBPID %in% g1) %>%
  pivot_longer(cols = c(n_up, n_d), names_to = "reg", values_to = "n") %>%
  group_by(GOBPID) %>%
  mutate(n_sum = sum(n)) %>%
  ungroup() %>%
  mutate(reg = factor(reg, c("n_up", "n_d"), c("Upregulated", "Downregulated")),
         term = fct_reorder(factor(str_wrap(term, 40)), n_sum),
         immresp = factor(GOBPID %in% gobp_immresp$GOBPID, labels = c("Not related\nto imm.resp.", "Related to\nimm.resp.")))
ggplot() +
  geom_bar(aes(y = term, x = n, fill = reg), gobp_time_plot, 
           stat = "identity", position = "stack", width = 0.6) +
  geom_hline(aes(yintercept = y), data.frame(y = c(1:20) + 0.5), color = "grey50") +
  geom_vline(xintercept = 0) +
  facet_grid(immresp ~ time, scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("#499894", "#EDC948")) +
  labs(y = element_blank(), fill = element_blank(),
       x = "Number of genes with big significant ATE",
       title = "Figure 9. Biological processes (GO), represented by the largest number of\ngenes with the largest significant ATEs",
       subtitle = "ATEs after linear mixed models") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.8, "lines"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 11),
        panel.grid.major.x = element_line(linewidth = .2, color = '#ebebebFF'),
panel.grid.minor.x = element_line(linewidth = .1, color = '#ebebebFF'),
panel.grid.major.y = element_blank(),
panel.grid.minor.y = element_blank(),
strip.background = element_blank(),
strip.text = element_text(size = 12, color = "black", face = "bold"),
strip.placement = "outside",
panel.spacing = unit(1, "lines"),
plot.title = element_text(size = 13, face = "bold"),
plot.title.position = "plot",
plot.subtitle = element_text(size = 10))
```

![](Project_ImmuneSpace_SDY_1260_files/figure-html/fig9_bps_big-1.png)<!-- -->

<br>
  
In addition, we filtered only immune-related processes and showed the most represented of them at Figure 10.


```r
gobp_imm_time_plot <- gobp_time %>%
  filter(GOBPID %in% gobp_immresp$GOBPID) %>%
  group_by(time, GOBPID) %>%
  summarise(term = unique(TERM), def = unique(def),
            n_up = n[reg == "+"], n_d = n[reg == "-"], n_sum = sum(n)) %>%
  arrange(-n_sum)
g1 <- map_dfr(levels(gobp_imm_time_plot$time), ~ gobp_imm_time_plot %>% filter(time == .x) %>% head(15)) %>%
  pull(GOBPID) %>% unique()
gobp_imm_time_plot <- gobp_imm_time_plot %>%
  filter(GOBPID %in% g1) %>%
  pivot_longer(cols = c(n_up, n_d), names_to = "reg", values_to = "n") %>%
  group_by(GOBPID) %>%
  mutate(n_sum = sum(n)) %>%
  ungroup() %>%
  mutate(reg = factor(reg, c("n_up", "n_d"), c("Upregulated", "Downregulated")),
         term = fct_reorder(factor(str_wrap(term, 50)), n_sum))
ggplot() +
  geom_bar(aes(y = term, x = n, fill = reg), gobp_imm_time_plot, 
           stat = "identity", position = "stack", width = 0.6) +
  geom_hline(aes(yintercept = y), data.frame(y = c(1:length(unique(gobp_imm_time_plot$GOBPID))) + 0.5), color = "grey50") +
  geom_vline(xintercept = 0) +
  facet_grid(~ time, scales = "fixed") +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("#499894", "#EDC948")) +
  labs(y = element_blank(), fill = element_blank(),
       x = "Number of genes with big significant ATE",
       title = "Figure 10. Immune-related biological processes (GO), represented by the largest number of\ngenes with the largest significant ATEs",
       subtitle = "ATEs after linear mixed models") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.8, "lines"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 11),
        panel.grid.major.x = element_line(linewidth = .2, color = '#ebebebFF'),
        panel.grid.minor.x = element_line(linewidth = .1, color = '#ebebebFF'),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, color = "black", face = "bold"),
        panel.spacing = unit(1, "lines"),
        plot.title = element_text(size = 13, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(size = 10))
```

![](Project_ImmuneSpace_SDY_1260_files/figure-html/fig10_bps_big_imm-1.png)<!-- -->

<br>
  
## **Probability of high response vs. gene expression: logistic regression**
  
<br>
  

```r
df_expr_wide_fin <- df_expr_long_fin %>%
  mutate(expr_rank = expr_rank * 100,
         response = as.numeric(response) - 1) %>%
  pivot_wider(id_cols = c(gene, participant_id, response), 
              names_from = "time", values_from = c("expr", "expr_rank"))
```

### **Model description**

<br>
  
  For every gene out of 5 thousand genes with the largest MAD we will estimate **logistic regression** of the following kind:
  
  $log(odds(response_{i})) = \gamma_0 + \gamma_1*expr_{ij}$, where:
  
  - $response_i$ - response to vaccination for the _i_-th participant (1 - high responder, 0 - low responder), 

- $expr_{ij}$ - gene expression for the _i_-th participant at the _j_-th time point, $j=0,1,3,7$. As for linear mixed models we will estimate separate specifications for:
  
  - the initially given data, 

- ranks of genes by expression, calculated for each participant at each time point as the quantile of the empirical distribution function for all genes' expression for this participant at this time point (for interpretability of results in logistic regressions we transformed ranks into the scale from 0 to 100), 

- $\gamma_0$ - regression intercept,

- $\gamma_1$ - regression coefficient.

Estimates of the regression coefficient for each gene and time point will be exponentiated, and after that we will be able to interpret them as the odds ratios of response in respect to the increase in the gene expression by 1 at the corresponding time point. 


```r
# Оценка логистических регрессий
glm_q <- purrr::quietly(.f = glm)
logreg_fits_SDY1260 <- df_expr_wide_fin %>%
  nest(-gene) %>%
  mutate(fit_init0 = map(data, ~ glm_q(response ~ expr_0, .x, family = "binomial")),
         fit_init3 = map(data, ~ glm_q(response ~ expr_3, .x, family = "binomial")),
         fit_init7 = map(data, ~ glm_q(response ~ expr_7, .x, family = "binomial")),
         fit_rank0 = map(data, ~ glm_q(response ~ expr_rank_0, .x, family = "binomial")),
         fit_rank3 = map(data, ~ glm_q(response ~ expr_rank_3, .x, family = "binomial")),
         fit_rank7 = map(data, ~ glm_q(response ~ expr_rank_7, .x, family = "binomial")))
logreg_res_SDY1260  <- logreg_fits_SDY1260  %>%
  mutate(
    gamma_p_init0 = map(fit_init0, ~ .x$result %>% broom.mixed::tidy(exp = TRUE)),
    gamma_p_init3 = map(fit_init3, ~ .x$result %>% broom.mixed::tidy(exp = TRUE)),
    gamma_p_init7 = map(fit_init7, ~ .x$result %>% broom.mixed::tidy(exp = TRUE)),
    gamma_p_rank0 = map(fit_rank0, ~ .x$result %>% broom.mixed::tidy(exp = TRUE)),
    gamma_p_rank3 = map(fit_rank3, ~ .x$result %>% broom.mixed::tidy(exp = TRUE)),
    gamma_p_rank7 = map(fit_rank7, ~ .x$result %>% broom.mixed::tidy(exp = TRUE))) %>%
  select(-data, -contains("fit"))
# # Сохраним результаты в rds
  saveRDS(logreg_res_SDY1260 , "logreg_res_SDY1260 .rds")
# # logreg_res <- readRDS("OlgaMironenko/res/logreg_res.rds")
# # logreg_res <- readRDS(file.path("..", "OlgaMironenko", "res", "logreg_res.rds"))
# Exp.gammas and p-values
logreg_gammas_SDY1260  <- logreg_res_SDY1260  %>%
  unnest(-gene, names_sep = "_") %>% 
  rename_with(function(x) gsub("gamma_p_", "", x), -gene)
  
# Сохраним результаты в rds, чтобы при формировании отчёта не ждать оценки регрессий
saveRDS(logreg_gammas_SDY1260 , "logreg_gammas_SDY1260 .rds")
```


```r
# logreg_gammas <- readRDS("OlgaMironenko/res/logreg_gammas.rds")
logreg_gammas_SDY1260  <- readRDS("logreg_gammas_SDY1260.rds")
critv <- -log10(0.05)
# Data frame for all exp.gammas and p-values for all LMMs
logreg_res_long <- logreg_gammas_SDY1260  %>%
  filter(grepl("expr", init0_term)) %>%
  select(gene, contains("estimate"), contains("p.value")) %>%
  pivot_longer(cols = -gene, names_pattern = "(init|rank)(\\d)_(estimate|p.value)",
               names_to = c("expr", "time", ".value")) %>%
  mutate(logp = -log10(p.value),
         sig = logp > critv,
         model_var = sprintf("log_%s_%s", expr, time),
         expr = factor(expr, c("init", "rank"), c("Expression as is", "Ranks by expression")),
         time = factor(as.numeric(time), timepoints, labels = ifelse(timepoints == 0, "Baseline", paste(timepoints, "d."))))
# Significant genes
genes_sig_logreg <- map(
  logreg_res_long %>%
    filter(sig) %>%
    split(f = as.factor(.$model_var)),
  ~ .x %>% pull(gene))
models_logreg <- names(genes_sig_logreg)
models_logreg_num <- as.numeric(str_extract(models_logreg, "\\d+"))
genes_sig_logreg_lbls <- sprintf("%s, sig.OR (%s)",
                                 ifelse(grepl("init", models_logreg), "Expr.", "Ranks"),
                                 ifelse(models_logreg_num == 0, "Baseline",
                                        sprintf("%d d.", models_logreg_num))) %>%
  as.list() %>% setNames(models_logreg)
```

<br>

### **Results for logistic regressions**

<br>

#### **Sets of genes with expression related to high response**

<br>

Figures 11(1) and 11(2) shows how much **the sets of genes with significant coefficients in logistic regressions are intersected**.


```r
# Whether each gene was significant in any model
genes_maxvar <- cbind(
  genes_maxvar,
  map_dfc(genes_sig_logreg, ~ genes_maxvar$gene %in% .x))
upset_plt_df_1 <- genes_maxvar %>%
  select(contains("log_init")) %>%
  filter(if_any(everything(), ~.)) %>%
  setNames(genes_sig_logreg_lbls[names(.)])
upset_plot(upset_plt_df_1, rev(colnames(upset_plt_df_1)), 
           rep("#F1CE63", 4),
           "Figure 11(1). Intersection of gene sets with significant ORs,\nExpression as is",
           "Logistic regressions for the (strong) response", TRUE, -0.8)
```

![](Project_ImmuneSpace_SDY_1260_files/figure-html/fig11_upset_logreg-1.png)<!-- -->

```r
upset_plt_df_2 <- genes_maxvar %>%
  select(contains("log_rank")) %>%
  filter(if_any(everything(), ~.)) %>%
  setNames(genes_sig_logreg_lbls[names(.)])
upset_plot(upset_plt_df_2, rev(colnames(upset_plt_df_2)), 
           rep("#A0CBE8", 4),
           "Figure 11(2). Intersection of gene sets with significant ORs,\nRanks by expression",
           "Logistic regressions for the (strong) response", TRUE, -0.8)
```

![](Project_ImmuneSpace_SDY_1260_files/figure-html/fig11_upset_logreg-2.png)<!-- -->

Also for every time point we can **compare gene sets with significant coefficient for expression in the logistic regression and gene sets with significant ATE at the corresponding time point as a result of the linear mixed model** - results are shown at Fugures 12.:


```r
upset_plt_df_1 <- genes_maxvar %>%
  select(init_int_ame0, log_init_0) %>%
  filter(if_any(everything(), ~.)) %>%
  setNames(c(genes_sig_lbls, genes_sig_logreg_lbls)[names(.)])
p1 <- upset_plot(upset_plt_df_1, rev(colnames(upset_plt_df_1)), 
                 rep("#F1CE63", 2),
                 element_blank(), "Expression as is", TRUE, 0.4, 1, FALSE)
upset_plt_df_2 <- genes_maxvar %>%
  select(rank_int_ame0, log_rank_0) %>%
  filter(if_any(everything(), ~.)) %>%
  setNames(c(genes_sig_lbls, genes_sig_logreg_lbls)[names(.)])
p2 <- upset_plot(upset_plt_df_2, rev(colnames(upset_plt_df_2)), 
                 rep("#A0CBE8", 2),
                 element_blank(), "Ranks by expression", TRUE, 0.4, 2, FALSE)
p12 <- ggarrange(p1, p2, nrow = 2)
annotate_figure(p12, top = text_grob("Figure 12(1). Intersection of gene sets with significant ATEs and ORs, Baseline",
                                     face = "bold", size = 12, hjust = 0.5))
```

![](Project_ImmuneSpace_SDY_1260_files/figure-html/fig12_compmodels-1.png)<!-- -->

```r
upset_plt_df_1 <- genes_maxvar %>%
  select(init_int_ame3, log_init_3) %>%
  filter(if_any(everything(), ~.)) %>%
  setNames(c(genes_sig_lbls, genes_sig_logreg_lbls)[names(.)])
p1 <- upset_plot(upset_plt_df_1, rev(colnames(upset_plt_df_1)), 
                 rep("#F1CE63", 2),
                 element_blank(), "Expression as is", TRUE, 0.2, 1, FALSE)
upset_plt_df_2 <- genes_maxvar %>%
  select(rank_int_ame3, log_rank_3) %>%
  filter(if_any(everything(), ~.)) %>%
  setNames(c(genes_sig_lbls, genes_sig_logreg_lbls)[names(.)])
p2 <- upset_plot(upset_plt_df_2, rev(colnames(upset_plt_df_2)), 
                 rep("#A0CBE8", 2),
                 element_blank(), "Ranks by expression", TRUE, 0.2, 2, FALSE)
p12 <- ggarrange(p1, p2, nrow = 2)
annotate_figure(p12, top = text_grob("Figure 12(3). Intersection of gene sets with significant ATEs and ORs, 3 days",
                                     face = "bold", size = 12, hjust = 0.5))
```

![](Project_ImmuneSpace_SDY_1260_files/figure-html/fig12_compmodels-2.png)<!-- -->

```r
upset_plt_df_1 <- genes_maxvar %>%
  select(init_int_ame7, log_init_7) %>%
  filter(if_any(everything(), ~.)) %>%
  setNames(c(genes_sig_lbls, genes_sig_logreg_lbls)[names(.)])
p1 <- upset_plot(upset_plt_df_1, rev(colnames(upset_plt_df_1)), 
                 rep("#F1CE63", 2),
                 element_blank(), "Expression as is", TRUE, 0.2, 1, FALSE)
upset_plt_df_2 <- genes_maxvar %>%
  select(rank_int_ame7, log_rank_7) %>%
  filter(if_any(everything(), ~.)) %>%
  setNames(c(genes_sig_lbls, genes_sig_logreg_lbls)[names(.)])
p2 <- upset_plot(upset_plt_df_2, rev(colnames(upset_plt_df_2)), 
                 rep("#A0CBE8", 2),
                 element_blank(), "Ranks by expression", TRUE, 0.2, 2, FALSE)
p12 <- ggarrange(p1, p2, nrow = 2)
annotate_figure(p12, top = text_grob("Figure 12(4). Intersection of gene sets with significant ATEs and ORs, 7 days",
                                     face = "bold", size = 12, hjust = 0.5))
```

![](Project_ImmuneSpace_SDY_1260_files/figure-html/fig12_compmodels-3.png)<!-- -->

We see that the gene sets are much more intersected at day 3 and number of intersected genes decreases to day 7.

<br>

#### **Significance vs. size of the response odds ratio**

<br>

We can estimate the size of the effect of the gene expression at the particular point of time on the probability of high response using **odds rations**, OR (exponentiated coefficients from the logistic regression). 

Figure 12 contains volcano plot matching the effect size ($log_{10}$ of the OR on the X asis) with its statistical significance ($-log_{10}$-transformed p-value on the Y axis). For the purpose of visualiazation in the results for regressions with the initial data on expressions we replaced $log_{10}$ values of ORs out of the interval [-10, 10] with the nearest border of this interval, we did the same with the $log_{10}$ values of ORs out of the interval [-4, 4] from the results in regressions with ranks of genes by expression. We colored points for 10 genes with the lowest and largest significant ORs in orange.


```r
logreg_v_plt <- logreg_res_long %>%
  mutate(estimate = log10(estimate),
         estimate = ifelse(expr == "Expression as is",
                           ifelse(estimate < -10, -10, 
                                  ifelse(estimate > 10, 10, estimate)),
                           ifelse(estimate < -4, -4, 
                                  ifelse(estimate > 4, 4, estimate)))) %>%
  group_by(expr, time, sig) %>%
  arrange(estimate) %>%
  mutate(sig_plot = row_number() %in% c(1:10, (n()-4):n())) %>%
  ungroup() %>%
  mutate(sig_plot = ifelse(!sig, FALSE, sig_plot))
ggplot(logreg_v_plt, aes(x = estimate, y = logp)) +
  geom_point(aes(color = sig_plot), alpha = 0.5, show.legend = FALSE) +
  geom_hline(aes(yintercept = 1), linewidth = 0.7, color = "#E15759", linetype = "dashed") +
  geom_hline(aes(yintercept = critv), linewidth = 0.7, color = "red") +
  geom_text_repel(aes(label = gene), logreg_v_plt %>% filter(sig_plot), size = 3) +
  scale_x_continuous(expand = c(0.02, 0)) +
  scale_y_continuous(expand = c(0.02, 0)) +
  scale_color_manual(values = c("#86BCB6", "#F28E2B")) +
  facet_rep_grid(time ~ expr, repeat.tick.labels = TRUE, scales = "free", switch = "y") +
  labs(x = bquote(log[10](OR)), y = bquote(-log[10](p[OR])), 
       title = "Figure 13. P-value vs. size of expression's OR of response",
       subtitle = bquote(log[10]~"odds ratios after logistic regressions for the (strong) response"),
       caption = "Solid red line is for p = 0.05, dashed red line is for p = 0.1.\nThe higher the Y axis value, the lower is the p-value.") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        panel.grid.major = element_line(linewidth = .2, color = '#ebebebFF'),
        panel.grid.minor = element_line(linewidth = .1, color = '#ebebebFF'),
        strip.background = element_blank(),
        strip.text = element_text(size = 12, color = "black", face = "bold"),
        strip.placement = "outside",
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10),
        plot.caption = element_text(size = 10, color = "black", face = "italic", hjust = 0))
```

![](Project_ImmuneSpace_SDY_1260_files/figure-html/fig13_volcano_or-1.png)<!-- -->

<br>

## **Functional analysis**

<br>

We have already used Gene Ontology (GO) to obtain annotations for biological processes (BPs) into which products of the differentially expressed genes are involved. However, there can be many findings (i.e.many DEGs), and it is both inconvenient and incorrect to look for annotations for each of them, because products of the same gene can participate in different paths, as well as products from different genes can be involved into the same process. Approaches from the **singular enrichment analysis (SEA)** [@tipney_introduction_2010] allows to determine which particular processes are overrepresented by the genes in some set in comparison with the full gene set (the latter comrises 5 thousand genes in our case). 

In our research we will try to reveal such overrepresented BPs using **p-values from the hypergeometric distribution**, i.e. probabilities of getting the same or even more representation of the given BP in the set of $n$ genes in comparison with the full set of $N$ genes [@boyle_gotermfinder_open_2004]. These p-values can be calculated with the `hyperGtest` function from the `Category` package in R [@gentleman_r_category_2022]. We will set the option `conditional` to TRUE in it. According to the package vignette, in this case "`hyperGTest` function uses the structure of the GO graph to estimate for each term whether or not there is evidence beyond that which is provided by the term's children to call the term in question statistically overrepresented. The algorithm conditions on all child terms that are themselves significant at the specified p-value, odds ratio, minimum or maximum gene set size cutoff. Given a subgraph of the GO ontology, the terms with no child categories are tested first. Next the nodes whose children have already been tested are tested. If any of a given node's children tested significant, the appropriate conditioning is performed)" [@falcon_using_2007]. 

We use 0.05 as a p-value cut-off and set minimum and maximum number of genes in BPs in the gene universe to 10 and 500, correspondingly. After estimating the hypergeometric test we will exclude BPs with less than 5 genes from the chosen gene set.

We will use two approaches to **control FDR** here: p-value adjustment under Benjamini & Hochberg method and q-values estimations [@storey_statistical_2003].

SEA will be held separately for each set of significant genes, revealed after linear mixed models and logistic regressions.


```r
library(org.Hs.eg.db)
library(GOstats)
# org.Hs.eg.db - Genome wide annotation for Human
# EntrezID for genes (exclude those with NA)
genes_universe <- genes_maxvar %>%
  left_join(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, genes_maxvar$gene, "ENTREZID", "SYMBOL") %>%
              distinct(SYMBOL, .keep_all = TRUE),
            by = c("gene" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID))
# Conditional hypergeometric test for each set of significant genes
# https://biocorecrg.github.io/PHINDaccess_RNAseq_2020/functional_analysis.html
hg_res_SDY1260 <- map(
  genes_universe %>% dplyr::select(dplyr::starts_with("init"), dplyr::starts_with("rank"), 
                                   dplyr::starts_with("log"), -contains("b2b3")),
  function(x) {
    hg_params <- new(getClassDef("GOHyperGParams", package = "Category"),
                     geneIds = genes_universe$ENTREZID[x],
                     universeGeneIds = genes_universe$ENTREZID,
                     annotation = "org.Hs.eg",
                     ontology = "BP",
                     pvalueCutoff = 0.05,
                     conditional = TRUE,
                     minSizeCutoff = 10, maxSizeCutoff = 500,
                     testDirection = "over")
    Category::hyperGTest(hg_params)
  })
# # Saving results into rds (in case of necessity, without sending to github)
saveRDS(hg_res_SDY1260, "hg_res_SDY1260.rds")
# 
# # hg_res <- readRDS("OlgaMironenko/res/hg_res.rds")
# hg_res <- readRDS(file.path("..", "OlgaMironenko", "res", "hg_res.rds"))
ann_select_q <- purrr::quietly(.f = AnnotationDbi::select)
# Genes in significant BPs
hg_res_sig_genes <- imap(
  hg_res_SDY1260, 
  function(x, y) {
    imap_dfr(Category::geneIdsByCategory(x, summary(x)$GOBPID),
             ~ ann_select_q(org.Hs.eg.db::org.Hs.eg.db, .x, "SYMBOL", "ENTREZID") %>%
               .$result %>% 
               left_join(lmer_res_long %>% filter(model_var == y),
                         by = c("SYMBOL" = "gene")) %>%
               transmute(GOID = .y, genes = paste(SYMBOL, collapse = "/"),
                         up = sum(estimate > 0), down = sum(estimate <= 0)) %>%
               unique())
  })
# Results for hypergeom.test for each model
hg_res_pq_SDY1260 <- imap(
  hg_res_SDY1260,
  ~ cbind(pvalue = Category::pvalues(.x),
          Count = Category::geneCounts(.x))%>%
    as_tibble(rownames = "GOID") %>%
    filter(Count > 4) %>%
    # adjusting p-values and calculate qvalues
    mutate(GeneRatio = Count/length(Category::geneIds(.x)),
           p.adjust = p.adjust(pvalue, method = "BH"),
           qvalue = qvalue(pvalue) %>% .$qvalues) %>%
    # joining BP terms and definitions
    left_join(
      AnnotationDbi::select(GO.db::GO.db, .$GOID, c("TERM", "DEFINITION"), keytype = "GOID"),
      by = "GOID") %>%
    # joining genes for significant BPs
    left_join(
      hg_res_sig_genes[[.y]], by = "GOID"))
# Saving the main results into rds for not waiting for results while knitting
saveRDS(hg_res_pq_SDY1260, "hg_res_pq_SDY1260.rds")
```

**The number of the overrepresented BPs** obtained in every model with different cut-offs for the adjusted p-values and q-values is listed in Table 4 (empty cells state for zero BPs).


```r
# hg_res_pq <- readRDS("OlgaMironenko/res/hg_res_pq.rds")
hg_res_pq_SDY1260 <- readRDS("hg_res_pq_SDY1260.rds")
imap_dfr(
  hg_res_pq_SDY1260[!grepl("noint", names(hg_res_pq_SDY1260))], 
  ~ tibble(model = c(genes_sig_lbls, genes_sig_logreg_lbls)[[.y]],
           p01 = sum(.x$p.adjust < 0.1),
           p005 = sum(.x$p.adjust < 0.05),
           p001 = sum(.x$p.adjust < 0.01),
           p0001 = sum(.x$p.adjust < 0.001),
           q01 = sum(.x$qvalue < 0.1),
           q005 = sum(.x$qvalue < 0.05),
           q001 = sum(.x$qvalue < 0.01),
           q0001 = sum(.x$qvalue < 0.001))) %>%
  mutate_if(is.numeric, ~ ifelse(. == 0, NA, .)) %>%
  kable(align = "lcccccccc",
        col.names = c("Gene set source", rep(paste("< ", c(0.1,0.05,0.01,0.001)), 2)),
        caption = "<b>Table 4. Number of significant biological processes (GO) found in hypergeometric tests</b>") %>%
  kableExtra::add_header_above(c(" " = 1, "Adjusted p-value" = 4, "q-value" = 4)) %>% 
  kableExtra::row_spec(0, bold = TRUE) %>%
  kableExtra::kable_classic(full_width = FALSE, position = "left", font_size = 14,
                            html_font = "\"Source Sans Pro\", helvetica, sans-serif")
```

<table class=" lightable-classic" style='font-size: 14px; font-family: "Source Sans Pro", helvetica, sans-serif; width: auto !important; '>
<caption style="font-size: initial !important;"><b>Table 4. Number of significant biological processes (GO) found in hypergeometric tests</b></caption>
 <thead>
<tr>
<th style="empty-cells: hide;border-bottom:hidden;" colspan="1"></th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="4"><div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">Adjusted p-value</div></th>
<th style="border-bottom:hidden;padding-bottom:0; padding-left:3px;padding-right:3px;text-align: center; " colspan="4"><div style="border-bottom: 1px solid #ddd; padding-bottom: 5px; ">q-value</div></th>
</tr>
  <tr>
   <th style="text-align:left;font-weight: bold;"> Gene set source </th>
   <th style="text-align:center;font-weight: bold;"> &lt;  0.1 </th>
   <th style="text-align:center;font-weight: bold;"> &lt;  0.05 </th>
   <th style="text-align:center;font-weight: bold;"> &lt;  0.01 </th>
   <th style="text-align:center;font-weight: bold;"> &lt;  0.001 </th>
   <th style="text-align:center;font-weight: bold;"> &lt;  0.1 </th>
   <th style="text-align:center;font-weight: bold;"> &lt;  0.05 </th>
   <th style="text-align:center;font-weight: bold;"> &lt;  0.01 </th>
   <th style="text-align:center;font-weight: bold;"> &lt;  0.001 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Expr., sig. ATE (Baseline) </td>
   <td style="text-align:center;"> 35 </td>
   <td style="text-align:center;"> 18 </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 374 </td>
   <td style="text-align:center;"> 95 </td>
   <td style="text-align:center;"> 14 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expr., sig. ATE (3 d.) </td>
   <td style="text-align:center;"> 24 </td>
   <td style="text-align:center;"> 13 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 256 </td>
   <td style="text-align:center;"> 58 </td>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expr., sig. ATE (7 d.) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 248 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expr., sig. b2 (Baseline) </td>
   <td style="text-align:center;"> 18 </td>
   <td style="text-align:center;"> 6 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 218 </td>
   <td style="text-align:center;"> 85 </td>
   <td style="text-align:center;"> 6 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expr., sig. b3 (Change) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ranks, sig. ATE (Baseline) </td>
   <td style="text-align:center;"> 13 </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 52 </td>
   <td style="text-align:center;"> 32 </td>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ranks, sig. ATE (3 d.) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 979 </td>
   <td style="text-align:center;"> 190 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ranks, sig. ATE (7 d.) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ranks, sig. b2 (Baseline) </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 50 </td>
   <td style="text-align:center;"> 21 </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ranks, sig. b3 (Change) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expr., sig.OR (Baseline) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 5 </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expr., sig.OR (3 d.) </td>
   <td style="text-align:center;"> 121 </td>
   <td style="text-align:center;"> 75 </td>
   <td style="text-align:center;"> 15 </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 645 </td>
   <td style="text-align:center;"> 358 </td>
   <td style="text-align:center;"> 78 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Expr., sig.OR (7 d.) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ranks, sig.OR (Baseline) </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ranks, sig.OR (3 d.) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;"> 774 </td>
   <td style="text-align:center;"> 34 </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Ranks, sig.OR (7 d.) </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
   <td style="text-align:center;">  </td>
  </tr>
</tbody>
</table>

Mostly all overrepresented biological processes are found for the Baseline and 1-day post-vaccination.


```r
gobp_overrep <- imap_dfr(
  hg_res_pq_SDY1260
  [!grepl("noint", names(hg_res_pq_SDY1260))],
  ~ .x %>% filter(qvalue < 0.05) %>% 
    transmute(model = c(genes_sig_lbls, genes_sig_logreg_lbls)[[.y]], 
              GOID, TERM, DEFINITION, Count, up, down, genes)) %>%
  mutate(time = str_remove_all(str_extract(model, "\\(.+\\)"), "[\\(\\)]"),
         imm = as.numeric(GOID %in% gobp_immresp$GOBPID)) %>%
  select(time, GOID, TERM, DEFINITION, imm) %>%
  unique() %>%
  group_by(GOID, TERM, DEFINITION, imm) %>%
  summarise(time = paste(time, collapse = ", ")) %>%
  dplyr::select(time, everything()) %>%
  arrange(-imm, time, GOID)
g1 <- which(gobp_overrep$imm == 1)
```

Table 5 contains all **BPs with q-value < 0.05**, grouped by time points when they were overrepresented at least in one model. There are 28 among them which are directly related to the immune system (particularly to the innate immunity), namely: cytokine production, positive regulation of cytokine production, activation of immune response, immune system development, regulation of immune system process, negative regulation of immune system process, regulation of T cell mediated immunity, inflammatory response, hemopoiesis, regulation of interleukin-10 production, positive regulation of tumor necrosis factor production, T cell proliferation, regulation of T cell proliferation, negative regulation of T cell proliferation, positive regulation of myeloid cell differentiation, regulation of alpha-beta T cell proliferation, lymphocyte proliferation, spleen development, negative regulation of T cell activation, regulation of lymphocyte activation, negative regulation of lymphocyte activation, myeloid cell development, antimicrobial humoral immune response mediated by antimicrobial peptide, granulocyte chemotaxis, myeloid leukocyte migration, positive regulation of leukocyte differentiation, positive regulation of tumor necrosis factor superfamily cytokine production, megakaryocyte differentiation - they are marked in bold in the table.


```r
gobp_overrep %>%
  select(-imm) %>%
  kable(align = "lcll",
        col.names = c("Time", "GOBPID", "BP (GO term)", "Definition (GO)"),
        caption = "<b>Table 5. Overrepresented BPs, by time</b>") %>%
  kableExtra::row_spec(c(0, g1), bold = TRUE) %>%
  kableExtra::column_spec(1:4, extra_css = "vertical-align:top;") %>%
  kableExtra::kable_classic(full_width = FALSE, position = "left", font_size = 14,
                            html_font = "\"Source Sans Pro\", helvetica, sans-serif")
```

<table class=" lightable-classic" style='font-size: 14px; font-family: "Source Sans Pro", helvetica, sans-serif; width: auto !important; '>
<caption style="font-size: initial !important;"><b>Table 5. Overrepresented BPs, by time</b></caption>
 <thead>
  <tr>
   <th style="text-align:left;font-weight: bold;"> Time </th>
   <th style="text-align:center;font-weight: bold;"> GOBPID </th>
   <th style="text-align:left;font-weight: bold;"> BP (GO term) </th>
   <th style="text-align:left;font-weight: bold;"> Definition (GO) </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0001816 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> cytokine production </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> The appearance of a cytokine due to biosynthesis or secretion following a cellular stimulus, resulting in an increase in its intracellular or extracellular levels. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0001819 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> positive regulation of cytokine production </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of production of a cytokine. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0002253 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> activation of immune response </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> Any process that initiates an immune response. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0002520 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> immune system development </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> The process whose specific outcome is the progression of an organismal system whose objective is to provide calibrated responses by an organism to a potential internal or invasive threat, over time, from its formation to the mature structure. A system is a regularly interacting or interdependent group of organs or tissues that work together to carry out a given biological process. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0002682 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> regulation of immune system process </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> Any process that modulates the frequency, rate, or extent of an immune system process. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0002683 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> negative regulation of immune system process </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate, or extent of an immune system process. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0002709 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> regulation of T cell mediated immunity </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> Any process that modulates the frequency, rate, or extent of T cell mediated immunity. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0006954 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> inflammatory response </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> The immediate defensive reaction (by vertebrate tissue) to infection or injury caused by chemical or physical agents. The process is characterized by local vasodilation, extravasation of plasma into intercellular spaces and accumulation of white blood cells and macrophages. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0030097 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> hemopoiesis </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> The process whose specific outcome is the progression of the myeloid and lymphoid derived organ/tissue systems of the blood and other parts of the body over time, from formation to the mature structure. The site of hemopoiesis is variable during development, but occurs primarily in bone marrow or kidney in many adult vertebrates. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0032653 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> regulation of interleukin-10 production </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> Any process that modulates the frequency, rate, or extent of interleukin-10 production. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0032760 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> positive regulation of tumor necrosis factor production </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of tumor necrosis factor production. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0042098 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> T cell proliferation </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> The expansion of a T cell population by cell division. Follows T cell activation. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0042129 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> regulation of T cell proliferation </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> Any process that modulates the frequency, rate or extent of T cell proliferation. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0042130 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> negative regulation of T cell proliferation </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> Any process that stops, prevents or reduces the rate or extent of T cell proliferation. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0045639 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> positive regulation of myeloid cell differentiation </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of myeloid cell differentiation. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0046640 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> regulation of alpha-beta T cell proliferation </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> Any process that modulates the frequency, rate or extent of alpha-beta T cell proliferation. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0046651 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> lymphocyte proliferation </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> The expansion of a lymphocyte population by cell division. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0048536 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> spleen development </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> The process whose specific outcome is the progression of the spleen over time, from its formation to the mature structure. The spleen is a large vascular lymphatic organ composed of white and red pulp, involved both in hemopoietic and immune system functions. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0050868 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> negative regulation of T cell activation </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of T cell activation. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0051249 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> regulation of lymphocyte activation </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> Any process that modulates the frequency, rate or extent of lymphocyte activation. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0051250 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> negative regulation of lymphocyte activation </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of lymphocyte activation. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0061515 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> myeloid cell development </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> The process whose specific outcome is the progression of a myeloid cell over time, from its formation to the mature structure. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0061844 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> antimicrobial humoral immune response mediated by antimicrobial peptide </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> An immune response against microbes mediated by anti-microbial peptides in body fluid. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0071621 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> granulocyte chemotaxis </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> The movement of a granulocyte in response to an external stimulus. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0097529 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> myeloid leukocyte migration </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> The movement of a myeloid leukocyte within or between different tissues and organs of the body. </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:1902107 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> positive regulation of leukocyte differentiation </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:1903557 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> positive regulation of tumor necrosis factor superfamily cytokine production </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> 3 d., Baseline </td>
   <td style="text-align:center;font-weight: bold;vertical-align:top;"> GO:0030219 </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> megakaryocyte differentiation </td>
   <td style="text-align:left;font-weight: bold;vertical-align:top;"> The process in which a myeloid precursor cell acquires specializes features of a megakaryocyte. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0000375 </td>
   <td style="text-align:left;vertical-align:top;"> RNA splicing, via transesterification reactions </td>
   <td style="text-align:left;vertical-align:top;"> Splicing of RNA via a series of two transesterification reactions. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0000398 </td>
   <td style="text-align:left;vertical-align:top;"> mRNA splicing, via spliceosome </td>
   <td style="text-align:left;vertical-align:top;"> The joining together of exons from one or more primary transcripts of messenger RNA (mRNA) and the excision of intron sequences, via a spliceosomal mechanism, so that mRNA consisting only of the joined exons is produced. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0000902 </td>
   <td style="text-align:left;vertical-align:top;"> cell morphogenesis </td>
   <td style="text-align:left;vertical-align:top;"> The developmental process in which the size or shape of a cell is generated and organized. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0000904 </td>
   <td style="text-align:left;vertical-align:top;"> cell morphogenesis involved in differentiation </td>
   <td style="text-align:left;vertical-align:top;"> The change in form (cell shape and size) that occurs when relatively unspecialized cells, e.g. embryonic or regenerative cells, acquire specialized structural and/or functional features that characterize the cells, tissues, or organs of the mature organism or some other relatively stable phase of the organism's life history. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0001501 </td>
   <td style="text-align:left;vertical-align:top;"> skeletal system development </td>
   <td style="text-align:left;vertical-align:top;"> The process whose specific outcome is the progression of the skeleton over time, from its formation to the mature structure. The skeleton is the bony framework of the body in vertebrates (endoskeleton) or the hard outer envelope of insects (exoskeleton or dermoskeleton). </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0001666 </td>
   <td style="text-align:left;vertical-align:top;"> response to hypoxia </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a stimulus indicating lowered oxygen tension. Hypoxia, defined as a decline in O2 levels below normoxic levels of 20.8 - 20.95%, results in metabolic adaptation at both the cellular and organismal level. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0001676 </td>
   <td style="text-align:left;vertical-align:top;"> long-chain fatty acid metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving long-chain fatty acids, A long-chain fatty acid is a fatty acid with a chain length between C13 and C22. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0001932 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of protein phosphorylation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of addition of phosphate groups into an amino acid in a protein. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0001933 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of protein phosphorylation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents or reduces the rate of addition of phosphate groups to amino acids within a protein. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0001936 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of endothelial cell proliferation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate, or extent of endothelial cell proliferation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0001953 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of cell-matrix adhesion </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the rate or extent of cell adhesion to the extracellular matrix. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0002027 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of heart rate </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency or rate of heart contraction. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0002040 </td>
   <td style="text-align:left;vertical-align:top;"> sprouting angiogenesis </td>
   <td style="text-align:left;vertical-align:top;"> The extension of new blood vessels from existing vessels into avascular tissues, this process includes the specialization of endothelial cells into leading tip and stalk cells, proliferation and migration of the endothelial cells and cell adhesion resulting in angiogenic sprout fusion or lumen formation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0002064 </td>
   <td style="text-align:left;vertical-align:top;"> epithelial cell development </td>
   <td style="text-align:left;vertical-align:top;"> The process whose specific outcome is the progression of an epithelial cell over time, from its formation to the mature structure. An epithelial cell is a cell usually found in a two-dimensional sheet with a free surface. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0003008 </td>
   <td style="text-align:left;vertical-align:top;"> system process </td>
   <td style="text-align:left;vertical-align:top;"> A multicellular organismal process carried out by any of the organs or tissues in an organ system. An organ system is a regularly interacting or interdependent group of organs or tissues that work together to carry out a biological objective. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0003012 </td>
   <td style="text-align:left;vertical-align:top;"> muscle system process </td>
   <td style="text-align:left;vertical-align:top;"> A organ system process carried out at the level of a muscle. Muscle tissue is composed of contractile cells or fibers. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0003015 </td>
   <td style="text-align:left;vertical-align:top;"> heart process </td>
   <td style="text-align:left;vertical-align:top;"> A circulatory system process carried out by the heart. The heart is a hollow, muscular organ, which, by contracting rhythmically, keeps up the circulation of the blood. The heart is a hollow, muscular organ, which, by contracting rhythmically, keeps up the circulation of the blood. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0003018 </td>
   <td style="text-align:left;vertical-align:top;"> vascular process in circulatory system </td>
   <td style="text-align:left;vertical-align:top;"> A circulatory process that occurs at the level of the vasculature. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0003158 </td>
   <td style="text-align:left;vertical-align:top;"> endothelium development </td>
   <td style="text-align:left;vertical-align:top;"> The process whose specific outcome is the progression of an endothelium over time, from its formation to the mature structure. Endothelium refers to the layer of cells lining blood vessels, lymphatics, the heart, and serous cavities, and is derived from bone marrow or mesoderm. Corneal endothelium is a special case, derived from neural crest cells. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006139 </td>
   <td style="text-align:left;vertical-align:top;"> nucleobase-containing compound metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> Any cellular metabolic process involving nucleobases, nucleosides, nucleotides and nucleic acids. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006631 </td>
   <td style="text-align:left;vertical-align:top;"> fatty acid metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving fatty acids, aliphatic monocarboxylic acids liberated from naturally occurring fats and oils by hydrolysis. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006636 </td>
   <td style="text-align:left;vertical-align:top;"> unsaturated fatty acid biosynthetic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways resulting in the formation of an unsaturated fatty acid, any fatty acid containing one or more double bonds between carbon atoms. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006725 </td>
   <td style="text-align:left;vertical-align:top;"> cellular aromatic compound metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving aromatic compounds, any organic compound characterized by one or more planar rings, each of which contains conjugated double bonds and delocalized pi electrons, as carried out by individual cells. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006749 </td>
   <td style="text-align:left;vertical-align:top;"> glutathione metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving glutathione, the tripeptide glutamylcysteinylglycine, which acts as a coenzyme for some enzymes and as an antioxidant in the protection of sulfhydryl groups in enzymes and other proteins; it has a specific role in the reduction of hydrogen peroxide (H2O2) and oxidized ascorbate, and it participates in the gamma-glutamyl cycle. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006753 </td>
   <td style="text-align:left;vertical-align:top;"> nucleoside phosphate metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving any phosphorylated nucleoside. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006793 </td>
   <td style="text-align:left;vertical-align:top;"> phosphorus metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving the nonmetallic element phosphorus or compounds that contain phosphorus, usually in the form of a phosphate group (PO4). </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006796 </td>
   <td style="text-align:left;vertical-align:top;"> phosphate-containing compound metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving the phosphate group, the anion or salt of any phosphoric acid. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006807 </td>
   <td style="text-align:left;vertical-align:top;"> nitrogen compound metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving organic or inorganic compounds that contain nitrogen. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006874 </td>
   <td style="text-align:left;vertical-align:top;"> cellular calcium ion homeostasis </td>
   <td style="text-align:left;vertical-align:top;"> Any process involved in the maintenance of an internal steady state of calcium ions at the level of a cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006886 </td>
   <td style="text-align:left;vertical-align:top;"> intracellular protein transport </td>
   <td style="text-align:left;vertical-align:top;"> The directed movement of proteins in a cell, including the movement of proteins between specific compartments or structures within a cell, such as organelles of a eukaryotic cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006892 </td>
   <td style="text-align:left;vertical-align:top;"> post-Golgi vesicle-mediated transport </td>
   <td style="text-align:left;vertical-align:top;"> The directed movement of substances from the Golgi to other parts of the cell, including organelles and the plasma membrane, mediated by small transport vesicles. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006897 </td>
   <td style="text-align:left;vertical-align:top;"> endocytosis </td>
   <td style="text-align:left;vertical-align:top;"> A vesicle-mediated transport process in which cells take up external materials or membrane constituents by the invagination of a small region of the plasma membrane to form a new membrane-bounded vesicle. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006909 </td>
   <td style="text-align:left;vertical-align:top;"> phagocytosis </td>
   <td style="text-align:left;vertical-align:top;"> A vesicle-mediated transport process that results in the engulfment of external particulate material by phagocytes and their delivery to the lysosome. The particles are initially contained within phagocytic vacuoles (phagosomes), which then fuse with primary lysosomes to effect digestion of the particles. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006935 </td>
   <td style="text-align:left;vertical-align:top;"> chemotaxis </td>
   <td style="text-align:left;vertical-align:top;"> The directed movement of a motile cell or organism, or the directed growth of a cell guided by a specific chemical concentration gradient. Movement may be towards a higher concentration (positive chemotaxis) or towards a lower concentration (negative chemotaxis). </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006936 </td>
   <td style="text-align:left;vertical-align:top;"> muscle contraction </td>
   <td style="text-align:left;vertical-align:top;"> A process in which force is generated within muscle tissue, resulting in a change in muscle geometry. Force generation involves a chemo-mechanical energy conversion step that is carried out by the actin/myosin complex activity, which generates force through ATP hydrolysis. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006937 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of muscle contraction </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of muscle contraction. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006939 </td>
   <td style="text-align:left;vertical-align:top;"> smooth muscle contraction </td>
   <td style="text-align:left;vertical-align:top;"> A process in which force is generated within smooth muscle tissue, resulting in a change in muscle geometry. Force generation involves a chemo-mechanical energy conversion step that is carried out by the actin/myosin complex activity, which generates force through ATP hydrolysis. Smooth muscle differs from striated muscle in the much higher actin/myosin ratio, the absence of conspicuous sarcomeres and the ability to contract to a much smaller fraction of its resting length. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006950 </td>
   <td style="text-align:left;vertical-align:top;"> response to stress </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a disturbance in organismal or cellular homeostasis, usually, but not necessarily, exogenous (e.g. temperature, humidity, ionizing radiation). </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007043 </td>
   <td style="text-align:left;vertical-align:top;"> cell-cell junction assembly </td>
   <td style="text-align:left;vertical-align:top;"> The aggregation, arrangement and bonding together of a set of components to form a junction between cells. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007093 </td>
   <td style="text-align:left;vertical-align:top;"> mitotic cell cycle checkpoint signaling </td>
   <td style="text-align:left;vertical-align:top;"> A signaling process that ensures accurate chromosome replication and segregation by preventing progression through a mitotic cell cycle until conditions are suitable for the cell to proceed to the next stage. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007154 </td>
   <td style="text-align:left;vertical-align:top;"> cell communication </td>
   <td style="text-align:left;vertical-align:top;"> Any process that mediates interactions between a cell and its surroundings. Encompasses interactions such as signaling or attachment between one cell and another cell, between a cell and an extracellular matrix, or between a cell and any other aspect of its environment. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007155 </td>
   <td style="text-align:left;vertical-align:top;"> cell adhesion </td>
   <td style="text-align:left;vertical-align:top;"> The attachment of a cell, either to another cell or to an underlying substrate such as the extracellular matrix, via cell adhesion molecules. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007160 </td>
   <td style="text-align:left;vertical-align:top;"> cell-matrix adhesion </td>
   <td style="text-align:left;vertical-align:top;"> The binding of a cell to the extracellular matrix via adhesion molecules. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007162 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of cell adhesion </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of cell adhesion. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007165 </td>
   <td style="text-align:left;vertical-align:top;"> signal transduction </td>
   <td style="text-align:left;vertical-align:top;"> The cellular process in which a signal is conveyed to trigger a change in the activity or state of a cell. Signal transduction begins with reception of a signal (e.g. a ligand binding to a receptor or receptor activation by a stimulus such as light), or for signal transduction in the absence of ligand, signal-withdrawal or the activity of a constitutively active receptor. Signal transduction ends with regulation of a downstream cellular process, e.g. regulation of transcription or regulation of a metabolic process. Signal transduction covers signaling from receptors located on the surface of the cell and signaling via molecules located within the cell. For signaling between cells, signal transduction is restricted to events at and within the receiving cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007179 </td>
   <td style="text-align:left;vertical-align:top;"> transforming growth factor beta receptor signaling pathway </td>
   <td style="text-align:left;vertical-align:top;"> The series of molecular signals initiated by an extracellular ligand binding to a transforming growth factor beta receptor on the surface of a target cell, and ending with the regulation of a downstream cellular process, e.g. transcription. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007267 </td>
   <td style="text-align:left;vertical-align:top;"> cell-cell signaling </td>
   <td style="text-align:left;vertical-align:top;"> Any process that mediates the transfer of information from one cell to another. This process includes signal transduction in the receiving cell and, where applicable, release of a ligand and any processes that actively facilitate its transport and presentation to the receiving cell.  Examples include signaling via soluble ligands, via cell adhesion molecules and via gap junctions. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007272 </td>
   <td style="text-align:left;vertical-align:top;"> ensheathment of neurons </td>
   <td style="text-align:left;vertical-align:top;"> The process in which glial cells envelop neuronal cell bodies and/or axons to form an insulating layer. This can take the form of myelinating or non-myelinating ensheathment. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007411 </td>
   <td style="text-align:left;vertical-align:top;"> axon guidance </td>
   <td style="text-align:left;vertical-align:top;"> The chemotaxis process that directs the migration of an axon growth cone to a specific target site in response to a combination of attractive and repulsive cues. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007517 </td>
   <td style="text-align:left;vertical-align:top;"> muscle organ development </td>
   <td style="text-align:left;vertical-align:top;"> The process whose specific outcome is the progression of the muscle over time, from its formation to the mature structure. The muscle is an organ consisting of a tissue made up of various elongated cells that are specialized to contract and thus to produce movement and mechanical work. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007586 </td>
   <td style="text-align:left;vertical-align:top;"> digestion </td>
   <td style="text-align:left;vertical-align:top;"> The whole of the physical, chemical, and biochemical processes carried out by multicellular organisms to break down ingested nutrients into components that may be easily absorbed and directed into metabolism. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007599 </td>
   <td style="text-align:left;vertical-align:top;"> hemostasis </td>
   <td style="text-align:left;vertical-align:top;"> The stopping of bleeding (loss of body fluid) or the arrest of the circulation to an organ or part. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007610 </td>
   <td style="text-align:left;vertical-align:top;"> behavior </td>
   <td style="text-align:left;vertical-align:top;"> The internally coordinated responses (actions or inactions) of animals (individuals or groups) to internal or external stimuli, via a mechanism that involves nervous system activity. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007612 </td>
   <td style="text-align:left;vertical-align:top;"> learning </td>
   <td style="text-align:left;vertical-align:top;"> Any process in an organism in which a relatively long-lasting adaptive behavioral change occurs as the result of experience. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007613 </td>
   <td style="text-align:left;vertical-align:top;"> memory </td>
   <td style="text-align:left;vertical-align:top;"> The activities involved in the mental information processing system that receives (registers), modifies, stores, and retrieves informational stimuli. The main stages involved in the formation and retrieval of memory are encoding (processing of received information by acquisition), storage (building a permanent record of received information as a result of consolidation) and retrieval (calling back the stored information and use it in a suitable way to execute a given task). </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0008152 </td>
   <td style="text-align:left;vertical-align:top;"> metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways, including anabolism and catabolism, by which living organisms transform chemical substances. Metabolic processes typically transform small molecules, but also include macromolecular processes such as DNA repair and replication, and protein synthesis and degradation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0008154 </td>
   <td style="text-align:left;vertical-align:top;"> actin polymerization or depolymerization </td>
   <td style="text-align:left;vertical-align:top;"> Assembly or disassembly of actin filaments by the addition or removal of actin monomers from a filament. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0008217 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of blood pressure </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the force with which blood travels through the circulatory system. The process is controlled by a balance of processes that increase pressure and decrease pressure. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0008284 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of cell population proliferation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the rate or extent of cell proliferation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0008610 </td>
   <td style="text-align:left;vertical-align:top;"> lipid biosynthetic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways resulting in the formation of lipids, compounds soluble in an organic solvent but not, or sparingly, in an aqueous solvent. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0008625 </td>
   <td style="text-align:left;vertical-align:top;"> extrinsic apoptotic signaling pathway via death domain receptors </td>
   <td style="text-align:left;vertical-align:top;"> The series of molecular signals in which a signal is conveyed from the cell surface to trigger the apoptotic death of a cell. The pathway starts with a ligand binding to a death domain receptor on the cell surface, and ends when the execution phase of apoptosis is triggered. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0008637 </td>
   <td style="text-align:left;vertical-align:top;"> apoptotic mitochondrial changes </td>
   <td style="text-align:left;vertical-align:top;"> The morphological and physiological alterations undergone by mitochondria during apoptosis. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0008643 </td>
   <td style="text-align:left;vertical-align:top;"> carbohydrate transport </td>
   <td style="text-align:left;vertical-align:top;"> The directed movement of carbohydrate into, out of or within a cell, or between cells, by means of some agent such as a transporter or pore. Carbohydrates are a group of organic compounds based of the general formula Cx(H2O)y. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0009064 </td>
   <td style="text-align:left;vertical-align:top;"> glutamine family amino acid metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving amino acids of the glutamine family, comprising arginine, glutamate, glutamine and proline. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0009165 </td>
   <td style="text-align:left;vertical-align:top;"> nucleotide biosynthetic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways resulting in the formation of nucleotides, any nucleoside that is esterified with (ortho)phosphate or an oligophosphate at any hydroxyl group on the glycose moiety; may be mono-, di- or triphosphate; this definition includes cyclic-nucleotides (nucleoside cyclic phosphates). </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0009605 </td>
   <td style="text-align:left;vertical-align:top;"> response to external stimulus </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of an external stimulus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0009612 </td>
   <td style="text-align:left;vertical-align:top;"> response to mechanical stimulus </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a mechanical stimulus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0009636 </td>
   <td style="text-align:left;vertical-align:top;"> response to toxic substance </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a toxic stimulus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0009880 </td>
   <td style="text-align:left;vertical-align:top;"> embryonic pattern specification </td>
   <td style="text-align:left;vertical-align:top;"> The process that results in the patterns of cell differentiation that will arise in an embryo. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0009887 </td>
   <td style="text-align:left;vertical-align:top;"> animal organ morphogenesis </td>
   <td style="text-align:left;vertical-align:top;"> Morphogenesis of an animal organ. An organ is defined as a tissue or set of tissues that work together to perform a specific function or functions. Morphogenesis is the process in which anatomical structures are generated and organized. Organs are commonly observed as visibly distinct structures, but may also exist as loosely associated clusters of cells that work together to perform a specific function or functions. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0009888 </td>
   <td style="text-align:left;vertical-align:top;"> tissue development </td>
   <td style="text-align:left;vertical-align:top;"> The process whose specific outcome is the progression of a tissue over time, from its formation to the mature structure. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0009892 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of the chemical reactions and pathways within a cell or an organism. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0009893 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of the chemical reactions and pathways within a cell or an organism. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0009967 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of signal transduction </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of signal transduction. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010033 </td>
   <td style="text-align:left;vertical-align:top;"> response to organic substance </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of an organic substance stimulus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010543 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of platelet activation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the rate or frequency of platelet activation. Platelet activation is a series of progressive, overlapping events triggered by exposure of the platelets to subendothelial tissue. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010562 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of phosphorus metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that increases the frequency, rate or extent of the chemical reactions and pathways involving phosphorus or compounds containing phosphorus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010594 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of endothelial cell migration </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the rate, frequency, or extent of the orderly movement of an endothelial cell into the extracellular matrix to form an endothelium. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010628 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of gene expression </td>
   <td style="text-align:left;vertical-align:top;"> Any process that increases the frequency, rate or extent of gene expression. Gene expression is the process in which a gene's coding sequence is converted into a mature gene product (protein or RNA). </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010634 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of epithelial cell migration </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of epithelial cell migration. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010638 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of organelle organization </td>
   <td style="text-align:left;vertical-align:top;"> Any process that increases the frequency, rate or extent of a process involved in the formation, arrangement of constituent parts, or disassembly of an organelle. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010639 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of organelle organization </td>
   <td style="text-align:left;vertical-align:top;"> Any process that decreases the frequency, rate or extent of a process involved in the formation, arrangement of constituent parts, or disassembly of an organelle. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010644 </td>
   <td style="text-align:left;vertical-align:top;"> cell communication by electrical coupling </td>
   <td style="text-align:left;vertical-align:top;"> The process that mediates signaling interactions between one cell and another cell by transfer of current between their adjacent cytoplasms via intercellular protein channels. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010646 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of cell communication </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of cell communication. Cell communication is the process that mediates interactions between a cell and its surroundings. Encompasses interactions such as signaling or attachment between one cell and another cell, between a cell and an extracellular matrix, or between a cell and any other aspect of its environment. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010647 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of cell communication </td>
   <td style="text-align:left;vertical-align:top;"> Any process that increases the frequency, rate or extent of cell communication. Cell communication is the process that mediates interactions between a cell and its surroundings. Encompasses interactions such as signaling or attachment between one cell and another cell, between a cell and an extracellular matrix, or between a cell and any other aspect of its environment. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010648 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of cell communication </td>
   <td style="text-align:left;vertical-align:top;"> Any process that decreases the frequency, rate or extent of cell communication. Cell communication is the process that mediates interactions between a cell and its surroundings. Encompasses interactions such as signaling or attachment between one cell and another cell, between a cell and an extracellular matrix, or between a cell and any other aspect of its environment. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010720 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of cell development </td>
   <td style="text-align:left;vertical-align:top;"> Any process that increases the rate, frequency or extent of the progression of the cell over time, from its formation to the mature structure. Cell development does not include the steps involved in committing a cell to a specific fate. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010721 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of cell development </td>
   <td style="text-align:left;vertical-align:top;"> Any process that decreases the rate, frequency or extent of the progression of the cell over time, from its formation to the mature structure. Cell development does not include the steps involved in committing a cell to a specific fate. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010799 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of peptidyl-threonine phosphorylation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of peptidyl-threonine phosphorylation. Peptidyl-threonine phosphorylation is the phosphorylation of peptidyl-threonine to form peptidyl-O-phospho-L-threonine. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010828 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of glucose transmembrane transport </td>
   <td style="text-align:left;vertical-align:top;"> Any process that increases the frequency, rate or extent of glucose transport across a membrane. Glucose transport is the directed movement of the hexose monosaccharide glucose into, out of or within a cell, or between cells, by means of some agent such as a transporter or pore. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010876 </td>
   <td style="text-align:left;vertical-align:top;"> lipid localization </td>
   <td style="text-align:left;vertical-align:top;"> Any process in which a lipid is transported to, or maintained in, a specific location. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010951 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of endopeptidase activity </td>
   <td style="text-align:left;vertical-align:top;"> Any process that decreases the frequency, rate or extent of endopeptidase activity, the endohydrolysis of peptide bonds within proteins. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010952 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of peptidase activity </td>
   <td style="text-align:left;vertical-align:top;"> Any process that increases the frequency, rate or extent of peptidase activity, the hydrolysis of peptide bonds within proteins. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010971 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of G2/M transition of mitotic cell cycle </td>
   <td style="text-align:left;vertical-align:top;"> Any signalling pathway that activates or increases the activity of a cell cycle cyclin-dependent protein kinase to modulate the switch from G2 phase to M phase of the mitotic cell cycle. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010977 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of neuron projection development </td>
   <td style="text-align:left;vertical-align:top;"> Any process that decreases the rate, frequency or extent of neuron projection development. Neuron projection development is the process whose specific outcome is the progression of a neuron projection over time, from its formation to the mature structure. A neuron projection is any process extending from a neural cell, such as axons or dendrites (collectively called neurites). </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0014068 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of phosphatidylinositol 3-kinase signaling </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of signal transduction mediated by the phosphatidylinositol 3-kinase cascade. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0014070 </td>
   <td style="text-align:left;vertical-align:top;"> response to organic cyclic compound </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of an organic cyclic compound stimulus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0014855 </td>
   <td style="text-align:left;vertical-align:top;"> striated muscle cell proliferation </td>
   <td style="text-align:left;vertical-align:top;"> The multiplication or reproduction of striated muscle cells, resulting in the expansion of a cell population. Striated muscles contain fibers that are divided by transverse bands into striations, and cardiac and skeletal muscle are types of striated muscle. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0014902 </td>
   <td style="text-align:left;vertical-align:top;"> myotube differentiation </td>
   <td style="text-align:left;vertical-align:top;"> The process in which a relatively unspecialized cell acquires specialized features of a myotube cell. Myotube differentiation starts with myoblast fusion and the appearance of specific cell markers (this is the cell development step). Then individual myotubes can fuse to form bigger myotubes and start to contract. Myotubes are multinucleated cells that are formed when proliferating myoblasts exit the cell cycle, differentiate and fuse. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0015908 </td>
   <td style="text-align:left;vertical-align:top;"> fatty acid transport </td>
   <td style="text-align:left;vertical-align:top;"> The directed movement of fatty acids into, out of or within a cell, or between cells, by means of some agent such as a transporter or pore. Fatty acids are aliphatic monocarboxylic acids liberated from naturally occurring fats and oils by hydrolysis. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0016310 </td>
   <td style="text-align:left;vertical-align:top;"> phosphorylation </td>
   <td style="text-align:left;vertical-align:top;"> The process of introducing a phosphate group into a molecule, usually with the formation of a phosphoric ester, a phosphoric anhydride or a phosphoric amide. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0016482 </td>
   <td style="text-align:left;vertical-align:top;"> cytosolic transport </td>
   <td style="text-align:left;vertical-align:top;"> The directed movement of substances or organelles within the cytosol. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0016525 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of angiogenesis </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of angiogenesis. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0018108 </td>
   <td style="text-align:left;vertical-align:top;"> peptidyl-tyrosine phosphorylation </td>
   <td style="text-align:left;vertical-align:top;"> The phosphorylation of peptidyl-tyrosine to form peptidyl-O4'-phospho-L-tyrosine. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0019220 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of phosphate metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of the chemical reactions and pathways involving phosphates. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0019222 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of the chemical reactions and pathways within a cell or an organism. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0019369 </td>
   <td style="text-align:left;vertical-align:top;"> arachidonic acid metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving arachidonic acid, a straight chain fatty acid with 20 carbon atoms and four double bonds per molecule. Arachidonic acid is the all-Z-(5,8,11,14)-isomer. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0019932 </td>
   <td style="text-align:left;vertical-align:top;"> second-messenger-mediated signaling </td>
   <td style="text-align:left;vertical-align:top;"> Any intracellular signal transduction in which the signal is passed on within the cell via a second messenger; a small molecule or ion that can be quickly generated or released from intracellular stores, and can diffuse within the cell. Second-messenger signaling includes production or release of the second messenger, and effectors downstream of the second messenger that further transmit the signal within the cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0021782 </td>
   <td style="text-align:left;vertical-align:top;"> glial cell development </td>
   <td style="text-align:left;vertical-align:top;"> The process aimed at the progression of a glial cell over time, from initial commitment of the cell to a specific fate, to the fully functional differentiated cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0022008 </td>
   <td style="text-align:left;vertical-align:top;"> neurogenesis </td>
   <td style="text-align:left;vertical-align:top;"> Generation of cells within the nervous system. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0022409 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of cell-cell adhesion </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the rate or extent of cell adhesion to another cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0022411 </td>
   <td style="text-align:left;vertical-align:top;"> cellular component disassembly </td>
   <td style="text-align:left;vertical-align:top;"> A cellular process that results in the breakdown of a cellular component. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0022600 </td>
   <td style="text-align:left;vertical-align:top;"> digestive system process </td>
   <td style="text-align:left;vertical-align:top;"> A physical, chemical, or biochemical process carried out by living organisms to break down ingested nutrients into components that may be easily absorbed and directed into metabolism. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0022603 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of anatomical structure morphogenesis </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of anatomical structure morphogenesis. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0022604 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of cell morphogenesis </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of cell morphogenesis. Cell morphogenesis is the developmental process in which the shape of a cell is generated and organized. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0022607 </td>
   <td style="text-align:left;vertical-align:top;"> cellular component assembly </td>
   <td style="text-align:left;vertical-align:top;"> The aggregation, arrangement and bonding together of a cellular component. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0023051 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of signaling </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of a signaling process. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0023052 </td>
   <td style="text-align:left;vertical-align:top;"> signaling </td>
   <td style="text-align:left;vertical-align:top;"> The entirety of a process in which information is transmitted within a biological system. This process begins with an active signal and ends when a cellular response has been triggered. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0023056 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of signaling </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates, maintains or increases the frequency, rate or extent of a signaling process. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0023057 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of signaling </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of a signaling process. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0030036 </td>
   <td style="text-align:left;vertical-align:top;"> actin cytoskeleton organization </td>
   <td style="text-align:left;vertical-align:top;"> A process that is carried out at the cellular level which results in the assembly, arrangement of constituent parts, or disassembly of cytoskeletal structures comprising actin filaments and their associated proteins. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0030038 </td>
   <td style="text-align:left;vertical-align:top;"> contractile actin filament bundle assembly </td>
   <td style="text-align:left;vertical-align:top;"> Assembly of actin filament bundles in which the filaments are loosely packed (approximately 30-60 nm apart) and arranged with opposing polarities; the loose packing allows myosin (usually myosin-II) to enter the bundle. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0030100 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of endocytosis </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of endocytosis. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0030155 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of cell adhesion </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of attachment of a cell to another cell or to the extracellular matrix. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0030168 </td>
   <td style="text-align:left;vertical-align:top;"> platelet activation </td>
   <td style="text-align:left;vertical-align:top;"> A series of progressive, overlapping events triggered by exposure of the platelets to subendothelial tissue. These events include shape change, adhesiveness, aggregation, and release reactions. When carried through to completion, these events lead to the formation of a stable hemostatic plug. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0030182 </td>
   <td style="text-align:left;vertical-align:top;"> neuron differentiation </td>
   <td style="text-align:left;vertical-align:top;"> The process in which a relatively unspecialized cell acquires specialized features of a neuron. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0030195 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of blood coagulation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of blood coagulation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0030323 </td>
   <td style="text-align:left;vertical-align:top;"> respiratory tube development </td>
   <td style="text-align:left;vertical-align:top;"> The process whose specific outcome is the progression of the respiratory tube over time, from its formation to the mature structure. The respiratory tube is assumed to mean any tube in the respiratory tract. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0030324 </td>
   <td style="text-align:left;vertical-align:top;"> lung development </td>
   <td style="text-align:left;vertical-align:top;"> The process whose specific outcome is the progression of the lung over time, from its formation to the mature structure. In all air-breathing vertebrates the lungs are developed from the ventral wall of the oesophagus as a pouch which divides into two sacs. In amphibians and many reptiles the lungs retain very nearly this primitive sac-like character, but in the higher forms the connection with the esophagus becomes elongated into the windpipe and the inner walls of the sacs become more and more divided, until, in the mammals, the air spaces become minutely divided into tubes ending in small air cells, in the walls of which the blood circulates in a fine network of capillaries. In mammals the lungs are more or less divided into lobes, and each lung occupies a separate cavity in the thorax. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0030335 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of cell migration </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of cell migration. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0030832 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of actin filament length </td>
   <td style="text-align:left;vertical-align:top;"> Any process that controls the length of actin filaments in a cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0030833 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of actin filament polymerization </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of the assembly of actin filaments by the addition of actin monomers to a filament. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0031324 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of cellular metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of the chemical reactions and pathways by which individual cells transform chemical substances. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0031325 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of cellular metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of the chemical reactions and pathways by which individual cells transform chemical substances. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0031334 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of protein-containing complex assembly </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of protein complex assembly. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0031589 </td>
   <td style="text-align:left;vertical-align:top;"> cell-substrate adhesion </td>
   <td style="text-align:left;vertical-align:top;"> The attachment of a cell to the underlying substrate via adhesion molecules. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0031667 </td>
   <td style="text-align:left;vertical-align:top;"> response to nutrient levels </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a stimulus reflecting the presence, absence, or concentration of nutrients. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0032355 </td>
   <td style="text-align:left;vertical-align:top;"> response to estradiol </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of stimulus by estradiol, a C18 steroid hormone hydroxylated at C3 and C17 that acts as a potent estrogen. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0032368 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of lipid transport </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of the directed movement of lipids into, out of or within a cell, or between cells, by means of some agent such as a transporter or pore. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0032409 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of transporter activity </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the activity of a transporter. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0032412 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of ion transmembrane transporter activity </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the activity of an ion transporter. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0032526 </td>
   <td style="text-align:left;vertical-align:top;"> response to retinoic acid </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a retinoic acid stimulus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0032879 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of localization </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of any process in which a cell, a substance, or a cellular entity is transported to, or maintained in, a specific location. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0032944 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of mononuclear cell proliferation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of mononuclear cell proliferation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0032970 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of actin filament-based process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of any cellular process that depends upon or alters the actin cytoskeleton. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0032989 </td>
   <td style="text-align:left;vertical-align:top;"> cellular component morphogenesis </td>
   <td style="text-align:left;vertical-align:top;"> The process in which cellular structures, including whole cells or cell parts, are generated and organized. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0032990 </td>
   <td style="text-align:left;vertical-align:top;"> cell part morphogenesis </td>
   <td style="text-align:left;vertical-align:top;"> The process in which the anatomical structures of a cell part are generated and organized. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0033043 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of organelle organization </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of a process involved in the formation, arrangement of constituent parts, or disassembly of an organelle. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0033046 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of sister chromatid segregation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of sister chromatid segregation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0033047 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of mitotic sister chromatid segregation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of sister chromatid segregation during mitosis. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0033273 </td>
   <td style="text-align:left;vertical-align:top;"> response to vitamin </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a vitamin stimulus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0034110 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of homotypic cell-cell adhesion </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate, or extent of homotypic cell-cell adhesion. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0034249 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of cellular amide metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of the chemical reactions and pathways involving amides. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0034330 </td>
   <td style="text-align:left;vertical-align:top;"> cell junction organization </td>
   <td style="text-align:left;vertical-align:top;"> A process that is carried out at the cellular level which results in the assembly, arrangement of constituent parts, or disassembly of a cell junction. A cell junction is a specialized region of connection between two cells or between a cell and the extracellular matrix. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0034332 </td>
   <td style="text-align:left;vertical-align:top;"> adherens junction organization </td>
   <td style="text-align:left;vertical-align:top;"> A process that is carried out at the cellular level which results in the assembly, arrangement of constituent parts, or disassembly of an adherens junction. An adherens junction is a cell-cell junction composed of the epithelial cadherin-catenin complex at which the cytoplasmic face of the plasma membrane is attached to actin filaments. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0034614 </td>
   <td style="text-align:left;vertical-align:top;"> cellular response to reactive oxygen species </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a reactive oxygen species stimulus. Reactive oxygen species include singlet oxygen, superoxide, and oxygen free radicals. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0034641 </td>
   <td style="text-align:left;vertical-align:top;"> cellular nitrogen compound metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving various organic and inorganic nitrogenous compounds, as carried out by individual cells. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0034655 </td>
   <td style="text-align:left;vertical-align:top;"> nucleobase-containing compound catabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways resulting in the breakdown of nucleobases, nucleosides, nucleotides and nucleic acids. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0034754 </td>
   <td style="text-align:left;vertical-align:top;"> cellular hormone metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving any hormone, naturally occurring substances secreted by specialized cells that affects the metabolism or behavior of other cells possessing functional receptors for the hormone, as carried out by individual cells. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0034762 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of transmembrane transport </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of the directed movement of a solute from one side of a membrane to the other. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0034767 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of ion transmembrane transport </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of the directed movement of ions from one side of a membrane to the other. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0035150 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of tube size </td>
   <td style="text-align:left;vertical-align:top;"> Ensuring that a tube is of the correct length and diameter. Tube size must be maintained not only during tube formation, but also throughout development and in some physiological processes. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0035265 </td>
   <td style="text-align:left;vertical-align:top;"> organ growth </td>
   <td style="text-align:left;vertical-align:top;"> The increase in size or mass of an organ. Organs are commonly observed as visibly distinct structures, but may also exist as loosely associated clusters of cells that function together as to perform a specific function. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0035295 </td>
   <td style="text-align:left;vertical-align:top;"> tube development </td>
   <td style="text-align:left;vertical-align:top;"> The process whose specific outcome is the progression of a tube over time, from its initial formation to a mature structure. Epithelial and endothelial tubes transport gases, liquids and cells from one site to another and form the basic structure of many organs and tissues including lung and trachea, kidney, the mammary gland, the vascular system and the gastrointestinal and urinary-genital tracts. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0035296 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of tube diameter </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the diameter of a tube. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0035556 </td>
   <td style="text-align:left;vertical-align:top;"> intracellular signal transduction </td>
   <td style="text-align:left;vertical-align:top;"> The process in which a signal is passed on to downstream components within the cell, which become activated themselves to further propagate the signal and finally trigger a change in the function or state of the cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0040007 </td>
   <td style="text-align:left;vertical-align:top;"> growth </td>
   <td style="text-align:left;vertical-align:top;"> The increase in size or mass of an entire organism, a part of an organism or a cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0040013 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of locomotion </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of locomotion of a cell or organism. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0040017 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of locomotion </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of locomotion of a cell or organism. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0042060 </td>
   <td style="text-align:left;vertical-align:top;"> wound healing </td>
   <td style="text-align:left;vertical-align:top;"> The series of events that restore integrity to a damaged tissue, following an injury. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0042063 </td>
   <td style="text-align:left;vertical-align:top;"> gliogenesis </td>
   <td style="text-align:left;vertical-align:top;"> The process that results in the generation of glial cells. This includes the production of glial progenitors and their differentiation into mature glia. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0042127 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of cell population proliferation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of cell proliferation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0042221 </td>
   <td style="text-align:left;vertical-align:top;"> response to chemical </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a chemical stimulus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0042310 </td>
   <td style="text-align:left;vertical-align:top;"> vasoconstriction </td>
   <td style="text-align:left;vertical-align:top;"> A decrease in the diameter of blood vessels, especially arteries, due to constriction of smooth muscle cells that line the vessels, and usually causing an increase in blood pressure. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0042325 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of phosphorylation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of addition of phosphate groups into a molecule. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0042327 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of phosphorylation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of addition of phosphate groups to a molecule. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0042330 </td>
   <td style="text-align:left;vertical-align:top;"> taxis </td>
   <td style="text-align:left;vertical-align:top;"> The directed movement of a motile cell or organism in response to an external stimulus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0042552 </td>
   <td style="text-align:left;vertical-align:top;"> myelination </td>
   <td style="text-align:left;vertical-align:top;"> The process in which myelin sheaths are formed and maintained around neurons. Oligodendrocytes in the brain and spinal cord and Schwann cells in the peripheral nervous system wrap axons with compact layers of their plasma membrane. Adjacent myelin segments are separated by a non-myelinated stretch of axon called a node of Ranvier. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0042592 </td>
   <td style="text-align:left;vertical-align:top;"> homeostatic process </td>
   <td style="text-align:left;vertical-align:top;"> Any biological process involved in the maintenance of an internal steady state. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0043066 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of apoptotic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of cell death by apoptotic process. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0043067 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of programmed cell death </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of programmed cell death, cell death resulting from activation of endogenous cellular processes. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0043170 </td>
   <td style="text-align:left;vertical-align:top;"> macromolecule metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving macromolecules, any molecule of high relative molecular mass, the structure of which essentially comprises the multiple repetition of units derived, actually or conceptually, from molecules of low relative molecular mass. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0043244 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of protein-containing complex disassembly </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of protein complex disassembly, the disaggregation of a protein complex into its constituent components. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0043254 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of protein-containing complex assembly </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of protein complex assembly. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0043277 </td>
   <td style="text-align:left;vertical-align:top;"> apoptotic cell clearance </td>
   <td style="text-align:left;vertical-align:top;"> The recognition and removal of an apoptotic cell by a neighboring cell or by a phagocyte. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0043280 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of cysteine-type endopeptidase activity involved in apoptotic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the activity of a cysteine-type endopeptidase involved in the apoptotic process. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0043281 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of cysteine-type endopeptidase activity involved in apoptotic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the activity of a cysteine-type endopeptidase involved in apoptosis. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0043406 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of MAP kinase activity </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of MAP kinase activity. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0043408 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of MAPK cascade </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of signal transduction mediated by the MAP kinase (MAPK) cascade. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0043412 </td>
   <td style="text-align:left;vertical-align:top;"> macromolecule modification </td>
   <td style="text-align:left;vertical-align:top;"> The covalent alteration of one or more monomeric units in a polypeptide, polynucleotide, polysaccharide, or other biological macromolecule, resulting in a change in its properties. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0043534 </td>
   <td style="text-align:left;vertical-align:top;"> blood vessel endothelial cell migration </td>
   <td style="text-align:left;vertical-align:top;"> The orderly movement of an endothelial cell into the extracellular matrix in order to form new blood vessels during angiogenesis. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0043549 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of kinase activity </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of kinase activity, the catalysis of the transfer of a phosphate group, usually from ATP, to a substrate molecule. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0044057 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of system process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of a system process, a multicellular organismal process carried out by any of the organs or tissues in an organ system. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0044085 </td>
   <td style="text-align:left;vertical-align:top;"> cellular component biogenesis </td>
   <td style="text-align:left;vertical-align:top;"> A process that results in the biosynthesis of constituent macromolecules, assembly, and arrangement of constituent parts of a cellular component. Includes biosynthesis of constituent macromolecules, and those macromolecular modifications that are involved in synthesis or assembly of the cellular component. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0044092 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of molecular function </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops or reduces the rate or extent of a molecular function, an elemental biological activity occurring at the molecular level, such as catalysis or binding. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0044093 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of molecular function </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the rate or extent of a molecular function, an elemental biological activity occurring at the molecular level, such as catalysis or binding. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0044237 </td>
   <td style="text-align:left;vertical-align:top;"> cellular metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways by which individual cells transform chemical substances. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0044238 </td>
   <td style="text-align:left;vertical-align:top;"> primary metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving those compounds which are formed as a part of the normal anabolic and catabolic processes. These processes take place in most, if not all, cells of the organism. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0044272 </td>
   <td style="text-align:left;vertical-align:top;"> sulfur compound biosynthetic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways resulting in the formation of compounds that contain sulfur, such as the amino acids methionine and cysteine or the tripeptide glutathione. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0044458 </td>
   <td style="text-align:left;vertical-align:top;"> motile cilium assembly </td>
   <td style="text-align:left;vertical-align:top;"> The aggregation, arrangement and bonding together of a set of components to form a motile cilium. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045185 </td>
   <td style="text-align:left;vertical-align:top;"> maintenance of protein location </td>
   <td style="text-align:left;vertical-align:top;"> Any process in which a protein is maintained in a location and prevented from moving elsewhere. These include sequestration, stabilization to prevent transport elsewhere and the active retrieval of proteins that do move away. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045444 </td>
   <td style="text-align:left;vertical-align:top;"> fat cell differentiation </td>
   <td style="text-align:left;vertical-align:top;"> The process in which a relatively unspecialized cell acquires specialized features of an adipocyte, an animal connective tissue cell specialized for the synthesis and storage of fat. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045446 </td>
   <td style="text-align:left;vertical-align:top;"> endothelial cell differentiation </td>
   <td style="text-align:left;vertical-align:top;"> The process in which a mesodermal, bone marrow or neural crest cell acquires specialized features of an endothelial cell, a thin flattened cell. A layer of such cells lines the inside surfaces of body cavities, blood vessels, and lymph vessels, making up the endothelium. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045471 </td>
   <td style="text-align:left;vertical-align:top;"> response to ethanol </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of an ethanol stimulus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045595 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of cell differentiation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of cell differentiation, the process in which relatively unspecialized cells acquire specialized structural and functional features. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045597 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of cell differentiation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of cell differentiation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045604 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of epidermal cell differentiation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of epidermal cell differentiation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045665 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of neuron differentiation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of neuron differentiation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045766 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of angiogenesis </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases angiogenesis. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045786 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of cell cycle </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents or reduces the rate or extent of progression through the cell cycle. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045807 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of endocytosis </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of endocytosis. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045834 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of lipid metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of the chemical reactions and pathways involving lipids. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045839 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of mitotic nuclear division </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents or reduces the rate or extent of mitosis. Mitosis is the division of the eukaryotic cell nucleus to produce two daughter nuclei that, usually, contain the identical chromosome complement to their mother. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045859 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of protein kinase activity </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of protein kinase activity. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045860 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of protein kinase activity </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of protein kinase activity. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045861 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of proteolysis </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of the hydrolysis of a peptide bond or bonds within a protein. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045927 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of growth </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the rate or extent of growth, the increase in size or mass of all or part of an organism. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045931 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of mitotic cell cycle </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the rate or extent of progression through the mitotic cell cycle. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045936 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of phosphate metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of the chemical reactions and pathways involving phosphates. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045944 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of transcription by RNA polymerase II </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of transcription from an RNA polymerase II promoter. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0046483 </td>
   <td style="text-align:left;vertical-align:top;"> heterocycle metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving heterocyclic compounds, those with a cyclic molecular structure and at least two different atoms in the ring (or rings). </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0046890 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of lipid biosynthetic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of the chemical reactions and pathways resulting in the formation of lipids. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048015 </td>
   <td style="text-align:left;vertical-align:top;"> phosphatidylinositol-mediated signaling </td>
   <td style="text-align:left;vertical-align:top;"> The series of molecular signals in which a cell uses a phosphatidylinositol-mediated signaling to convert a signal into a response. Phosphatidylinositols include phosphatidylinositol (PtdIns) and its phosphorylated derivatives. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048144 </td>
   <td style="text-align:left;vertical-align:top;"> fibroblast proliferation </td>
   <td style="text-align:left;vertical-align:top;"> The multiplication or reproduction of fibroblast cells, resulting in the expansion of the fibroblast population. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048514 </td>
   <td style="text-align:left;vertical-align:top;"> blood vessel morphogenesis </td>
   <td style="text-align:left;vertical-align:top;"> The process in which the anatomical structures of blood vessels are generated and organized. The blood vessel is the vasculature carrying blood. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048518 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of biological process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of a biological process. Biological processes are regulated by many means; examples include the control of gene expression, protein modification or interaction with a protein or substrate molecule. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048585 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of response to stimulus </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of a response to a stimulus. Response to stimulus is a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a stimulus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048588 </td>
   <td style="text-align:left;vertical-align:top;"> developmental cell growth </td>
   <td style="text-align:left;vertical-align:top;"> The growth of a cell, where growth contributes to the progression of the cell over time from one condition to another. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048589 </td>
   <td style="text-align:left;vertical-align:top;"> developmental growth </td>
   <td style="text-align:left;vertical-align:top;"> The increase in size or mass of an entire organism, a part of an organism or a cell, where the increase in size or mass has the specific outcome of the progression of the organism over time from one condition to another. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048638 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of developmental growth </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of developmental growth. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048646 </td>
   <td style="text-align:left;vertical-align:top;"> anatomical structure formation involved in morphogenesis </td>
   <td style="text-align:left;vertical-align:top;"> The developmental process pertaining to the initial formation of an anatomical structure from unspecified parts. This process begins with the specific processes that contribute to the appearance of the discrete structure and ends when the structural rudiment is recognizable. An anatomical structure is any biological entity that occupies space and is distinguished from its surroundings. Anatomical structures can be macroscopic such as a carpel, or microscopic such as an acrosome. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048667 </td>
   <td style="text-align:left;vertical-align:top;"> cell morphogenesis involved in neuron differentiation </td>
   <td style="text-align:left;vertical-align:top;"> The process in which the structures of a neuron are generated and organized. This process occurs while the initially relatively unspecialized cell is acquiring the specialized features of a neuron. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048675 </td>
   <td style="text-align:left;vertical-align:top;"> axon extension </td>
   <td style="text-align:left;vertical-align:top;"> Long distance growth of a single axon process involved in cellular development. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048706 </td>
   <td style="text-align:left;vertical-align:top;"> embryonic skeletal system development </td>
   <td style="text-align:left;vertical-align:top;"> The process, occurring during the embryonic phase, whose specific outcome is the progression of the skeleton over time, from its formation to the mature structure. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048729 </td>
   <td style="text-align:left;vertical-align:top;"> tissue morphogenesis </td>
   <td style="text-align:left;vertical-align:top;"> The process in which the anatomical structures of a tissue are generated and organized. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048732 </td>
   <td style="text-align:left;vertical-align:top;"> gland development </td>
   <td style="text-align:left;vertical-align:top;"> The process whose specific outcome is the progression of a gland over time, from its formation to the mature structure. A gland is an organ specialised for secretion. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048762 </td>
   <td style="text-align:left;vertical-align:top;"> mesenchymal cell differentiation </td>
   <td style="text-align:left;vertical-align:top;"> The process in which a relatively unspecialized cell acquires specialized features of a mesenchymal cell. A mesenchymal cell is a loosely associated cell that is part of the connective tissue in an organism. Mesenchymal cells give rise to more mature connective tissue cell types. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048812 </td>
   <td style="text-align:left;vertical-align:top;"> neuron projection morphogenesis </td>
   <td style="text-align:left;vertical-align:top;"> The process in which the anatomical structures of a neuron projection are generated and organized. A neuron projection is any process extending from a neural cell, such as axons or dendrites. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048839 </td>
   <td style="text-align:left;vertical-align:top;"> inner ear development </td>
   <td style="text-align:left;vertical-align:top;"> The process whose specific outcome is the progression of the inner ear over time, from its formation to the mature structure. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048871 </td>
   <td style="text-align:left;vertical-align:top;"> multicellular organismal homeostasis </td>
   <td style="text-align:left;vertical-align:top;"> Any process involved in the maintenance of an internal steady state at the level of the multicellular organism. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0050679 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of epithelial cell proliferation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the rate or extent of epithelial cell proliferation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0050684 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of mRNA processing </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of mRNA processing, those processes involved in the conversion of a primary mRNA transcript into a mature mRNA prior to its translation into polypeptide. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0050730 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of peptidyl-tyrosine phosphorylation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of the phosphorylation of peptidyl-tyrosine. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0050770 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of axonogenesis </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of axonogenesis, the generation of an axon, the long process of a neuron. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0050789 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of biological process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of a biological process. Biological processes are regulated by many means; examples include the control of gene expression, protein modification or interaction with a protein or substrate molecule. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0050790 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of catalytic activity </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the activity of an enzyme. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0050793 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of developmental process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of development, the biological process whose specific outcome is the progression of a multicellular organism over time from an initial condition (e.g. a zygote, or a young adult) to a later condition (e.g. a multicellular animal or an aged adult). </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0050817 </td>
   <td style="text-align:left;vertical-align:top;"> coagulation </td>
   <td style="text-align:left;vertical-align:top;"> The process in which a fluid solution, or part of it, changes into a solid or semisolid mass. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0050818 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of coagulation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of coagulation, the process in which a fluid solution, or part of it, changes into a solid or semisolid mass. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0050848 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of calcium-mediated signaling </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of calcium-mediated signaling, the process in which a cell uses calcium ions to convert an extracellular signal into a response. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0050865 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of cell activation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of cell activation, the change in the morphology or behavior of a cell resulting from exposure to an activating factor such as a cellular or soluble ligand. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0050866 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of cell activation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of cell activation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0050867 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of cell activation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of activation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0050890 </td>
   <td style="text-align:left;vertical-align:top;"> cognition </td>
   <td style="text-align:left;vertical-align:top;"> The operation of the mind by which an organism becomes aware of objects of thought or perception; it includes the mental activities associated with thinking, learning, and memory. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0050896 </td>
   <td style="text-align:left;vertical-align:top;"> response to stimulus </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a stimulus. The process begins with detection of the stimulus and ends with a change in state or activity or the cell or organism. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0050920 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of chemotaxis </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of the directed movement of a motile cell or organism in response to a specific chemical concentration gradient. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051028 </td>
   <td style="text-align:left;vertical-align:top;"> mRNA transport </td>
   <td style="text-align:left;vertical-align:top;"> The directed movement of mRNA, messenger ribonucleic acid, into, out of or within a cell, or between cells, by means of some agent such as a transporter or pore. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051049 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of transport </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of the directed movement of substances (such as macromolecules, small molecules, ions) into, out of or within a cell, or between cells, by means of some agent such as a transporter or pore. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051093 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of developmental process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents or reduces the rate or extent of development, the biological process whose specific outcome is the progression of an organism over time from an initial condition (e.g. a zygote, or a young adult) to a later condition (e.g. a multicellular animal or an aged adult). </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051094 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of developmental process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the rate or extent of development, the biological process whose specific outcome is the progression of an organism over time from an initial condition (e.g. a zygote, or a young adult) to a later condition (e.g. a multicellular animal or an aged adult). </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051098 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of binding </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of binding, the selective interaction of a molecule with one or more specific sites on another molecule. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051129 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of cellular component organization </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of a process involved in the formation, arrangement of constituent parts, or disassembly of cell structures, including the plasma membrane and any external encapsulating structures such as the cell wall and cell envelope. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051173 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of nitrogen compound metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of the chemical reactions and pathways involving nitrogen or nitrogenous compounds. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051174 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of phosphorus metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of the chemical reactions and pathways involving phosphorus or compounds containing phosphorus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051209 </td>
   <td style="text-align:left;vertical-align:top;"> release of sequestered calcium ion into cytosol </td>
   <td style="text-align:left;vertical-align:top;"> The process in which calcium ions sequestered in the endoplasmic reticulum, Golgi apparatus or mitochondria are released into the cytosolic compartment. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051216 </td>
   <td style="text-align:left;vertical-align:top;"> cartilage development </td>
   <td style="text-align:left;vertical-align:top;"> The process whose specific outcome is the progression of a cartilage element over time, from its formation to the mature structure. Cartilage elements are skeletal elements that consist of connective tissue dominated by extracellular matrix containing collagen type II and large amounts of proteoglycan, particularly chondroitin sulfate. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051239 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of multicellular organismal process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of a multicellular organismal process, the processes pertinent to the function of a multicellular organism above the cellular level; includes the integrated processes of tissues and organs. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051240 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of multicellular organismal process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of an organismal process, any of the processes pertinent to the function of an organism above the cellular level; includes the integrated processes of tissues and organs. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051282 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of sequestering of calcium ion </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of the binding or confining calcium ions such that they are separated from other components of a biological system. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051302 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of cell division </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of the physical partitioning and separation of a cell into daughter cells. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051346 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of hydrolase activity </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops or reduces the rate of hydrolase activity, the catalysis of the hydrolysis of various bonds. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051347 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of transferase activity </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of transferase activity, the catalysis of the transfer of a group, e.g. a methyl group, glycosyl group, acyl group, phosphorus-containing, or other groups, from a donor compound to an acceptor. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051402 </td>
   <td style="text-align:left;vertical-align:top;"> neuron apoptotic process </td>
   <td style="text-align:left;vertical-align:top;"> Any apoptotic process in a neuron, the basic cellular unit of nervous tissue. Each neuron consists of a body, an axon, and dendrites. Their purpose is to receive, conduct, and transmit impulses in the nervous system. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051493 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of cytoskeleton organization </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of the formation, arrangement of constituent parts, or disassembly of cytoskeletal structures. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051604 </td>
   <td style="text-align:left;vertical-align:top;"> protein maturation </td>
   <td style="text-align:left;vertical-align:top;"> Any process leading to the attainment of the full functional capacity of a protein. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051668 </td>
   <td style="text-align:left;vertical-align:top;"> localization within membrane </td>
   <td style="text-align:left;vertical-align:top;"> Any process in which a substance or cellular entity, such as a protein complex or organelle, is transported to, and/or maintained in, a specific location within a membrane. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051716 </td>
   <td style="text-align:left;vertical-align:top;"> cellular response to stimulus </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a stimulus. The process begins with detection of the stimulus by a cell and ends with a change in state or activity or the cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051897 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of protein kinase B signaling </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of protein kinase B signaling, a series of reactions mediated by the intracellular serine/threonine kinase protein kinase B. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051928 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of calcium ion transport </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of the directed movement of calcium ions into, out of or within a cell, or between cells, by means of some agent such as a transporter or pore. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0052547 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of peptidase activity </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of peptidase activity, the hydrolysis of peptide bonds within proteins. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0055001 </td>
   <td style="text-align:left;vertical-align:top;"> muscle cell development </td>
   <td style="text-align:left;vertical-align:top;"> The process whose specific outcome is the progression of a muscle cell over time, from its formation to the mature structure. Muscle cell development does not include the steps involved in committing an unspecified cell to the muscle cell fate. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0055076 </td>
   <td style="text-align:left;vertical-align:top;"> transition metal ion homeostasis </td>
   <td style="text-align:left;vertical-align:top;"> Any process involved in the maintenance of an internal steady state of transition metal ions within an organism or cell. A transition metal is an element whose atom has an incomplete d-subshell of extranuclear electrons, or which gives rise to a cation or cations with an incomplete d-subshell. Transition metals often have more than one valency state. Biologically relevant transition metals include vanadium, manganese, iron, copper, cobalt, nickel, molybdenum and silver. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0055082 </td>
   <td style="text-align:left;vertical-align:top;"> cellular chemical homeostasis </td>
   <td style="text-align:left;vertical-align:top;"> Any biological process involved in the maintenance of an internal steady state of a chemical at the level of the cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0060021 </td>
   <td style="text-align:left;vertical-align:top;"> roof of mouth development </td>
   <td style="text-align:left;vertical-align:top;"> The biological process whose specific outcome is the progression of the roof of the mouth from an initial condition to its mature state. This process begins with the formation of the structure and ends with the mature structure. The roof of the mouth is the partition that separates the nasal and oral cavities. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0060047 </td>
   <td style="text-align:left;vertical-align:top;"> heart contraction </td>
   <td style="text-align:left;vertical-align:top;"> The multicellular organismal process in which the heart decreases in volume in a characteristic way to propel blood through the body. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0060249 </td>
   <td style="text-align:left;vertical-align:top;"> anatomical structure homeostasis </td>
   <td style="text-align:left;vertical-align:top;"> A homeostatic process involved in the maintenance of an internal steady state within a defined anatomical structure of an organism, including control of cellular proliferation and death and control of metabolic function. An anatomical structure is any biological entity that occupies space and is distinguished from its surroundings. Anatomical structures can be macroscopic such as a carpel, or microscopic such as an acrosome. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0060255 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of macromolecule metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of the chemical reactions and pathways involving macromolecules, any molecule of high relative molecular mass, the structure of which essentially comprises the multiple repetition of units derived, actually or conceptually, from molecules of low relative molecular mass. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0060314 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of ryanodine-sensitive calcium-release channel activity </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the activity of a ryanodine-sensitive calcium-release channel. The ryanodine-sensitive calcium-release channel catalyzes the transmembrane transfer of a calcium ion by a channel that opens when a ryanodine class ligand has been bound by the channel complex or one of its constituent parts. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0060322 </td>
   <td style="text-align:left;vertical-align:top;"> head development </td>
   <td style="text-align:left;vertical-align:top;"> The biological process whose specific outcome is the progression of a head from an initial condition to its mature state. The head is the anterior-most division of the body. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0060348 </td>
   <td style="text-align:left;vertical-align:top;"> bone development </td>
   <td style="text-align:left;vertical-align:top;"> The process whose specific outcome is the progression of bone over time, from its formation to the mature structure. Bone is the hard skeletal connective tissue consisting of both mineral and cellular components. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0060548 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of cell death </td>
   <td style="text-align:left;vertical-align:top;"> Any process that decreases the rate or frequency of cell death. Cell death is the specific activation or halting of processes within a cell so that its vital functions markedly cease, rather than simply deteriorating gradually over time, which culminates in cell death. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0060560 </td>
   <td style="text-align:left;vertical-align:top;"> developmental growth involved in morphogenesis </td>
   <td style="text-align:left;vertical-align:top;"> The increase in size or mass of an anatomical structure that contributes to the structure attaining its shape. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0061025 </td>
   <td style="text-align:left;vertical-align:top;"> membrane fusion </td>
   <td style="text-align:left;vertical-align:top;"> The membrane organization process that joins two lipid bilayers to form a single membrane. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0061028 </td>
   <td style="text-align:left;vertical-align:top;"> establishment of endothelial barrier </td>
   <td style="text-align:left;vertical-align:top;"> The establishment of a barrier between endothelial cell layers, such as those in the brain, lung or intestine, to exert specific and selective control over the passage of water and solutes, thus allowing formation and maintenance of compartments that differ in fluid and solute composition. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0061041 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of wound healing </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the rate, frequency, or extent of the series of events that restore integrity to a damaged tissue, following an injury. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0061061 </td>
   <td style="text-align:left;vertical-align:top;"> muscle structure development </td>
   <td style="text-align:left;vertical-align:top;"> The progression of a muscle structure over time, from its formation to its mature state. Muscle structures are contractile cells, tissues or organs that are found in multicellular organisms. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0061180 </td>
   <td style="text-align:left;vertical-align:top;"> mammary gland epithelium development </td>
   <td style="text-align:left;vertical-align:top;"> The process whose specific outcome is the progression of the mammary gland epithelium over time, from its formation to the mature structure. The mammary gland is a large compound sebaceous gland that in female mammals is modified to secrete milk. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0061337 </td>
   <td style="text-align:left;vertical-align:top;"> cardiac conduction </td>
   <td style="text-align:left;vertical-align:top;"> Transfer of an organized electrical impulse across the heart to coordinate the contraction of cardiac muscles. The process begins with generation of an action potential (in the sinoatrial node (SA) in humans) and ends with a change in the rate, frequency, or extent of the contraction of the heart muscles. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0061387 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of extent of cell growth </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the extent of cell growth. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0061564 </td>
   <td style="text-align:left;vertical-align:top;"> axon development </td>
   <td style="text-align:left;vertical-align:top;"> The progression of an axon over time. Covers axonogenesis (de novo generation of an axon) and axon regeneration (regrowth), as well as processes pertaining to the progression of the axon over time (fasciculation and defasciculation). </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0061572 </td>
   <td style="text-align:left;vertical-align:top;"> actin filament bundle organization </td>
   <td style="text-align:left;vertical-align:top;"> A process that results in the assembly, arrangement of constituent parts, or disassembly of an actin filament bundle. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0062197 </td>
   <td style="text-align:left;vertical-align:top;"> cellular response to chemical stress </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a chemical stimulus indicating the organism is under stress. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0065003 </td>
   <td style="text-align:left;vertical-align:top;"> protein-containing complex assembly </td>
   <td style="text-align:left;vertical-align:top;"> The aggregation, arrangement and bonding together of a set of macromolecules to form a protein-containing complex. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0065007 </td>
   <td style="text-align:left;vertical-align:top;"> biological regulation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates a measurable attribute of any biological process, quality or function. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0065008 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of biological quality </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates a qualitative or quantitative trait of a biological quality. A biological quality is a measurable attribute of an organism or part of an organism, such as size, mass, shape, color, etc. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0070252 </td>
   <td style="text-align:left;vertical-align:top;"> actin-mediated cell contraction </td>
   <td style="text-align:left;vertical-align:top;"> The actin filament-based process in which cytoplasmic actin filaments slide past one another resulting in contraction of all or part of the cell body. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0070304 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of stress-activated protein kinase signaling cascade </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of signaling via the stress-activated protein kinase signaling cascade. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0070374 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of ERK1 and ERK2 cascade </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of signal transduction mediated by the ERK1 and ERK2 cascade. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0070482 </td>
   <td style="text-align:left;vertical-align:top;"> response to oxygen levels </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a stimulus reflecting the presence, absence, or concentration of oxygen. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0070585 </td>
   <td style="text-align:left;vertical-align:top;"> protein localization to mitochondrion </td>
   <td style="text-align:left;vertical-align:top;"> A process in which a protein is transported to, or maintained in, a location within the mitochondrion. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0070661 </td>
   <td style="text-align:left;vertical-align:top;"> leukocyte proliferation </td>
   <td style="text-align:left;vertical-align:top;"> The expansion of a leukocyte population by cell division. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0070665 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of leukocyte proliferation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of leukocyte proliferation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0071363 </td>
   <td style="text-align:left;vertical-align:top;"> cellular response to growth factor stimulus </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a growth factor stimulus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0071375 </td>
   <td style="text-align:left;vertical-align:top;"> cellular response to peptide hormone stimulus </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a peptide hormone stimulus. A peptide hormone is any of a class of peptides that are secreted into the blood stream and have endocrine functions in living animals. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0071396 </td>
   <td style="text-align:left;vertical-align:top;"> cellular response to lipid </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a lipid stimulus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0071407 </td>
   <td style="text-align:left;vertical-align:top;"> cellular response to organic cyclic compound </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of an organic cyclic compound stimulus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0071417 </td>
   <td style="text-align:left;vertical-align:top;"> cellular response to organonitrogen compound </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of an organonitrogen stimulus. An organonitrogen compound is formally a compound containing at least one carbon-nitrogen bond. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0071496 </td>
   <td style="text-align:left;vertical-align:top;"> cellular response to external stimulus </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of an external stimulus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0071559 </td>
   <td style="text-align:left;vertical-align:top;"> response to transforming growth factor beta </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a transforming growth factor beta stimulus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0071704 </td>
   <td style="text-align:left;vertical-align:top;"> organic substance metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving an organic substance, any molecular entity containing carbon. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0071805 </td>
   <td style="text-align:left;vertical-align:top;"> potassium ion transmembrane transport </td>
   <td style="text-align:left;vertical-align:top;"> A process in which a potassium ion is transported from one side of a membrane to the other. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0071840 </td>
   <td style="text-align:left;vertical-align:top;"> cellular component organization or biogenesis </td>
   <td style="text-align:left;vertical-align:top;"> A process that results in the biosynthesis of constituent macromolecules, assembly, arrangement of constituent parts, or disassembly of a cellular component. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0071900 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of protein serine/threonine kinase activity </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the rate, frequency, or extent of protein serine/threonine kinase activity. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0071901 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of protein serine/threonine kinase activity </td>
   <td style="text-align:left;vertical-align:top;"> Any process that decreases the rate, frequency, or extent of protein serine/threonine kinase activity. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0072330 </td>
   <td style="text-align:left;vertical-align:top;"> monocarboxylic acid biosynthetic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways resulting in the formation of monocarboxylic acids, any organic acid containing one carboxyl (-COOH) group. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0072359 </td>
   <td style="text-align:left;vertical-align:top;"> circulatory system development </td>
   <td style="text-align:left;vertical-align:top;"> The process whose specific outcome is the progression of the circulatory system over time, from its formation to the mature structure. The circulatory system is the organ system that passes nutrients (such as amino acids and electrolytes), gases, hormones, blood cells, etc. to and from cells in the body to help fight diseases and help stabilize body temperature and pH to maintain homeostasis. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0072507 </td>
   <td style="text-align:left;vertical-align:top;"> divalent inorganic cation homeostasis </td>
   <td style="text-align:left;vertical-align:top;"> Any process involved in the maintenance of an internal steady state of divalent cations within an organism or cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0072522 </td>
   <td style="text-align:left;vertical-align:top;"> purine-containing compound biosynthetic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways resulting in the formation of a purine-containing compound, i.e. any compound that contains purine or a formal derivative thereof. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0072594 </td>
   <td style="text-align:left;vertical-align:top;"> establishment of protein localization to organelle </td>
   <td style="text-align:left;vertical-align:top;"> The directed movement of a protein to a specific location on or in an organelle. Encompasses establishment of localization in the membrane or lumen of a membrane-bounded organelle. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0072659 </td>
   <td style="text-align:left;vertical-align:top;"> protein localization to plasma membrane </td>
   <td style="text-align:left;vertical-align:top;"> A process in which a protein is transported to, or maintained in, a specific location in the plasma membrane. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0080134 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of response to stress </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of a response to stress. Response to stress is a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a disturbance in organismal or cellular homeostasis, usually, but not necessarily, exogenous (e.g. temperature, humidity, ionizing radiation). </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0090066 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of anatomical structure size </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the size of an anatomical structure. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0090100 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of transmembrane receptor protein serine/threonine kinase signaling pathway </td>
   <td style="text-align:left;vertical-align:top;"> Any process that increases the rate, frequency, or extent of the series of molecular signals generated as a consequence of a transmembrane receptor serine/threonine kinase binding to its physiological ligand. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0090304 </td>
   <td style="text-align:left;vertical-align:top;"> nucleic acid metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> Any cellular metabolic process involving nucleic acids. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0090502 </td>
   <td style="text-align:left;vertical-align:top;"> RNA phosphodiester bond hydrolysis, endonucleolytic </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving the hydrolysis of internal 3',5'-phosphodiester bonds in one or two strands of ribonucleotides. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0097237 </td>
   <td style="text-align:left;vertical-align:top;"> cellular response to toxic substance </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a toxic stimulus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0097746 </td>
   <td style="text-align:left;vertical-align:top;"> blood vessel diameter maintenance </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the diameter of blood vessels. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0098609 </td>
   <td style="text-align:left;vertical-align:top;"> cell-cell adhesion </td>
   <td style="text-align:left;vertical-align:top;"> The attachment of one cell to another cell via adhesion molecules. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0098754 </td>
   <td style="text-align:left;vertical-align:top;"> detoxification </td>
   <td style="text-align:left;vertical-align:top;"> Any process that reduces or removes the toxicity of a toxic substance. These may include transport of the toxic substance away from sensitive areas and to compartments or complexes whose purpose is sequestration of the toxic substance. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0098869 </td>
   <td style="text-align:left;vertical-align:top;"> cellular oxidant detoxification </td>
   <td style="text-align:left;vertical-align:top;"> Any process carried out at the cellular level that reduces or removes the toxicity superoxide radicals or hydrogen peroxide. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0098876 </td>
   <td style="text-align:left;vertical-align:top;"> vesicle-mediated transport to the plasma membrane </td>
   <td style="text-align:left;vertical-align:top;"> The directed movement of substances to the plasma membrane in transport vesicles that fuse with the plasma membrane by exocytosis. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0110053 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of actin filament organization </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of actin filament organization. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0120036 </td>
   <td style="text-align:left;vertical-align:top;"> plasma membrane bounded cell projection organization </td>
   <td style="text-align:left;vertical-align:top;"> A process that is carried out at the cellular level which results in the assembly, arrangement of constituent parts, or disassembly of a plasma membrane bounded prolongation or process extending from a cell, e.g. a cilium or axon. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1900046 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of hemostasis </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1901342 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of vasculature development </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1901343 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of vasculature development </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1901360 </td>
   <td style="text-align:left;vertical-align:top;"> organic cyclic compound metabolic process </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1901564 </td>
   <td style="text-align:left;vertical-align:top;"> organonitrogen compound metabolic process </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1901570 </td>
   <td style="text-align:left;vertical-align:top;"> fatty acid derivative biosynthetic process </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1901653 </td>
   <td style="text-align:left;vertical-align:top;"> cellular response to peptide </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1901698 </td>
   <td style="text-align:left;vertical-align:top;"> response to nitrogen compound </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1901701 </td>
   <td style="text-align:left;vertical-align:top;"> cellular response to oxygen-containing compound </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1901880 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of protein depolymerization </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1901888 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of cell junction assembly </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1901991 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of mitotic cell cycle phase transition </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1902680 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of RNA biosynthetic process </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1902749 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of cell cycle G2/M phase transition </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1902905 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of supramolecular fiber organization </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1903008 </td>
   <td style="text-align:left;vertical-align:top;"> organelle disassembly </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1903035 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of response to wounding </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1903036 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of response to wounding </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1903037 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of leukocyte cell-cell adhesion </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1903039 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of leukocyte cell-cell adhesion </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1903078 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of protein localization to plasma membrane </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1903169 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of calcium ion transmembrane transport </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1903317 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of protein maturation </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1903508 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of nucleic acid-templated transcription </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1903844 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of cellular response to transforming growth factor beta stimulus </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1904062 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of cation transmembrane transport </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1904064 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of cation transmembrane transport </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1905475 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of protein localization to membrane </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1990000 </td>
   <td style="text-align:left;vertical-align:top;"> amyloid fibril formation </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:2000026 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of multicellular organismal development </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:2000116 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of cysteine-type endopeptidase activity </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:2000145 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of cell motility </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:2000816 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of mitotic sister chromatid separation </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:2001056 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of cysteine-type endopeptidase activity </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:2001233 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of apoptotic signaling pathway </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:2001243 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of intrinsic apoptotic signaling pathway </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:2001258 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of cation channel activity </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:2001259 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of cation channel activity </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d., Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006812 </td>
   <td style="text-align:left;vertical-align:top;"> cation transport </td>
   <td style="text-align:left;vertical-align:top;"> The directed movement of cations, atoms or small molecules with a net positive charge, into, out of or within a cell, or between cells, by means of some agent such as a transporter or pore. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d., Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0008015 </td>
   <td style="text-align:left;vertical-align:top;"> blood circulation </td>
   <td style="text-align:left;vertical-align:top;"> The flow of blood through the body of an animal, enabling the transport of nutrients to the tissues and the removal of waste products. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d., Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010563 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of phosphorus metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that decreases the frequency, rate or extent of the chemical reactions and pathways involving phosphorus or compounds containing phosphorus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d., Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010959 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of metal ion transport </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate, or extent of metal ion transport. Metal ion transport is the directed movement of metal ions, any metal ion with an electric charge, into, out of or within a cell, or between cells, by means of some agent such as a transporter or pore. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d., Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0030003 </td>
   <td style="text-align:left;vertical-align:top;"> cellular cation homeostasis </td>
   <td style="text-align:left;vertical-align:top;"> Any process involved in the maintenance of an internal steady state of cations at the level of a cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d., Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048523 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of cellular process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of a cellular process, any of those that are carried out at the cellular level, but are not necessarily restricted to a single cell. For example, cell communication occurs among more than one cell, but occurs at the cellular level. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d., Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0050801 </td>
   <td style="text-align:left;vertical-align:top;"> ion homeostasis </td>
   <td style="text-align:left;vertical-align:top;"> Any process involved in the maintenance of an internal steady state of ions within an organism or cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d., Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0055065 </td>
   <td style="text-align:left;vertical-align:top;"> metal ion homeostasis </td>
   <td style="text-align:left;vertical-align:top;"> Any process involved in the maintenance of an internal steady state of metal ions within an organism or cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> 3 d., Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0098662 </td>
   <td style="text-align:left;vertical-align:top;"> inorganic cation transmembrane transport </td>
   <td style="text-align:left;vertical-align:top;"> A process in which an inorganic cation is transported from one side of a membrane to the other by means of some agent such as a transporter or pore. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0001525 </td>
   <td style="text-align:left;vertical-align:top;"> angiogenesis </td>
   <td style="text-align:left;vertical-align:top;"> Blood vessel formation when new vessels emerge from the proliferation of pre-existing blood vessels. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0002042 </td>
   <td style="text-align:left;vertical-align:top;"> cell migration involved in sprouting angiogenesis </td>
   <td style="text-align:left;vertical-align:top;"> The orderly movement of endothelial cells into the extracellular matrix in order to form new blood vessels involved in sprouting angiogenesis. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006637 </td>
   <td style="text-align:left;vertical-align:top;"> acyl-CoA metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving acyl-CoA, any derivative of coenzyme A in which the sulfhydryl group is in thiolester linkage with an acyl group. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006779 </td>
   <td style="text-align:left;vertical-align:top;"> porphyrin-containing compound biosynthetic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways resulting in the formation of any member of a large group of derivatives or analogs of porphyrin. Porphyrin consists of a ring of four pyrrole nuclei linked each to the next at their alpha positions through a methine group. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006783 </td>
   <td style="text-align:left;vertical-align:top;"> heme biosynthetic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways resulting in the formation of heme, any compound of iron complexed in a porphyrin (tetrapyrrole) ring, from less complex precursors. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006790 </td>
   <td style="text-align:left;vertical-align:top;"> sulfur compound metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving the nonmetallic element sulfur or compounds that contain sulfur, such as the amino acids methionine and cysteine or the tripeptide glutathione. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006816 </td>
   <td style="text-align:left;vertical-align:top;"> calcium ion transport </td>
   <td style="text-align:left;vertical-align:top;"> The directed movement of calcium (Ca) ions into, out of or within a cell, or between cells, by means of some agent such as a transporter or pore. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006919 </td>
   <td style="text-align:left;vertical-align:top;"> activation of cysteine-type endopeptidase activity involved in apoptotic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that initiates the activity of the inactive enzyme cysteine-type endopeptidase in the context of an apoptotic process. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007015 </td>
   <td style="text-align:left;vertical-align:top;"> actin filament organization </td>
   <td style="text-align:left;vertical-align:top;"> A process that is carried out at the cellular level which results in the assembly, arrangement of constituent parts, or disassembly of cytoskeletal structures comprising actin filaments. Includes processes that control the spatial distribution of actin filaments, such as organizing filaments into meshworks, bundles, or other structures, as by cross-linking. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007163 </td>
   <td style="text-align:left;vertical-align:top;"> establishment or maintenance of cell polarity </td>
   <td style="text-align:left;vertical-align:top;"> Any cellular process that results in the specification, formation or maintenance of anisotropic intracellular organization or cell growth patterns. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0009719 </td>
   <td style="text-align:left;vertical-align:top;"> response to endogenous stimulus </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a stimulus arising within the organism. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0009890 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of biosynthetic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the rate of the chemical reactions and pathways resulting in the formation of substances. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0009895 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of catabolic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of the chemical reactions and pathways resulting in the breakdown of substances. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010631 </td>
   <td style="text-align:left;vertical-align:top;"> epithelial cell migration </td>
   <td style="text-align:left;vertical-align:top;"> The orderly movement of an epithelial cell from one site to another, often during the development of a multicellular organism. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010770 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of cell morphogenesis involved in differentiation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that increases the frequency, rate or extent of cell morphogenesis contributing to cell differentiation. Cell morphogenesis involved in differentiation is the change in form (cell shape and size) that occurs when relatively unspecialized cells acquire specialized structural and/or functional features that characterize the cells, tissues, or organs of the mature organism or some other relatively stable phase of the organism's life history. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0014075 </td>
   <td style="text-align:left;vertical-align:top;"> response to amine </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of an amine stimulus. An amine is a compound formally derived from ammonia by replacing one, two or three hydrogen atoms by hydrocarbyl groups. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0019098 </td>
   <td style="text-align:left;vertical-align:top;"> reproductive behavior </td>
   <td style="text-align:left;vertical-align:top;"> The specific behavior of an organism that is associated with reproduction. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0030029 </td>
   <td style="text-align:left;vertical-align:top;"> actin filament-based process </td>
   <td style="text-align:left;vertical-align:top;"> Any cellular process that depends upon or alters the actin cytoskeleton, that part of the cytoskeleton comprising actin filaments and their associated proteins. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0030041 </td>
   <td style="text-align:left;vertical-align:top;"> actin filament polymerization </td>
   <td style="text-align:left;vertical-align:top;"> Assembly of actin filaments by the addition of actin monomers to a filament. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0030511 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of transforming growth factor beta receptor signaling pathway </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of TGF-beta receptor signaling pathway activity. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0031099 </td>
   <td style="text-align:left;vertical-align:top;"> regeneration </td>
   <td style="text-align:left;vertical-align:top;"> The regrowth of a lost or destroyed body part, such as an organ or tissue. This process may occur via renewal, repair, and/or growth alone (i.e. increase in size or mass). </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0033013 </td>
   <td style="text-align:left;vertical-align:top;"> tetrapyrrole metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving tetrapyrroles, natural pigments containing four pyrrole rings joined by one-carbon units linking position 2 of one pyrrole ring to position 5 of the next. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0033014 </td>
   <td style="text-align:left;vertical-align:top;"> tetrapyrrole biosynthetic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways leading to the formation of tetrapyrroles, natural pigments containing four pyrrole rings joined by one-carbon units linking position 2 of one pyrrole ring to position 5 of the next. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0033598 </td>
   <td style="text-align:left;vertical-align:top;"> mammary gland epithelial cell proliferation </td>
   <td style="text-align:left;vertical-align:top;"> The multiplication or reproduction of mammary gland epithelial cells, resulting in the expansion of a cell population. Mammary gland epithelial cells make up the covering of surfaces of the mammary gland. The mammary gland is a large compound sebaceous gland that in female mammals is modified to secrete milk. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0033628 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of cell adhesion mediated by integrin </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate, or extent of cell adhesion mediated by integrin. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0033865 </td>
   <td style="text-align:left;vertical-align:top;"> nucleoside bisphosphate metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving a nucleoside bisphosphate, a compound consisting of a nucleobase linked to a deoxyribose or ribose sugar esterified with one phosphate group attached to each of two different hydroxyl groups on the sugar. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0034030 </td>
   <td style="text-align:left;vertical-align:top;"> ribonucleoside bisphosphate biosynthetic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways resulting in the formation of a ribonucleoside bisphosphate, a compound consisting of a nucleobase linked to a ribose sugar esterified with one phosphate group attached to each of two different hydroxyl groups on the sugar. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0034033 </td>
   <td style="text-align:left;vertical-align:top;"> purine nucleoside bisphosphate biosynthetic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways resulting in the formation of a purine nucleoside bisphosphate, a compound consisting of a purine base linked to a deoxyribose or ribose sugar esterified with one phosphate group attached to each of two different hydroxyl groups on the sugar. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0034112 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of homotypic cell-cell adhesion </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate, or extent of homotypic cell-cell adhesion. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0042168 </td>
   <td style="text-align:left;vertical-align:top;"> heme metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving heme, any compound of iron complexed in a porphyrin (tetrapyrrole) ring. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0043270 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of ion transport </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of the directed movement of charged atoms or small charged molecules into, out of or within a cell, or between cells, by means of some agent such as a transporter or pore. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045596 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of cell differentiation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of cell differentiation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045667 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of osteoblast differentiation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of osteoblast differentiation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045933 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of muscle contraction </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of muscle contraction. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0045934 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of nucleobase-containing compound metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> Any cellular process that stops, prevents, or reduces the frequency, rate or extent of the chemical reactions and pathways involving nucleobases, nucleosides, nucleotides and nucleic acids. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0046148 </td>
   <td style="text-align:left;vertical-align:top;"> pigment biosynthetic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways resulting in the formation of a pigment, any general or particular coloring matter in living organisms, e.g. melanin. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048870 </td>
   <td style="text-align:left;vertical-align:top;"> cell motility </td>
   <td style="text-align:left;vertical-align:top;"> Any process involved in the controlled self-propelled movement of a cell that results in translocation of the cell from one place to another. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0050678 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of epithelial cell proliferation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of epithelial cell proliferation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051234 </td>
   <td style="text-align:left;vertical-align:top;"> establishment of localization </td>
   <td style="text-align:left;vertical-align:top;"> Any process that localizes a substance or cellular component. This may occur via movement, tethering or selective degradation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051279 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of release of sequestered calcium ion into cytosol </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of the release into the cytosolic compartment of calcium ions sequestered in the endoplasmic reticulum or mitochondria. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0052548 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of endopeptidase activity </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of endopeptidase activity, the endohydrolysis of peptide bonds within proteins. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0055086 </td>
   <td style="text-align:left;vertical-align:top;"> nucleobase-containing small molecule metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The cellular chemical reactions and pathways involving a nucleobase-containing small molecule: a nucleobase, a nucleoside, or a nucleotide. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0071479 </td>
   <td style="text-align:left;vertical-align:top;"> cellular response to ionizing radiation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a ionizing radiation stimulus. Ionizing radiation is radiation with sufficient energy to remove electrons from atoms and may arise from spontaneous decay of unstable isotopes, resulting in alpha and beta particles and gamma rays. Ionizing radiation also includes X-rays. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0071560 </td>
   <td style="text-align:left;vertical-align:top;"> cellular response to transforming growth factor beta stimulus </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a transforming growth factor beta stimulus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0071702 </td>
   <td style="text-align:left;vertical-align:top;"> organic substance transport </td>
   <td style="text-align:left;vertical-align:top;"> The directed movement of organic substances into, out of or within a cell, or between cells, or within a multicellular organism by means of some agent such as a transporter or pore. An organic substance is a molecular entity that contains carbon. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0090130 </td>
   <td style="text-align:left;vertical-align:top;"> tissue migration </td>
   <td style="text-align:left;vertical-align:top;"> The process in which the population of cells that make up a tissue undergo directed movement. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0090330 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of platelet aggregation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the rate, frequency or extent of platelet aggregation. Platelet aggregation is the adhesion of one platelet to one or more other platelets via adhesion molecules. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:0098660 </td>
   <td style="text-align:left;vertical-align:top;"> inorganic ion transmembrane transport </td>
   <td style="text-align:left;vertical-align:top;"> The process in which an inorganic ion is transported across a membrane. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:1901385 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of voltage-gated calcium channel activity </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:1901731 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of platelet aggregation </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:1904752 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of vascular associated smooth muscle cell migration </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline </td>
   <td style="text-align:center;vertical-align:top;"> GO:1990542 </td>
   <td style="text-align:left;vertical-align:top;"> mitochondrial transmembrane transport </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0001667 </td>
   <td style="text-align:left;vertical-align:top;"> ameboidal-type cell migration </td>
   <td style="text-align:left;vertical-align:top;"> Cell migration that is accomplished by extension and retraction of a pseudopodium. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006163 </td>
   <td style="text-align:left;vertical-align:top;"> purine nucleotide metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving a purine nucleotide, a compound consisting of nucleoside (a purine base linked to a deoxyribose or ribose sugar) esterified with a phosphate group at either the 3' or 5'-hydroxyl group of the sugar. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006810 </td>
   <td style="text-align:left;vertical-align:top;"> transport </td>
   <td style="text-align:left;vertical-align:top;"> The directed movement of substances (such as macromolecules, small molecules, ions) or cellular components (such as complexes and organelles) into, out of or within a cell, or between cells, or within a multicellular organism by means of some agent such as a transporter or a transporter complex, a pore or a motor protein. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006814 </td>
   <td style="text-align:left;vertical-align:top;"> sodium ion transport </td>
   <td style="text-align:left;vertical-align:top;"> The directed movement of sodium ions (Na+) into, out of or within a cell, or between cells, by means of some agent such as a transporter or pore. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006942 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of striated muscle contraction </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of striated muscle contraction. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0006979 </td>
   <td style="text-align:left;vertical-align:top;"> response to oxidative stress </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of oxidative stress, a state often resulting from exposure to high levels of reactive oxygen species, e.g. superoxide anions, hydrogen peroxide (H2O2), and hydroxyl radicals. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007009 </td>
   <td style="text-align:left;vertical-align:top;"> plasma membrane organization </td>
   <td style="text-align:left;vertical-align:top;"> A process that is carried out at the cellular level which results in the assembly, arrangement of constituent parts, or disassembly of the plasma membrane. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007044 </td>
   <td style="text-align:left;vertical-align:top;"> cell-substrate junction assembly </td>
   <td style="text-align:left;vertical-align:top;"> The aggregation, arrangement and bonding together of a set of components to form a junction between a cell and its substrate. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007167 </td>
   <td style="text-align:left;vertical-align:top;"> enzyme-linked receptor protein signaling pathway </td>
   <td style="text-align:left;vertical-align:top;"> The series of molecular signals initiated by an extracellular ligand binding to a receptor on the surface of the target cell, where the receptor possesses catalytic activity or is closely associated with an enzyme such as a protein kinase, and ending with the regulation of a downstream cellular process, e.g. transcription. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007229 </td>
   <td style="text-align:left;vertical-align:top;"> integrin-mediated signaling pathway </td>
   <td style="text-align:left;vertical-align:top;"> The series of molecular signals initiated by an extracellular ligand binding to an integrin on the surface of a target cell, and ending with the regulation of a downstream cellular process, e.g. transcription. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007263 </td>
   <td style="text-align:left;vertical-align:top;"> nitric oxide mediated signal transduction </td>
   <td style="text-align:left;vertical-align:top;"> Any intracellular signal transduction in which the signal is passed on within the cell via nitric oxide (NO). Includes synthesis of nitric oxide, receptors/sensors for nitric oxide (such as soluble guanylyl cyclase/sGC) and downstream effectors that further transmit the signal within the cell. Nitric oxide transmits its downstream effects through either cyclic GMP (cGMP)-dependent or independent mechanisms. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007275 </td>
   <td style="text-align:left;vertical-align:top;"> multicellular organism development </td>
   <td style="text-align:left;vertical-align:top;"> The biological process whose specific outcome is the progression of a multicellular organism over time from an initial condition (e.g. a zygote or a young adult) to a later condition (e.g. a multicellular animal or an aged adult). </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007417 </td>
   <td style="text-align:left;vertical-align:top;"> central nervous system development </td>
   <td style="text-align:left;vertical-align:top;"> The process whose specific outcome is the progression of the central nervous system over time, from its formation to the mature structure. The central nervous system is the core nervous system that serves an integrating and coordinating function. In vertebrates it consists of the brain and spinal cord. In those invertebrates with a central nervous system it typically consists of a brain, cerebral ganglia and a nerve cord. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0007596 </td>
   <td style="text-align:left;vertical-align:top;"> blood coagulation </td>
   <td style="text-align:left;vertical-align:top;"> The sequential process in which the multiple coagulation factors of the blood interact, ultimately resulting in the formation of an insoluble fibrin clot; it may be divided into three stages: stage 1, the formation of intrinsic and extrinsic prothrombin converting principle; stage 2, the formation of thrombin; stage 3, the formation of stable fibrin polymers. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0008283 </td>
   <td style="text-align:left;vertical-align:top;"> cell population proliferation </td>
   <td style="text-align:left;vertical-align:top;"> The multiplication or reproduction of cells, resulting in the expansion of a cell population. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0008360 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of cell shape </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the surface configuration of a cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0009152 </td>
   <td style="text-align:left;vertical-align:top;"> purine ribonucleotide biosynthetic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways resulting in the formation of a purine ribonucleotide, a compound consisting of ribonucleoside (a purine base linked to a ribose sugar) esterified with a phosphate group at either the 3' or 5'-hydroxyl group of the sugar. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0009259 </td>
   <td style="text-align:left;vertical-align:top;"> ribonucleotide metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways involving a ribonucleotide, a compound consisting of ribonucleoside (a base linked to a ribose sugar) esterified with a phosphate group at either the 3' or 5'-hydroxyl group of the sugar. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0009611 </td>
   <td style="text-align:left;vertical-align:top;"> response to wounding </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a stimulus indicating damage to the organism. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0009653 </td>
   <td style="text-align:left;vertical-align:top;"> anatomical structure morphogenesis </td>
   <td style="text-align:left;vertical-align:top;"> The process in which anatomical structures are generated and organized. Morphogenesis pertains to the creation of form. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010595 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of endothelial cell migration </td>
   <td style="text-align:left;vertical-align:top;"> Any process that increases the rate, frequency, or extent of the orderly movement of an endothelial cell into the extracellular matrix to form an endothelium. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010632 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of epithelial cell migration </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of epithelial cell migration. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010718 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of epithelial to mesenchymal transition </td>
   <td style="text-align:left;vertical-align:top;"> Any process that increases the rate, frequency, or extent of epithelial to mesenchymal transition. Epithelial to mesenchymal transition is where an epithelial cell loses apical/basolateral polarity, severs intercellular adhesive junctions, degrades basement membrane components and becomes a migratory mesenchymal cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0010810 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of cell-substrate adhesion </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of cell-substrate adhesion. Cell-substrate adhesion is the attachment of a cell to the underlying substrate via adhesion molecules. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0014823 </td>
   <td style="text-align:left;vertical-align:top;"> response to activity </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of an activity stimulus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0022408 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of cell-cell adhesion </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents or reduces the rate or extent of cell adhesion to another cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0030032 </td>
   <td style="text-align:left;vertical-align:top;"> lamellipodium assembly </td>
   <td style="text-align:left;vertical-align:top;"> Formation of a lamellipodium, a thin sheetlike extension of the surface of a migrating cell. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0030193 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of blood coagulation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of blood coagulation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0030834 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of actin filament depolymerization </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of the disassembly of actin filaments by the removal of actin monomers from a filament. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0032272 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of protein polymerization </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of the process of creating protein polymers. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0032502 </td>
   <td style="text-align:left;vertical-align:top;"> developmental process </td>
   <td style="text-align:left;vertical-align:top;"> A biological process whose specific outcome is the progression of an integrated living unit: an anatomical structure (which may be a subcellular structure, cell, tissue, or organ), or organism over time from an initial condition to a later condition. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0034220 </td>
   <td style="text-align:left;vertical-align:top;"> ion transmembrane transport </td>
   <td style="text-align:left;vertical-align:top;"> A process in which an ion is transported across a membrane. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0034764 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of transmembrane transport </td>
   <td style="text-align:left;vertical-align:top;"> Any process that activates or increases the frequency, rate or extent of the directed movement of a solute from one side of a membrane to the other. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0042326 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of phosphorylation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents or decreases the rate of addition of phosphate groups to a molecule. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0042542 </td>
   <td style="text-align:left;vertical-align:top;"> response to hydrogen peroxide </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a hydrogen peroxide (H2O2) stimulus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0043149 </td>
   <td style="text-align:left;vertical-align:top;"> stress fiber assembly </td>
   <td style="text-align:left;vertical-align:top;"> The aggregation, arrangement and bonding together of a set of components to form a stress fiber. A stress fiber is a contractile actin filament bundle that consists of short actin filaments with alternating polarity. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0043154 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of cysteine-type endopeptidase activity involved in apoptotic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of a cysteine-type endopeptidase activity involved in the apoptotic process. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0043542 </td>
   <td style="text-align:left;vertical-align:top;"> endothelial cell migration </td>
   <td style="text-align:left;vertical-align:top;"> The orderly movement of an endothelial cell into the extracellular matrix to form an endothelium. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0046390 </td>
   <td style="text-align:left;vertical-align:top;"> ribose phosphate biosynthetic process </td>
   <td style="text-align:left;vertical-align:top;"> The chemical reactions and pathways resulting in the formation of ribose phosphate, any phosphorylated ribose sugar. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048259 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of receptor-mediated endocytosis </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of receptor mediated endocytosis, the uptake of external materials by cells, utilizing receptors to ensure specificity of transport. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048519 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of biological process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of a biological process. Biological processes are regulated by many means; examples include the control of gene expression, protein modification or interaction with a protein or substrate molecule. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048856 </td>
   <td style="text-align:left;vertical-align:top;"> anatomical structure development </td>
   <td style="text-align:left;vertical-align:top;"> The biological process whose specific outcome is the progression of an anatomical structure from an initial condition to its mature state. This process begins with the formation of the structure and ends with the mature structure, whatever form that may be including its natural destruction. An anatomical structure is any biological entity that occupies space and is distinguished from its surroundings. Anatomical structures can be macroscopic such as a carpel, or microscopic such as an acrosome. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0048869 </td>
   <td style="text-align:left;vertical-align:top;"> cellular developmental process </td>
   <td style="text-align:left;vertical-align:top;"> A biological process whose specific outcome is the progression of a cell over time from an initial condition to a later condition. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0050819 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of coagulation </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of coagulation. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0050878 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of body fluid levels </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the levels of body fluids. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051017 </td>
   <td style="text-align:left;vertical-align:top;"> actin filament bundle assembly </td>
   <td style="text-align:left;vertical-align:top;"> The assembly of actin filament bundles; actin filaments are on the same axis but may be oriented with the same or opposite polarities and may be packed with different levels of tightness. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051248 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of protein metabolic process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of chemical reactions and pathways involving a protein. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051283 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of sequestering of calcium ion </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of the binding or confining calcium ions such that they are separated from other components of a biological system. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051301 </td>
   <td style="text-align:left;vertical-align:top;"> cell division </td>
   <td style="text-align:left;vertical-align:top;"> The process resulting in division and partitioning of components of a cell to form more cells; may or may not be accompanied by the physical separation of a cell into distinct, individually membrane-bounded daughter cells. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051494 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of cytoskeleton organization </td>
   <td style="text-align:left;vertical-align:top;"> Any process that stops, prevents, or reduces the frequency, rate or extent of the formation, arrangement of constituent parts, or disassembly of cytoskeletal structures. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051651 </td>
   <td style="text-align:left;vertical-align:top;"> maintenance of location in cell </td>
   <td style="text-align:left;vertical-align:top;"> Any process in which a substance or cellular entity, such as a protein complex or organelle, is maintained in a specific location within, or in the membrane of, a cell, and is prevented from moving elsewhere. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051693 </td>
   <td style="text-align:left;vertical-align:top;"> actin filament capping </td>
   <td style="text-align:left;vertical-align:top;"> The binding of a protein or protein complex to the end of an actin filament, thus preventing the addition, exchange or removal of further actin subunits. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0051893 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of focal adhesion assembly </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of focal adhesion formation, the establishment and maturation of focal adhesions. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0060048 </td>
   <td style="text-align:left;vertical-align:top;"> cardiac muscle contraction </td>
   <td style="text-align:left;vertical-align:top;"> Muscle contraction of cardiac muscle tissue. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0061045 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of wound healing </td>
   <td style="text-align:left;vertical-align:top;"> Any process that decreases the rate, frequency, or extent of the series of events that restore integrity to a damaged tissue, following an injury. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0070527 </td>
   <td style="text-align:left;vertical-align:top;"> platelet aggregation </td>
   <td style="text-align:left;vertical-align:top;"> The adhesion of one platelet to one or more other platelets via adhesion molecules. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0071495 </td>
   <td style="text-align:left;vertical-align:top;"> cellular response to endogenous stimulus </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a stimulus arising within the organism. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0090092 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of transmembrane receptor protein serine/threonine kinase signaling pathway </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the rate, frequency, or extent of the series of molecular signals generated as a consequence of a transmembrane receptor serine/threonine kinase binding to its physiological ligand. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0090132 </td>
   <td style="text-align:left;vertical-align:top;"> epithelium migration </td>
   <td style="text-align:left;vertical-align:top;"> The process in which the population of cells that make up an epithelium undergo directed movement. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0090257 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of muscle system process </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the frequency, rate or extent of a muscle system process, a multicellular organismal process carried out by any of the organs or tissues in a muscle system. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0090287 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of cellular response to growth factor stimulus </td>
   <td style="text-align:left;vertical-align:top;"> Any process that modulates the rate, frequency, or extent of a change in state or activity of a cell (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a growth factor stimulus. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0090303 </td>
   <td style="text-align:left;vertical-align:top;"> positive regulation of wound healing </td>
   <td style="text-align:left;vertical-align:top;"> Any process that increases the rate, frequency, or extent of the series of events that restore integrity to a damaged tissue, following an injury. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0097553 </td>
   <td style="text-align:left;vertical-align:top;"> calcium ion transmembrane import into cytosol </td>
   <td style="text-align:left;vertical-align:top;"> A process in which a calcium ion is transported from one side of a membrane to the other into the cytosol by means of some agent such as a transporter or pore. </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:0150116 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of cell-substrate junction organization </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1900024 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of substrate adhesion-dependent cell spreading </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1900047 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of hemostasis </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1902904 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of supramolecular fiber organization </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1903034 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of response to wounding </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:1903522 </td>
   <td style="text-align:left;vertical-align:top;"> regulation of blood circulation </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Baseline, 3 d. </td>
   <td style="text-align:center;vertical-align:top;"> GO:2001237 </td>
   <td style="text-align:left;vertical-align:top;"> negative regulation of extrinsic apoptotic signaling pathway </td>
   <td style="text-align:left;vertical-align:top;">  </td>
  </tr>
  <tr>
   <td style="text-align:left;vertical-align:top;"> Change </td>
   <td style="text-align:center;vertical-align:top;"> GO:0070301 </td>
   <td style="text-align:left;vertical-align:top;"> cellular response to hydrogen peroxide </td>
   <td style="text-align:left;vertical-align:top;"> Any process that results in a change in state or activity of a cell (in terms of movement, secretion, enzyme production, gene expression, etc.) as a result of a hydrogen peroxide (H2O2) stimulus. </td>
  </tr>
</tbody>
</table>

<br>
  
  ## **References**
  
  <br>
