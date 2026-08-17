#!/bin/bash

r_file=$(basename $0 .sh).R
R=/home/hilmar/shared/R/4.3_no_X/bin/R

(cat <<"EOF"
.libPaths("/home/hilmar/R/4.3")
library(knitr)
rmarkdown::render("Preprocess_raw_files.Rmd")
rmarkdown::render("ExpressionAnalysis.Rmd")
rmarkdown::render("Preprocessing_Harmony.Rmd")
rmarkdown::render("GSEA_Analysis.Rmd")
rmarkdown::render("GSEA_plots.Rmd")
rmarkdown::render("Final_Volcano.Rmd")
quit("no")
EOF
) > $r_file

$R CMD BATCH $r_file

