#!/bin/bash

r_file=$(basename $0 .sh).R
R=/home/hilmar/shared/R/4.3_no_X/bin/R

(cat <<"EOF"
.libPaths("/home/hilmar/R/4.3")
library(knitr)
rmarkdown::render("Data_import.Rmd")
rmarkdown::render("Doublet_Detection.Rmd")
rmarkdown::render("Preprocessing_Harmony.Rmd")
rmarkdown::render("Cluster_identification_all.Rmd")
rmarkdown::render("Preprocessing_Harmony_Epithelial_v4_CCadj.Rmd")
rmarkdown::render("Cluster_identification_Epithelial_v4.Rmd")
rmarkdown::render("Cluster_identification_Epithelial_v4.2.Rmd")
rmarkdown::render("Preprocessing_Harmony_Epithelial_v5_antral_CCadj.Rmd")
rmarkdown::render("Cluster_identification_Epithelial_v5.1.Rmd")
rmarkdown::render("DGE_Epithelial_v4.2_CCadj.Rmd")
rmarkdown::render("DGE_Epithelial_v5.1_CCadj.Rmd")
rmarkdown::render("Final_visualizations_v4.2.Rmd")
rmarkdown::render("Final_visualizations_antral_only_v5.1.Rmd")
quit("no")
EOF
) > $r_file

$R CMD BATCH $r_file

