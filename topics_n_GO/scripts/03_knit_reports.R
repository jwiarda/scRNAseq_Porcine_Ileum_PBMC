library(rmarkdown)
library(knitr)

# Knit all Rmd files to html
rmarkdown::render(input = './RMD/Topic_Models_Ileum_BCells.Rmd')
rmarkdown::render(input = './RMD/Topic_Models_Ileum_CD4T.Rmd')
rmarkdown::render(input = './RMD/Topic_Models_Ileum_gdCD8T.Rmd')
rmarkdown::render(input = './RMD/Topic_Models_Ileum_ILC.Rmd')
rmarkdown::render(input = './RMD/Topic_Models_Ileum_Myeloid.Rmd')
rmarkdown::render(input = './RMD/Topic_Models_Ileum_NonImmune.Rmd')
rmarkdown::render(input = './RMD/Topic_Models_Ileum_TILC.Rmd')
rmarkdown::render(input = './RMD/Topic_Models_Ileum_Tnaive.Rmd')
rmarkdown::render(input = './RMD/Topic_Models_Ileum_AtlasALL.Rmd')
rmarkdown::render(input = './RMD/Topic_Models_ILC_bloodGut.Rmd')