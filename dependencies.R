
# select mirrors
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")


# install packages from CRAN
Packages <- c('DT', 'shiny', 'reshape2', 'stringi', 'ggplot2', 'markdown',
              'magrittr', 'shinydashboard', 'shinycssloaders', 'ggseqlogo')
install.packages(Packages)

# install packages from bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ShortRead")
