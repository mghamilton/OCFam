# 1. Introduction

*OCFam* is an R package implementing optimal contribution selection
(OCS) for highly fecund species (e.g., fish), where large family sizes
and rounding constraints renders classical individual-level OCS
impractical. The package provides two primary functions:

1.  *OCFamPrep()*: Prepares pedigree data and summarizes historical
    trends in inbreeding and coancestry.

2.  *OCFam()*: Performs optimal contribution selection and generates
    integer-feasible allocations of parents.

*OCFam* addresses the gap between the continuous mathematical optimum
and operational constraints on parental contributions, ensuring that
selected contributions are close to optimal and practically
implementable. The package supports both non-overlapping (discrete) and
overlapping generation breeding systems.

# 

# 2. Installation of OCFam

## R installation

1.  Install R from CRAN (<https://cran.r-project.org/>)

2.  Install the corresponding version of Rtools from CRAN Rtools
    (required for building packages from source on Windows)
    (<https://cran.r-project.org/bin/windows/Rtools/>).

3.  (Optional) Install RStudio Desktop from Posit for an enhanced
    development environment
    (<https://posit.co/download/rstudio-desktop/>).

## OCFam installation

*OCFam* is installed from GitHub. Start or restart R, after installing
Rtools, and then enter the following R code to install *OCFam* and its
dependencies:

> if (!requireNamespace(\"optiSel\", quietly = TRUE))
> install.packages(\"optiSel\")
>
> if (!requireNamespace(\"AGHmatrix\", quietly = TRUE))
> install.packages(\"AGHmatrix\")
>
> if (!requireNamespace(\"dplyr\", quietly = TRUE))
> install.packages(\"dplyr\")
>
> if (!requireNamespace(\"ggplot2\", quietly = TRUE))
> install.packages(\"ggplot2\")
>
> if (!requireNamespace(\"devtools\", quietly = TRUE))
> install.packages(\"devtools\")
>
> \# Install OCFam from GitHub
>
> devtools::install_github(\"mghamilton/OCFam\", dependencies = TRUE,
> force = TRUE)
>
> library(OCFam)
