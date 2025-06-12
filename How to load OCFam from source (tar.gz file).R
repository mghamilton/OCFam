#Loading OCFam from source (tar.gz file)

if("package:OCFam" %in% search()) {detach("package:OCFam", unload = TRUE)}
remove.packages("OCFam")
install.packages("path/to/package_name.tar.gz", repos = NULL, type = "source")
library(OCFam)