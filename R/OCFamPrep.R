
##########################################################################################
#OCFam_fam_K_matrix
##########################################################################################



rm(list=ls())
setwd("C:/Users/MHamilton/CGIAR/WorldFish DocuShare - Genetics/Data/GIFT/Spawn_rounds/G202402_TiLV_G20/08 OC")
ped <- read.csv("ped2.csv")

ped <- ped[,c("INDIV", "SIRE", "DAM", "FAM", "BORN", "EBV")]

ped[ped[,"SIRE"] == "","SIRE"] <- 0
ped[ped[,"DAM"] == "","DAM"] <- 0









##########################################################################################
#OCFamPrep
##########################################################################################


#' @param ped is a data frame with the following columns (class in parentheses):
#' \itemize{
#'  \item{'INDIV' is the individual identifier (character).}
#'  \item{'SIRE' is the male parent identifier (character).}
#'  \item{'DAM' is the female parent identifier (character).}
#'  \item{'FAM' is a full-sibling family identifier (character).}
#'  \item{'SEX' is the sex of the individual. "male" (or 1 = male), "female" (or 2 = female), NA if unknown.  If one animal is of unknown sex, sex will be ignored (integer or character).}
#'  \item{'BORN' integer indicating age class.  May be the year of birth if one age class per year or an integer indicating the sequence of age classes (integer).}
#'  \item{'EBV' is the estimated breeding value (numeric).}
#'  \item{'N_AS_PARENT_PREV' is the number of families in the next age class previously contributed to (integer).}
#'  \item{'AVAIL_BROOD' is TRUE if the individual is candidate parent of next age class (logical).}
#'  }
#'
#'
#' @param gene_flow_vector is a applicable to overlapping generations (if NA discrete generations is assumed).  It is vector representing parental contributions by age class to the next age class. For example, gene_flow_vector = c(0.2, 0.8, 0, 0) - oldest age class to youngest age class.

#Required packages: AGHmatrix, dplyr, ggplot2


OCFam_fam_K_matrix <- function(family_dat) {
  #data
  # class: data.frame
  # fields:
  #   FAM
  #   SIRE
  #   FAM_SIRE
  #   DAM
  #   FAM_DAM

  #sires
  sires <- matrix(unique(family_dat$SIRE), ncol = 1)
  colnames(sires)[1] <- "SIRE"
  sires <- merge(sires, family_dat[,c("SIRE","FAM_SIRE")], by = "SIRE", all.x = TRUE)
  colnames(sires) <- c("Indiv_id","FAM")
  sires <- merge(sires, family_dat[,c("FAM","SIRE","DAM")], by = "FAM", all.x = TRUE)
  sires <- subset(sires, select=-c(FAM)) #remove FAM column

  #dams
  dams <- matrix(unique(family_dat$DAM), ncol = 1)
  colnames(dams)[1] <- "DAM"
  dams <- merge(dams, family_dat[,c("DAM","FAM_DAM")], by = "DAM", all.x = TRUE)
  colnames(dams) <- c("Indiv_id","FAM")
  dams <- merge(dams, family_dat[,c("FAM","SIRE","DAM")], by = "FAM", all.x = TRUE)
  dams <- subset(dams, select=-c(FAM)) #remove FAM column

  parents <- rbind(sires,dams)
  parents <- unique(parents)
  parents <- parents[order(parents[,"Indiv_id"], decreasing = FALSE), ]

  families            <- unique(family_dat)
  tmp1 <- families
  tmp2 <- families
  tmp1$FAM <- paste0("1_", tmp1$FAM) #needs a unique id
  families <- rbind(tmp1, tmp2)

  fam_pedigree        <- rbind(as.matrix(parents),as.matrix(families[,c("FAM","SIRE","DAM")]))

  #Generate K Matrix
  # K_matrix_families      <- nadiv::makeA(fam_pedigree[,1:3])/2
  tmp <- fam_pedigree
  tmp[is.na(tmp)] <- 0
  K_matrix_families      <- AGHmatrix::Amatrix(tmp[,1:3], ploidy=2)/2
  rm(tmp)

  family_inbreeding <-  round((diag(as.matrix(K_matrix_families)) * 2) - 1, 10)
  family_inbreeding <-  family_inbreeding[(nrow(parents)+1):(nrow(parents)+nrow(families)/2)]

  K_matrix_families <- K_matrix_families[(nrow(parents)+1):(nrow(parents)+nrow(families)/2), (nrow(parents)+nrow(families)/2+1):nrow(K_matrix_families)]

  rownames(K_matrix_families) <- colnames(K_matrix_families)
  names(family_inbreeding) <- colnames(K_matrix_families)

  return(list(K_matrix_families = K_matrix_families,
              family_inbreeding = family_inbreeding))
}



OCFamPrep <- function(ped, gene_flow_vector = NA) {
  library(ggplot2)
  #library(dplyr)
  #library(AGHmatrix)

  #  family_dat <- family_dat[,c("FAM", "SIRE", "DAM", "FAM_SIRE", "FAM_DAM", "SPAWN_ROUND", "SPECIES", "FAM_NEST")]

  #create family_dat
  family_dat <- unique(ped[,c("FAM", "SIRE", "DAM", 'BORN')])
  sires <- ped[,c('INDIV', 'FAM', "EBV")]
  colnames(sires) <- c("SIRE", "FAM_SIRE", "EBV_SIRE")
  family_dat <- dplyr::left_join(family_dat, sires, by = "SIRE")
  rm(sires)
  dams <- ped[,c('INDIV', 'FAM', "EBV")]
  colnames(dams) <- c("DAM", "FAM_DAM", "EBV_DAM")
  family_dat <- dplyr::left_join(family_dat, dams, by = "DAM")
  rm(dams)

  family_dat$EBV_FAM <- (family_dat$EBV_SIRE + family_dat$EBV_DAM)/2

  family_dat[family_dat$FAM == "", "FAM"] <- NA
  family_dat <- family_dat[!is.na(family_dat$FAM),]

  ped <- dplyr::left_join(ped, family_dat, by = c("FAM", "SIRE", "DAM", 'BORN'))

  fam_K_matrix <- OCFam_fam_K_matrix(family_dat)

  family_inbreeding <- data.frame(FAM = names(fam_K_matrix$family_inbreeding),
                                  F    = fam_K_matrix$family_inbreeding)
  fam_K_matrix <- fam_K_matrix$K_matrix_families
  fam_K_matrix <- as.data.frame(as.matrix(fam_K_matrix))

  fam_K_matrix$FAM <- rownames(fam_K_matrix)

  ####################################################################################
  #Generate plots
  ####################################################################################

  age_classs <- unique(family_dat$BORN)

  age_class_K <- data.frame(BORN = NA, VARIABLE = NA, MEAN = NA)[-1,]
  age_class_F <- age_class_K
  #get means
  for(sr in age_classs) {
    fams <- family_dat[family_dat$BORN == sr,"FAM"]

    tmp_mean <- mean(as.matrix(fam_K_matrix[fams,fams]))
    tmp <- data.frame(BORN = sr, VARIABLE = "Coancestry (k)", MEAN = tmp_mean)
    age_class_K <- rbind(age_class_K, tmp)

    tmp_F <- mean(family_inbreeding[family_inbreeding$FAM %in% fams, "F"])
    tmp <- data.frame(BORN = sr, VARIABLE = "Inbreeding (F)", MEAN = tmp_F)
    age_class_F <- rbind(age_class_F, tmp)

  }

  age_class_means <- rbind(age_class_K, age_class_F)

  if(is.na(gene_flow_vector)) {
    p <- ggplot(age_class_means, aes(x=BORN, y=MEAN, group=VARIABLE)) +
      geom_line(aes(color=VARIABLE))+
      geom_point(aes(color=VARIABLE)) +
      #  theme_bw()  +
      theme_classic() +
      theme( # remove the vertical grid lines
        panel.grid.major.x = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.y = element_line(linewidth=.1, color="grey70" )
      )
    p
  }
  age_class_means <- data.frame(BORN = unique(age_class_means$BORN),
                                COANCESTRY = age_class_means[age_class_means$VARIABLE == "Coancestry (k)","MEAN"],
                                INBREEDING = age_class_means[age_class_means$VARIABLE == "Inbreeding (F)","MEAN"])

  return(list(plot = p,
              fam_K_matrix = fam_K_matrix,
              family_inbreeding = family_inbreeding,
              age_class_means = age_class_means))
}




##Export plot
#file_name <- paste0(tempdir(),"/",
#                    gsub(":","",Sys.time()), "_",
#                    "_mean_coancestry_and_inbreeding.png")
#ggsave(filename = file_name, plot = out$plot, device = "png")
#shell.exec(file_name)

#Export as Excel file

#https://trinkerrstuff.wordpress.com/2018/02/14/easily-make-multi-tabbed-xlsx-files-with-openxlsx/

#fam_K_matrix <- bind_cols(data.frame(Family_K_matrix = rownames(out$fam_K_matrix)),as.data.frame(out$fam_K_matrix))

#dat <- list(out$fam_K_matrix,
#            out$family_inbreeding,
#            out$age_class_means)

#names(dat) <- c(paste0("fam_K_matrix"),
#                paste0("family_inbreeding"),
#                paste0("age_class_means"))

## Create a blank workbook
#wb <- createWorkbook()

## Loop through the list of split tables as well as their names
##   and add each one as a sheet to the workbook
#Map(function(data, name){
#  addWorksheet(wb, name)
#  writeData(wb, name, data)
#}, dat, names(dat))

#file_name <- paste0(tempdir(),"/",
#                    gsub(":","",Sys.time()),
#                    " K matrices.xlsx")
#saveWorkbook(wb, file_name, overwrite = TRUE) #save Excel workbook
#shell.exec(file_name)






