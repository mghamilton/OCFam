
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
#'#' @param age_class_names is an data frame specifying names for age classes (class in parentheses):
#' \itemize{
#'  \item{'AGE_CLASS_NAME' is a name for the age class (character).}
#'  \item{'BORN' integer indicating age class.  May be the year of birth if one age class per year or an integer indicating the sequence of age classes (integer).}
#'  }
#'
#' @param gene_flow_vector is a applicable to overlapping generations (if NA discrete generations is assumed).  It is vector representing parental contributions by age class to the next age class. For example, gene_flow_vector = c(0.2, 0.8, 0, 0) - oldest age class to youngest age class.

#Required packages: AGHmatrix, dplyr, ggplot2


OCFam_fam_K_matrix <- function(family_ped) {
  #data
  # class: data.frame
  # fields:
  #   FAM
  #   SIRE
  #   FAM_SIRE
  #   DAM
  #   FAM_DAM

  #sires
  sires <- matrix(unique(family_ped$SIRE), ncol = 1)
  colnames(sires)[1] <- "SIRE"
  sires <- merge(sires, family_ped[,c("SIRE","FAM_SIRE")], by = "SIRE", all.x = TRUE)
  colnames(sires) <- c("Indiv_id","FAM")
  sires <- merge(sires, family_ped[,c("FAM","SIRE","DAM")], by = "FAM", all.x = TRUE)
  sires <- subset(sires, select=-c(FAM)) #remove FAM column

  #dams
  dams <- matrix(unique(family_ped$DAM), ncol = 1)
  colnames(dams)[1] <- "DAM"
  dams <- merge(dams, family_ped[,c("DAM","FAM_DAM")], by = "DAM", all.x = TRUE)
  colnames(dams) <- c("Indiv_id","FAM")
  dams <- merge(dams, family_ped[,c("FAM","SIRE","DAM")], by = "FAM", all.x = TRUE)
  dams <- subset(dams, select=-c(FAM)) #remove FAM column

  parents <- rbind(sires,dams)
  parents <- unique(parents)
  parents <- parents[order(parents[,"Indiv_id"], decreasing = FALSE), ]

  families            <- unique(family_ped)
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



OCFamPrep <- function(ped, age_class_names = NULL, gene_flow_vector = NULL) {
  library(ggplot2)
  if(length(gene_flow_vector) == 1) {gene_flow_vector <- NULL}
  if(length(age_class_names) == 1) {age_class_names <- NULL}
  #library(dplyr)
  #library(AGHmatrix)

  #  family_ped <- family_ped[,c("FAM", "SIRE", "DAM", "FAM_SIRE", "FAM_DAM", "SPAWN_ROUND", "SPECIES", "FAM_NEST")]

  #create family_ped
  family_ped <- unique(ped[,c("FAM", "SIRE", "DAM", 'BORN')])
  sires <- ped[,c('INDIV', 'FAM', "EBV")]
  colnames(sires) <- c("SIRE", "FAM_SIRE", "EBV_SIRE")
  family_ped <- dplyr::left_join(family_ped, sires, by = "SIRE")
  rm(sires)
  dams <- ped[,c('INDIV', 'FAM', "EBV")]
  colnames(dams) <- c("DAM", "FAM_DAM", "EBV_DAM")
  family_ped <- dplyr::left_join(family_ped, dams, by = "DAM")
  rm(dams)

  family_ped$EBV_FAM <- (family_ped$EBV_SIRE + family_ped$EBV_DAM)/2

  family_ped[family_ped$FAM == "" | is.na(family_ped$FAM), "FAM"] <- NA
  family_ped <- family_ped[!is.na(family_ped$FAM),]

  ped <- dplyr::left_join(ped, family_ped, by = c("FAM", "SIRE", "DAM", 'BORN'))

  fam_K_matrix <- OCFam_fam_K_matrix(family_ped)

  family_inbreeding <- data.frame(FAM = names(fam_K_matrix$family_inbreeding),
                                  F    = fam_K_matrix$family_inbreeding)
  fam_K_matrix <- fam_K_matrix$K_matrix_families
  fam_K_matrix <- as.data.frame(as.matrix(fam_K_matrix))

  fam_K_matrix$FAM <- rownames(fam_K_matrix)

  ####################################################################################
  #Generate plots
  ####################################################################################

  age_classs <- unique(family_ped$BORN)

  age_class_K_mat <- matrix(data = NA, nrow = length(age_classs), ncol = length(age_classs))
  rownames(age_class_K_mat) <- age_classs
  colnames(age_class_K_mat) <- age_classs

  #get means

  age_class_F <- data.frame(BORN = NA, VARIABLE = NA, MEAN = NA)[-1,]
  for(ac in age_classs) {
    fams <- family_ped[family_ped$BORN == as.character(ac),"FAM"]
    tmp_F <- mean(family_inbreeding[family_inbreeding$FAM %in% fams, "F"])
    tmp <- data.frame(BORN = ac, VARIABLE = "Inbreeding (F)", MEAN = tmp_F)
    age_class_F <- rbind(age_class_F, tmp)
  }

  #coancestry matrix
  for(ac1 in age_classs) {
    for(ac2 in age_classs) {

      if(ac1 >= ac2){
        fams1 <- family_ped[family_ped$BORN == as.character(ac1),"FAM"]
        fams2 <- family_ped[family_ped$BORN == as.character(ac2),"FAM"]

        tmp_mean <- mean(as.matrix(fam_K_matrix[fams1,fams2]))
        age_class_K_mat[as.character(ac1), as.character(ac2)] <-  tmp_mean
      }
    }
  }
  age_class_K_mat[upper.tri(age_class_K_mat)] <- age_class_K_mat[lower.tri(age_class_K_mat)]
  age_class_K <- data.frame(BORN = rownames(age_class_K_mat), VARIABLE = "Coancestry (k)", MEAN = diag(age_class_K_mat))

  if(is.null(gene_flow_vector)) {
    age_class_means <- rbind(age_class_F, age_class_K)

  }

    age_class_K_overlap <-  data.frame(BORN = age_classs, VARIABLE = "Coancestry (overlapping; k)", MEAN = NA)

  if(!is.null(gene_flow_vector)) {
    for(ac in age_classs[length(gene_flow_vector):length(age_classs)]) {
      age_class_K_overlap[age_class_K_overlap$BORN == ac,"MEAN"] <-  t(gene_flow_vector) %*%
        age_class_K_mat[as.character((ac-length(gene_flow_vector)+1):ac),
                        as.character((ac-length(gene_flow_vector)+1):ac)] %*%
        gene_flow_vector
    }
    age_class_means <- rbind(age_class_F, age_class_K, age_class_K_overlap)
  }

    age_class_means$BORN <-  as.numeric(age_class_means$BORN)

   if(is.null(age_class_names)) {
  p <- ggplot(age_class_means, aes(x=BORN, y=MEAN, group=VARIABLE))
   } else {
     age_class_means <- dplyr::left_join(age_class_means, age_class_names, by = "BORN")
     ped <- dplyr::left_join(ped, age_class_names, by = "BORN")
     family_ped <- dplyr::left_join(family_ped, age_class_names, by = "BORN")

     # Ensure AGE_CLASS_NAME is ordered by BORN
     age_class_means <- age_class_means[order(age_class_means$BORN), ]

     # Reassign AGE_CLASS_NAME as a factor, now ordered by BORN
     age_class_means$AGE_CLASS_NAME <- factor(age_class_means$AGE_CLASS_NAME,
                                              levels = unique(age_class_means$AGE_CLASS_NAME))

     p <- ggplot(age_class_means, aes(x=AGE_CLASS_NAME, y=MEAN, group=VARIABLE))

   }
    p <- p +
      geom_line(aes(color=VARIABLE))+
      geom_point(aes(color=VARIABLE)) +
      #  theme_bw()  +
      theme_classic() +
      theme( # remove the vertical grid lines
        panel.grid.major.x = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.y = element_line(linewidth=.1, color="grey70" ),
        axis.text.x = element_text(angle = 90, vjust = 0.5)
      )
  p

  age_class_means <- data.frame(BORN = age_classs,
                                INBREEDING = age_class_F$MEAN,
                                COANCESTRY = age_class_K$MEAN,
                                OVERLAPPING_COANCESTRY = age_class_K_overlap$MEAN)

  return(list(plot = p,
              fam_K_matrix = fam_K_matrix,
              age_class_K_mat = age_class_K_mat,
              family_inbreeding = family_inbreeding,
              age_class_means = age_class_means,
              ped = ped,
              family_ped = family_ped,
              gene_flow_vector = gene_flow_vector))
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






