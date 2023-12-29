# ============================== CODE METADATA =============================== #
# DATE CREATED: 2022-03-31
# ADAPTATION: Catalina Manduta
# DESCRIPTION:
              # Installation of the needed packages
              # Cleaning, organize the data prior to BayesPGLS regression
              # BayesPGLS regression implementation
              # Plotting implementation for visualization 

# ==== INSTALL R PACKAGES: ====

# System information:
sysinf <- Sys.info()

setwd("/work/CatalinaThesis")

install.packages("RhpcBLASctl")
library(RhpcBLASctl)
blas_set_num_threads(1)

# Install libraries:
if (sysinf["sysname"] == "Linux") {
  instPacks <- installed.packages()[, 1]
  if (!"snowfall" %in% instPacks) {
    install.packages("snowfall")
  }
  if (!"caper" %in% instPacks) {
    install.packages("caper")
  }
  if (!"phytools" %in% instPacks) {
    install.packages("phytools")
  }
  if (!"RColorBrewer" %in% instPacks) {
    install.packages("RColorBrewer")
  }
  
  if (!"BayesPGLS" %in% instPacks) {
    install.packages("/work/CatalinaThesis/analysis/packages/BayesPGLS_1.0.1.tar.gz", 
                     type = "source", repos = NULL)
  }
  
  # Set one thread (core) for BLAS operations: 
  Sys.setenv(OMP_NUM_THREADS = 1)
  Sys.setenv(MKL_NUM_THREADS = 1)
  
  # Main directory:
  analysisDir <- "/work/CatalinaThesis/analysis/packages/BayesPGLS_1.0.1.tar.gz"
  setwd(analysisDir)
  
} else {
  # Main directory:
  analysisDir <- "/work/CatalinaThesis/analysis/packages/BayesPGLS_1.0.1.tar.gz"
  
}

# Install additional packeges
install.packages("flextable")
install.packages("ggplot2")
install.packages("xlsx")
install.packages("writexl")


# Libraries:
library(BayesPGLS)
library(snowfall)
library(caper)
library(phytools)
library(RColorBrewer)
library(flextable)
library(readxl)
library(readr)
library(ggplot2)
library(xlsx)
library(writexl)

# Logical for saving results:
saveResulst <- FALSE

# ==== NAMING CONSENSUS: ====
vertlife <-  rep(NA, nrow(lifeExpBaSTA))
lifeExpBaSTA_new <- cbind(lifeExpBaSTA[,1:1], vertlife, lifeExpBaSTA[,2:ncol(lifeExpBaSTA)])

# Delete all the species which are not mammals
for (i in nrow(mammals_vertlife_taxonomies):1){
  if (! (mammals_vertlife_taxonomies[i, "group"] == "Mammals")) {
    mammals_vertlife_taxonomies <- mammals_vertlife_taxonomies[-i,]
  }
}

# Check to see if the species name from the lifeExpBaSTA_new is present also in the 
# mammals_vertlife_taxonomies
for (i in 1:nrow(lifeExpBaSTA_new)){
  if (lifeExpBaSTA_new[i, "species"] %in% mammals_vertlife_taxonomies$scientificname) {
    lifeExpBaSTA_new$vertlife[i] <- lifeExpBaSTA_new[i, "species"]
  }
}

# Select only the values that are NA in vertlife section in order to check only those manually 
lifeExpBaSTA_new_na <- lifeExpBaSTA_new[is.na(lifeExpBaSTA_new$vertlife), ]
lifeExpBaSTA_new$vertlife[9] <- "Aonyx cinerea"
lifeExpBaSTA_new$vertlife[20] <- "Bubalus arnee"
lifeExpBaSTA_new$vertlife[37] <- "Cervus elaphus"
lifeExpBaSTA_new$vertlife[54] <- "Equus africanus"

lifeExpBaSTA_new_na <- lifeExpBaSTA_new[is.na(lifeExpBaSTA_new$vertlife), ] # check to see if there is anythings else that must be checked
vert <- as.data.frame(lifeExpBaSTA_new$vertlife)



write.xlsx(vert,file='Species.xlsx')
write.csv(lifeExpBaSTA_new,file='lifeExpBaSTA_new.csv')


# ==== DATA: ====
taxa <- "Mammalia"
# Load life history variables:

lifeExpBaSTA_new <- read.csv("/work/CatalinaThesis/analysis/data/tables/lifeExpBaSTA_new.csv")

predDat <- read.csv(file = "/work/CatalinaThesis/analysis/data/tables/LifeHistTraitsMammals_02_wide.csv", # replace file with your life history data file / match taxonomy
                    header = TRUE, stringsAsFactors = FALSE)

# Load phylogeny:
load(sprintf("/work/CatalinaThesis/analysis/data/rdata/trees/maxCredTree%s.RData", taxa)) 

phylo <- maxCred

# Calculate castration difference:
respDat <- as.data.frame((lifeExpBaSTA_new$Male_None_Mean - lifeExpBaSTA_new$Male_Surgical_Mean)/ lifeExpBaSTA_new$Male_None_Mean) # calculate the differences between castrated and non castrated
colnames(respDat)[1] <- "castration_difference"


# ====DATA PREPARATION: ====
# Number of species with predictor data:
nPred <- nrow(predDat)

# Sexes:
sexes <- c(f = "Female", m = "Male")
nsex <- length(sexes) #[1] 2

# Data bases:
dbs <- c("ZIMS", "amni", "jm", "anage", "panth", "datlife")
ndbs <- length(dbs) # [1] 6

# Variables:
varNames <- c("bm.f", "bm.m", "asm.f", "asm.m", "lcs", "lc.y", "ilibi")
varLables <- c("BM_Female", "BM_Male", "ASM_Female", "ASM_Male", 
               "LCS", "LCY", "IBI")
nvars <- length(varNames) # [1] 7

# Colnames of predictors:
allPredNames <- colnames(predDat)
# [1] "gbif.species"          "zims.species"          "vertlife.treename"     "class"                
# [5] "matsys.socialmonogamy" "testes.mass.combined"  "parental.care"         "bm.f.ZIMS"            
# [9] "bm.m.ZIMS"             "sdbm.f.ZIMS"           "sdbm.m.ZIMS"           "bm.f.amni"            
# [13] "bm.m.amni"             "bm.m.jm"               "bm.f.jm"               "bm.unk.anage"         
# [17] "bm.unk.amni"           "bm.unk.panth"          "asm.f.ZIMS"            "asm.m.ZIMS"           
# [21] "asm.f.jm"              "asm.m.jm"              "asm.f.anage"           "asm.m.anage"          
# [25] "asm.f.amni"            "asm.m.amni"            "asm.unk.amni"          "asm.f.datlife"        
# [29] "asm.m.datlife"         "asm.unk.datlife"       "asm.unk.panth"         "lcs.anage"            
# [33] "lcs.amni"              "lcs.panth"             "lc.y.anage"            "lc.y.amni"            
# [37] "lc.y.panth"            "ilibi.anage"           "ilibi.amni"            "ilibi.panth"          
# [41] "long.amni"             "maxlong.anage"         "maxlong.amni"          "maxlong.panth"   


lhMat <- matrix(NA, nPred, nvars, 
                dimnames = list(predDat$gbif.species, varLables))

# fill-in missing body mass:
for (iv in 1:nvars) {
  varCols <- sprintf("%s.%s", varNames[iv], dbs)
  varCols <- varCols[which(varCols %in% allPredNames)]
  temp <- apply(as.matrix(predDat[, varCols]), 1, mean, na.rm = TRUE)
  lhMat[, varLables[iv]] <- temp
}

# Sexual size dimorphism:
SSD <- (log(lhMat[, "BM_Male"]) - log(lhMat[, "BM_Female"]))

# ------------------ #
# 2) ALL PREDICTORS:
# ------------------ #

predDf <- data.frame(species = predDat$gbif.species, 
                     MS = predDat$matsys.socialmonogamy,
                     TESTES = predDat$testes.mass.combined, 
                     CARE = predDat$parental.care)

predDf <- cbind(predDf, lhMat, SSD = SSD)
if (taxa == "Mammalia") {
  predDf$MS[grep("UNDESCRIBED", predDf$MS)] <- NA
  predDf$MS[grep("EVIDENCE", predDf$MS)] <- "No"
}

# Weights for life expectancy:
# (Note: take the inverse if you used the variance in the estimate)
varex <- lifeExpBaSTA_new$Male_None_SE + lifeExpBaSTA_new$Male_Surgical_SE

weightsEx <- varex / max(varex, na.rm = TRUE)

#print(weightsEx)


# Reduce the life history data and keep only species present in the life expectancy data
complete_data <- merge(lifeExpBaSTA_new, predDf, by.x ="vertlife", 
                       by.y = "species", all.x = TRUE)

reduced_predDf <- complete_data[,44:ncol(complete_data)]
reduced_predDf <- cbind( complete_data[,1], reduced_predDf)
colnames(reduced_predDf)[1]<- "species"

# Add up all the relevant columns in one single fullDat file
fullDat <- cbind(reduced_predDf, respDat, weightsEx = weightsEx)


# Check to see if the species name are correct
fullDat2 <-fullDat
for (i in 1:nrow(fullDat2)) {
  if (fullDat2[i, "species"] %in% lifeExpBaSTA_new[i, "vertlife"]) {
    next
  } else {
    fullDat2$species[i] <- NA
  }
}

# Check for duplicates in the data
check_duplicates <- duplicated(fullDat$species)

# Print the duplicate names
print(fullDat$species[check_duplicates])

# Handle duplication 
duplicates <- rbind(fullDat[37,], fullDat[38,])
cervus <- as.data.frame(t(apply(duplicates[13:14], 2, mean, na.rm=TRUE)))
fullDat[37,13] <- cervus[1,1]
fullDat[37,14] <- cervus[1,2]
fullDat <- fullDat[-38, ] # delete duplicate


# Residuals from testis mass
resiDat <- fullDat[which(!(is.na(fullDat$TESTES)) & !(is.na(fullDat$BM_Male))), ]
mod.testis <- lm(log(resiDat$BM_Male) ~ log(resiDat$TESTES))
summary(mod.testis)
resiDat$TESTES_res <- mod.testis$residuals
resiDat <- resiDat[, c("species", "TESTES_res")]
fullDat <- merge(fullDat, resiDat, by = c("species"), all.x = T)


# ==== ADDITIONAL LIFE HISOTORY VARIABLES: ====
Species_TotalNew_18_April <- read_excel("/work/CatalinaThesis/analysis/data/tables/Species_TotalNew.18.May.xlsx")

fullDat <- merge(fullDat, Species_TotalNew_18_April, by = c("species"), all.x = T)
fullDat[fullDat == "NaN"] <- NA

for (i in 1:nrow(fullDat)){
  if (!is.na(fullDat$Parental[i])){
    if (fullDat$Parental[i] == "Female"){
      fullDat$Parental[i] <- "female-only"
    }else if (fullDat$Parental[i] == "Both"){
      fullDat$Parental[i] <- "biparental"
    }else{
      next
    }
  }
  
}


for (i in 1:nrow(fullDat)){
  if (is.na(fullDat$CARE[i])){
    val <- fullDat$Parental[i]
    fullDat$CARE[i] <- val
  }
}


names(fullDat)[16] <-"MatingSimplified"

polygamy <- c("polygynous", "polyandrous","polygamous", "polygynous, polygynandrous (promiscuous)",
              "monogamous, polygynous, cooperative breeder",
              "polyandrous, cooperative breeder", 
              "polygynous, polygynandrous (promiscuous), cooperative breeder",
              "polygynandrous (promiscuous), cooperative breeder", 
              "polygynous, cooperative breeder", 
              "monogamous, polyandrous, polygynous")

promiscuity <- c("polygynandrous (promiscuous)")

monogamy <- c("monogamous", 
              "monogamous, polygynous", 
              "monogamous, polyandrous, cooperative breeder", 
              "monogamous, cooperative breeder", 
              "monogamous, polyandrous")

fullDat$MatingSimplified <- NA

for (i in 1:nrow(fullDat)){
  if (!is.na(fullDat$MatingSystem[i])){
    if (fullDat$MatingSystem[i] %in% polygamy){
      fullDat$MatingSimplified[i] <- "polygamy"
    }else if (fullDat$MatingSystem[i] %in% promiscuity){
      fullDat$MatingSimplified[i] <- "promiscuity"
    }else if (fullDat$MatingSystem[i] %in% monogamy){
      fullDat$MatingSimplified[i] <- "monogamy"
    }
  }
}

# ==== CROSS-CHECKING DATA: ====
# Cross-checking the mating systems collected from https://animaldiversity.org/
# to the mating systems collected from https://animalia.bio/

Species_TotalNew_18_April <- read_excel("/work/CatalinaThesis/analysis/data/tables/Species_TotalNew.18.April.xlsx")
Species_simplified <- read_excel("CatalinaThesis/analysis/data/tables/Species_simplified.xlsx")
View(Species_TotalNew_18_April)

Species_simplified <- Species_simplified[-38, ]
species <- as.data.frame(fullDat$MatingSimplified)
species2 <- as.data.frame(fullDat$MatingSystem)
family <- as.data.frame(Species_simplified$Matingbehaviour)
family <- cbind(species, species2, family)
family

fullDat$MatingSimplified[20] <- "polygamy"
fullDat$MatingSimplified[28] <- "polygamy"
fullDat$MatingSimplified[33] <- "polygamy"
fullDat$MatingSimplified[34] <- "promiscuity"
fullDat$MatingSimplified[39] <- "monogamy"
fullDat$MatingSimplified[48] <- "polygamy"
fullDat$MatingSimplified[52] <- "promiscuity"  #article
fullDat$MatingSimplified[53] <- "polygamy"  #article
fullDat$MatingSimplified[55] <- "polygamy"
fullDat$MatingSimplified[71] <- "polygamy"
fullDat$MatingSimplified[72] <- "polygamy" #Mueller_et_al_Mating_system_Supplements-V
fullDat$MatingSimplified[77] <- "polygamy" #article
fullDat$MatingSimplified[85] <- "promiscuity" #article
fullDat$MatingSimplified[90] <- "promiscuity" #article
fullDat$MatingSimplified[92] <- "polygamy"
fullDat$MatingSimplified[96] <- "polygamy"
fullDat$MatingSimplified[97] <- "promiscuity"
fullDat$MatingSimplified[100] <- "promiscuity"
fullDat$MatingSimplified[102] <- "polygamy"
fullDat$MatingSimplified[103] <- "polygamy"  #article
fullDat$MatingSimplified[104] <- "polygamy"  #article
fullDat$MatingSimplified[106] <- "polygamy"  #article
fullDat$MatingSimplified[107] <- "polygamy"
fullDat$MatingSimplified[112] <- "polygamy" #Mueller_et_al_Mating_system_Supplements-V
fullDat$MatingSimplified[114] <- "polygamy"
fullDat$MatingSimplified[117] <- "promiscuity"
fullDat$MatingSimplified[117] <- "promiscuity" #article
fullDat$MatingSimplified[122] <- "polygamy"
fullDat$MatingSimplified[127] <- "polygamy" #web article
fullDat$MatingSimplified[128] <- "polygamy"
fullDat$MatingSimplified[139] <- "monogamy" # article
fullDat$MatingSimplified[143] <- "promiscuity"
fullDat$MatingSimplified[145] <- "polygamy"
fullDat$MatingSimplified[148] <- "promiscuity"
fullDat$MatingSimplified[150] <- "polygamy"
fullDat$MatingSimplified[151] <- "promiscuity"
fullDat$MatingSimplified[152] <- "monogamy"  #article
fullDat$MatingSimplified[161] <- "polygamy"
fullDat$MatingSimplified[162] <- "polygamy"
fullDat$MatingSimplified[164] <- "polygamy" # article
fullDat$MatingSimplified[167] <- "polygamy"
fullDat$MatingSimplified[172] <- "polygamy"



# Delete the outlier since it has a value of -40 as compared to the other species
fullDat <- fullDat[fullDat$species != "Pseudocheirus peregrinus",]

# Write the final data to an Excel file
write_xlsx(fullDat, path = "DataSet.xlsx")


# ==== EXPLORATORY DATA ANALYSIS: ====

# Check the values and the minimum and maximum values in the dataset 
min(fullDat$SSD, na.rm=TRUE)
max(fullDat$SSD, na.rm=TRUE)

min(fullDat$castration_difference, na.rm=TRUE)
max(fullDat$castration_difference, na.rm=TRUE)

min(fullDat$TESTES_res, na.rm=TRUE)
max(fullDat$TESTES_res, na.rm=TRUE)


min(fullDat$ASM_Male, na.rm=TRUE)
max(fullDat$ASM_Male, na.rm=TRUE)

min(fullDat$BM_Male, na.rm=TRUE)
max(fullDat$BM_Male, na.rm=TRUE)

min(log(fullDat$BM_Male), na.rm=TRUE)
max(log(fullDat$BM_Male), na.rm=TRUE)

mean_by_mating_system <- function(data, mating_system) {
  subset_data <- data[data$MatingSimplified == mating_system,]
  mean_castration_difference <- mean(subset_data$castration_difference, na.rm = TRUE)
  return(mean_castration_difference)
}

mean_monogamy <- mean_by_mating_system(fullDat, "monogamy")
mean_polygamy <- mean_by_mating_system(fullDat, "polygamy")
mean_promiscuity <- mean_by_mating_system(fullDat, "promiscuity")

mean_monogamy
mean_polygamy
mean_promiscuity

# Count occurrences of each level in CARE
care_counts <- table(fullDat$CARE)
print(care_counts)

# Count occurrences of each level in MatingSimplified
mating_counts <- table(fullDat$MatingSimplified)
print(mating_counts)



# Filter columns where both castration_difference and ASM_Male have values
filtered_ASM <- fullDat[complete.cases(fullDat$castration_difference, fullDat$ASM_Male), ]
filtered_BM <- fullDat[complete.cases(fullDat$castration_difference, fullDat$BM_Male), ]
filtered_SSD <- fullDat[complete.cases(fullDat$castration_difference, fullDat$SSD), ]
filtered_Testies <- fullDat[complete.cases(fullDat$castration_difference, fullDat$TESTES_res), ]
filtered_ASM$log <- log(filtered_ASM$ASM_Male)

# Count the number of positive and negative values in the filtered_ASM column
positive_count <- sum(filtered_ASM$log > 0)
negative_count <- sum(filtered_ASM$log < 0)

# Display the counts
print(paste("Positive count:", positive_count))
print(paste("Negative count:", negative_count))

# Look at log(filtered_BM$BM_Male)
filtered_BM$log <- log(filtered_BM$BM_Male)

# Count the number of positive and negative values in the filtered_SSD column
positive_count <- sum(filtered_SSD$SSD > 0)
negative_count <- sum(filtered_SSD$SSD < 0)
print(paste("Positive count:", positive_count))
print(paste("Negative count:", negative_count))


# Count the number of positive and negative values in the filtered_Testes column
positive_count <- sum(filtered_Testies$TESTES_res > 0)
negative_count <- sum(filtered_Testies$TESTES_res < 0)
print(paste("Positive count:", positive_count))
print(paste("Negative count:", negative_count))


families <- species_data_with_family
families <- families[, -c(4:43)]
families <- families[, -c(1:2)]
names(families)[names(families) == "vertlife"] <- "species"

families <- families[-38, ] # delete duplicate
families <- families[families$species != "Pseudocheirus peregrinus",]
families <- merge(fullDat, families, by = "species", all.x = TRUE)

# Extract species, families, and castration difference into a new dataframe
extracted_data <- data.frame(species = families$species, family = families$family, castration_difference = fullDat$castration_difference)
extracted_data <- extracted_data[extracted_data$castration_difference > 0, ]

# Count the number of species for each family
species_counts <- table(extracted_data$family)


families2 <- species_data_with_family
families2 <- families2[, -c(4:43)]
families2 <- families2[, -c(1:2)]
names(families2)[names(families2) == "vertlife"] <- "species"

families2 <- families2[-38, ] # delete duplicate
families2 <- families2[families2$species != "Pseudocheirus peregrinus",]
families2 <- merge(fullDat, families2, by = "species", all.x = TRUE)

# Extract species, families, and castration difference into a new dataframe
extracted_data2 <- data.frame(species = families2$species, family = families2$family, castration_difference = fullDat$castration_difference)
extracted_data2 <- extracted_data2[extracted_data2$castration_difference < -1, ]
# Count the number of species for each family
species_counts2 <- table(extracted_data2$family)

# Display the species counts
print(species_counts2)
barplot(species_counts2)

# ==== PHYLOGENETIC REGRESSIONS: ====

# Load phylogeny:
load(sprintf("/work/CatalinaThesis/analysis/data/rdata/trees/maxCredTree%s.RData", taxa)) 
phylo <- maxCred
phylo$node.label<-NULL # problems with error: Labels duplicated between tips and nodes in phylogeny 

# List with the life traits considere
vars <- c("castration_difference", "CARE", "log(BM_Male)", "log(ASM_Male)", "SSD", 
          "TESTES_res", "MatingSimplified")

N <- list(1,2)
COMB <- sapply(N, function(m) combn(unique(vars[2:7]), m))

COMB2 <- list()
k=0
for(i in seq(COMB)){
  tmp <- COMB[[i]]
  for(j in seq(ncol(tmp))){
    k <- k + 1
    COMB2[[k]] <- formula(paste("castration_difference", "~", paste(tmp[,j], collapse=" + ")))
  }
}
COMB2

# Run different models
number <- 1
for (i in 1:length(COMB2)){
  form = COMB2[[i]]
  print(form)
  name <- paste("out", number, sep = "")
  # wrap analysis in a tryCatch block to catch any errors
  analysis <- tryCatch({
    RunBayesPGLS(formula = form, data = fullDat, phylo = phylo, weights = "weightsEx", nsim = 6, ncpus = 6)
  }, error = function(e) {
    message("Error in analysis: ", e$message)
    return(NULL)
  })
  
  # check if analysis result is NULL or not
  if(!is.null(analysis)){
    assign(name, analysis)
  }
  
  # increment counter
  number <- number + 1
  print(number)
}

# ==== SAVE RESULTS: ====
# get a list of all objects in the current environment that start with "out"
obj_list <- ls(pattern = "^out")

# save the objects to a file
save(list = obj_list, file = "/work/CatalinaThesis/analysis/results/rdata/LifeExpDifferencesRegrAnalysesMasterThesis.RData")

# ==== IMPLEMENT AND SAVE TABLE: ====
rm(list = ls())
sysinf <- Sys.info()
load("/work/CatalinaThesis/analysis/results/rdata/LifeExpDifferencesRegrAnalysesMammaliaMasterThesis.RData")

list <- list()
o <- objects()
o <- o[grep("out", o)]
for(i in 1:length(o)){
  model <- paste0("out", i)
  if(exists(model)){
    form <- get(model)$form
    vars <- rownames(get(model)[[1]][1])
    mean <- get(model)[[1]][1]
    zeroCov <- get(model)[[1]][5]
    DIC <- round(get(model)$DIC[4], 2)
    SD <- get(model)[[1]][2]
    Lower <- get(model)[[1]][3]
    Upper <- get(model)[[1]][4]
    
    df <- data.frame(model = c(model, rep("", length(vars)-1)), 
                     form = c(form, rep("", length(vars)-1)), 
                     vars = vars, 
                     mean = round(mean, 2),
                     zeroCov = zeroCov,
                     SD = round(SD, 2),
                     Lower = round(Lower, 2),
                     Upper = round(Upper, 2),
                     DIC = c(DIC, rep("", length(vars)-1)))
                     
    list[[i]] <- df
  }
}
df <- do.call("rbind", list)
rownames(df) <- NULL

ft <- theme_vanilla(flextable(df))
# Adjust column widths
ft <- width(ft, j = "model", width = 0.5)
ft <- width(ft, j = "form", width = 0.5)
ft <- width(ft, j = "vars", width = 0.5)
ft <- width(ft, j = "Mean", width = 0.5)
ft <- width(ft, j = "SD", width = 0.5)
ft <- width(ft, j = "Lower", width = 0.5)
ft <- width(ft, j = "Upper", width = 0.5)
ft <- width(ft, j = "zeroCoverage", width = 0.5)
ft <- width(ft, j = "DIC", width = 0.5)
summary(out1)
head(df)

ft <- highlight(ft, j = "zeroCoverage", i = ~ zeroCoverage <= 0.05, color = "yellow")
ft <- autofit(ft)
ft <- fontsize(ft, size = 10)
ft


# Save flextable
save_as_docx(ft, path = "/work/CatalinaThesis/analysis/results/tables/LifeExpDifferencesRegrAnalysesMammaliaSpeciesMasterThesis.RData.docx")


# ==== PLOT DIAGNOSTIC PLOTS: ====

# Save the diagnostics plot
summary(out2)
pdf("Care.pdf")
plot(out2)
plot(out2, plot.type = "density")
plot(out2, plot.type = "diagnostics")
dev.off()

summary(out3)
pdf("BM_Male.pdf")
plot(out3)
plot(out3, plot.type = "density")
plot(out3, plot.type = "diagnostics")
dev.off()

summary(out4)
pdf("ASM_Male.pdf")
plot(out4)
plot(out4, plot.type = "density")
plot(out4, plot.type = "diagnostics")
dev.off()

summary(out5)
pdf("SSD.pdf")
plot(out5)
plot(out5, plot.type = "density")
plot(out5, plot.type = "diagnostics")
dev.off()

summary(out6)
pdf("Testes.pdf")
plot(out6)
plot(out6, plot.type = "density")
plot(out6, plot.type = "diagnostics")
dev.off()

summary(out7)
pdf("Mating systems.pdf")
plot(out7)
plot(out7, plot.type = "density")
plot(out7, plot.type = "diagnostics")
dev.off()

# ==== VISUALISATION PLOTS: ====

ASM_M <-fullDat[complete.cases(fullDat[c("castration_difference", "ASM_Male")]),]
ASM_plot <- ggplot(data = ASM_M, aes(x = log(ASM_Male), y = castration_difference, color = castration_difference < 0)) +
  geom_hline(yintercept=0, color = "grey")+
  geom_abline(intercept = out4$coefficients$Mean[1], slope = out4$coefficients$Mean[2], color = "black", linetype = "dashed")+
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("red", "#0099f9"), name = "Lifespan Difference < 0") +
  guides(size = guide_legend(title = "Age at sexual maturity"), color = guide_legend(title = "Lifespan increase"))+
  labs(
    title = "Effect of Age at Sexual Maturity on Lifespan Difference",
    x = "Age at sexual maturity",
    y = "Lifespan Difference",
    caption = paste("----Intercept =", round(out4$coefficients$Mean[1], 3), "\n",
                    "       Slope =", round(out4$coefficients$Mean[2], 3)))+
  theme(
    plot.title = element_text(color = "#0099f9", size = 10, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 9, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "right",
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0, size=12))
print(ASM_plot)
ggsave("ASM_Male.pdf", ASM_plot)
summary(out4)

Testies <-fullDat[complete.cases(fullDat[c("castration_difference", "TESTES_res")]),]
summary(out6)
tes <- ggplot(data = Testies, aes(x = TESTES_res, y = castration_difference, color = castration_difference < 0)) +
  geom_hline(yintercept=0, color = "grey")+
  geom_point(alpha = 0.6)+
  geom_abline(intercept = out6$coefficients$Mean[1], slope = out6$coefficients$Mean[2], color = "black", linetype = "dashed")+
  scale_color_manual(values = c("red", "#0099f9"), name = "Lifespan Difference < 0") +
  guides(size = guide_legend(title = "Testicle Size"), color = guide_legend(title = "Lifespan increase"))+
  labs(
    title = "Effect of Testicles' residuals on Lifespan Difference",
    x = "Testicles' residuals",
    y = "Lifespan Difference",
    caption = paste("----Intercept =", round(out6$coefficients$Mean[1], 3), "\n",
                    "       Slope =", round(out6$coefficients$Mean[2], 3)))+
  theme(
    plot.title = element_text(color = "#0099f9", size = 11, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "right",
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0, size=12))
print(tes)
ggsave("Testicles.pdf", tes)

BM <-fullDat[complete.cases(fullDat[c("castration_difference", "BM_Male")]),]
summary(out3)
BM2 <- ggplot(data = BM, aes(x = log(BM_Male), y = castration_difference, color = castration_difference < 0)) +
  geom_hline(yintercept=0, color = "grey")+
  geom_point(alpha = 0.6) +
  geom_abline(intercept = out3$coefficients$Mean[1], slope = out3$coefficients$Mean[2], color = "black", linetype = "dashed")+
  scale_color_manual(values = c("red", "#0099f9"), name = "Castration Difference < 0") +
  guides(size = guide_legend(title = "Body Mass Males"), color = guide_legend(title = "Lifespan increase"))+
  labs(
    title = "Effect of Body Mass on Lifespan Difference",
    x = "Body Mass Males",
    y = "Lifespan Difference",
    caption = paste("----Intercept =", round(out3$coefficients$Mean[1], 3), "\n",
                    "       Slope =", round(out3$coefficients$Mean[2], 3)))+
  theme(
    plot.title = element_text(color = "#0099f9", size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "right",
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0, size=12))
print(BM2)
ggsave("BM.pdf", BM2)

SSD <-fullDat[complete.cases(fullDat[c("castration_difference", "SSD")]),]
summary(out5)
SSD2 <- ggplot(data = SSD, aes(x = SSD, y = castration_difference, color = castration_difference < 0)) +
  geom_hline(yintercept=0, color = "grey")+
  geom_point(alpha = 0.6)+
  geom_abline(intercept = out5$coefficients$Mean[1], slope = out5$coefficients$Mean[2], color = "black", linetype = "dashed")+
  scale_color_manual(values = c("red", "#0099f9"), name = "Lifespan Difference < 0") +
  guides(size = guide_legend(title = "Sexual Size Dimorphism"), color = guide_legend(title = "Lifespan increase"))+
  labs(
    title = "Effect of Sexual Size Dimorphism on Lifespan Difference",
    x = "Sexual Size Dimorphism",
    y = "Lifespan Difference",
    caption = paste("----Intercept =", round(out5$coefficients$Mean[1], 3), "\n",
                    "       Slope =", round(out5$coefficients$Mean[2], 3)))+
  theme(
    plot.title = element_text(color = "#0099f9", size = 10, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "right",
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0, size=12))
print(SSD2)
ggsave("SSD.pdf", SSD2)

mating <-fullDat[complete.cases(fullDat[c("castration_difference", "MatingSimplified")]),]
summary(out7)
mating3 <- ggplot(data = mating, aes(x = MatingSimplified, y = castration_difference)) +
  geom_hline(yintercept=0, color = "grey") +
  geom_boxplot(aes(color = MatingSimplified), width = 0.3) +
  scale_color_manual(values = c("red", "#0099f9", "green"), name = "Mating System") +
  labs(
    title = "Effect of Mating System on Lifespan Difference",
    x = "Mating System",
    y = "Lifespan Difference",
    caption = paste("----Intercept =", round(out7$coefficients$Mean[1], 3), "\n",
                    "       Polygamy =", round(out7$coefficients$Mean[2], 3), "\n",
                    "       Promiscuity =", round(out7$coefficients$Mean[3], 3)))+
  theme(
    plot.title = element_text(color = "#0099f9", size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "right",
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0, size=12))
print(mating3)
ggsave("Mating2.pdf", mating3)

Care <-fullDat[complete.cases(fullDat[c("castration_difference", "CARE")]),]
summary(out2)
Care3 <- ggplot(data = Care, aes(x = CARE, y = castration_difference)) +
  geom_hline(yintercept=0, color = "grey") +
  geom_boxplot(aes(color = CARE), width = 0.3) +
  scale_color_manual(values = c("red", "#0099f9"), name = "Parental Care") +
  labs(
    title = "Effect of Parental Care on Lifespan Difference",
    x = "Parental Care",
    y = "Lifespan Difference",
    caption = paste("----Intercept =", round(out2$coefficients$Mean[1], 3), "\n",
                    "       Female only =", round(out2$coefficients$Mean[2], 3)))+
  theme(
    plot.title = element_text(color = "#0099f9", size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "right",
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0, size=12))
print(Care3)
ggsave("Care2.pdf", Care3)
