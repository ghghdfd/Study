 ####################################################conect mysql in r##########################################################
#install.packages('RMySQL')
library(RMySQL)
drv <- dbDriver("MySQL")
con <- dbConnect(MySQL(), 
                 dbname = "tcga", 
                 user = "hwanghou", 
                 host = "103.22.220.149", 
                 port = 13306,
                 password = "yonsei2020!")

### tumor samples
r_sql <- c("SELECT `icgc_donor_id`,`project_code`,`submitted_sample_id`,`gene_stable_id`,`normalized_expression_level`, `platform` FROM `gene_expression` WHERE `submitted_sample_id`like'%-01A%'")
platform <- dbGetQuery(con, r_sql) #running 20min
nrow(platform)

head(platform)
table(platform$platform)

#Affymetrix HT Human Genome U133A     Affymetrix Human Exon 1.0 ST
#                        13452032                         17615457
#Agilent 244K Custom Gene Express           Illumina GA sequencing
#                        30376281                          1204467
#Illumina HiSeq
#     102325608


##### exclude "ENSG"
library(dplyr)
ENSG_na<-platform %>% 
  filter(!grepl('ENSG', gene_stable_id))
head(ENSG_na, 20)
ENSG_na[grep("ENSG", ENSG_na$gene_stable_id),]
nrow(ENSG_na) #103530075
head(ENSG_na)
table(ENSG_na$platform)
#Illumina GA sequencing         Illumina HiSeq
#               1204467              102325608


#####  include "ENSG"
ENSG_y<-platform[grep("ENSG", platform$gene_stable_id),]
ENSG_y %>% 
  filter(!grepl('ENSG', gene_stable_id))
nrow(ENSG_y) #61443770
head(ENSG_y)
table(ENSG_y$platform)
#Affymetrix HT Human Genome U133A     Affymetrix Human Exon 1.0 ST
#                        13452032                         17615457
#Agilent 244K Custom Gene Express
#                        30376281



########################################################################################################################
library(RMySQL)
drv <- dbDriver("MySQL")
con <- dbConnect(MySQL(), 
                 dbname = "tcga", 
                 user = "hwanghou", 
                 host = "103.22.220.149", 
                 port = 13306,
                 password = "yonsei2020!")


r_sql <- c("SELECT `icgc_donor_id`,`project_code`,`submitted_sample_id`,`gene_stable_id`,`normalized_expression_level`, `platform` FROM `gene_expression` WHERE `submitted_sample_id`like'%-11A%'")
data_ <- dbGetQuery(con, r_sql) #running 20min
data_SolidA<-data_
r_sql <- c("SELECT `icgc_donor_id`,`project_code`,`submitted_sample_id`,`gene_stable_id`,`normalized_expression_level`, `platform` FROM `gene_expression` WHERE `submitted_sample_id`like'%-11B%'")
data_ <- dbGetQuery(con, r_sql) #running 20min
data_SolidB<-data_
r_sql <- c("SELECT `icgc_donor_id`,`project_code`,`submitted_sample_id`,`gene_stable_id`,`normalized_expression_level`, `platform` FROM `gene_expression` WHERE `submitted_sample_id`like'%-11C%'")
data_ <- dbGetQuery(con, r_sql) #running 20min
data_SolidC<-data_

data_Solid<-rbind(data_SolidA, data_SolidB, data_SolidC)
head(data_Solid)

table(data_Solid$platform)



##### exclude "ENSG"
library(dplyr)
ENSG_na<-data_Solid %>% 
  filter(!grepl('ENSG', gene_stable_id))
head(ENSG_na, 20)
ENSG_na[grep("ENSG", ENSG_na$gene_stable_id),]
nrow(ENSG_na) #103530075
head(ENSG_na)
table(ENSG_na$platform)


ENSG_y<-data_Solid[grep("ENSG", data_Solid$gene_stable_id),]
ENSG_y %>% 
  filter(!grepl('ENSG', gene_stable_id))
nrow(ENSG_y) #61443770
head(ENSG_y)
table(ENSG_y$platform)


##############################################################################################################
####################################     features  ###########################################################
##############################################################################################################

##########################
# clinical data
clinical<-read.csv("/home/hwanghou/20200311/clinical (2).csv")
head(clinical)
colnames(clinical)
table(clinical$tumour_stage)
clinical2<-subset(clinical, select=c(icgc_donor_id , submitted_specimen_id, project_code, donor_sex, donor_age_at_diagnosis))

clinical2<-clinical[c("icgc_donor_id" , "submitted_specimen_id", "project_code", "donor_sex", "donor_age_at_diagnosis")]
head(clinical2)
clinical2_pyrhon<-clinical2
colnames(clinical2_pyrhon) = c("icgc_donor_id" ,"submitted_donor_id", "project_code_rename", "donor_sex_rename", "donor_age_at_diagnosis_rename")
head(clinical2_pyrhon)

####################################
## stage data
stage <- read.table('/home/hwanghou/20200312/Final_clinical.tsv',header=TRUE,sep='\t')
table(stage$Stage)
stage$bcr_patient_barcode

stage <- stage[,c(2,4,5,7,9)]
nrow(stage)
head(stage)
stage<-unique(stage)
names(stage)[1] <- c('bcr_patient_barcode')
stage_python<-stage
head(stage_python)
colnames(stage_python) = c( "submitted_sample_id", "Stage_rename", "Status_rename")
head(stage_python)
stage_python <- unique(stage_python)
sapply(stage_python, function(x) sum(is.na(x))) # count each columns 


table(stage_python$Stage)



table(stage_python$Stage_rename)

clinical2_pyrhon<-merge(x=clinical2_pyrhon, y=stage_python, by = c("submitted_sample_id"), all.x = TRUE)

#####################################
## match samples feature 'patient'

tcga_patient_id <- data[,c('icgc_donor_id','submitted_sample_id')]
tcga_patient_id<-unique(tcga_patient_id)
tcga_patient_id$submitted_sample_id <- substr(tcga_patient_id$submitted_sample_id,1,12)
names(tcga_patient_id)[2] <- c('bcr_patient_barcode')
write.csv(tcga_patient_id, file = "/home/hwanghou/20200311/tcgadatabase_patient_feature.csv", row.names = FALSE)
tcga_patient_id<-read.csv("/home/hwanghou/20200311/tcgadatabase_patient_feature.csv")
head(tcga_patient_id)
str(unique(tcga_patient_id$icgc_donor_id))
str(unique(tcga_patient_id$bcr_patient_barcode))



tcga_patient<-read.csv("/home/hwanghou/20200311/tcga_patient.csv")
head(tcga_patient)
tcga_patient2<-subset(tcga_patient, select=c(bcr_patient_barcode , race))
head(tcga_patient2)



tcga_patient3<-merge(x=tcga_patient_id, y=tcga_patient2, by = "bcr_patient_barcode", all.x = TRUE)
str(unique(tcga_patient3$bcr_patient_barcode))

tcga_patient3 <- unique(tcga_patient3)
table(tcga_patient3$race)
sapply(tcga_patient3, function(x) sum(is.na(x))) # count each columns 
# bcr_patient_barcode       icgc_donor_id submitted_sample_id                race 
#                   0                   0                   0                 689 
head(clinical2_pyrhon)
clinical2_pyrhon <- clinical2_pyrhon[,-c(2)]
head(tcga_patient3)

tcgabatabase_feature<-merge(x=clinical2_pyrhon, y=tcga_patient3, by = c("icgc_donor_id"), all.x = TRUE)
tcgabatabase_feature <- unique(tcgabatabase_feature)
sapply(tcgabatabase_feature, function(x) sum(is.na(x))) # count each columns 
#icgc_donor_id            submitted_donor_id           project_code_rename              donor_sex_rename donor_age_at_diagnosis_rename 
#            0                             0                             0                             0                             0 
#bcr_patient_barcode           submitted_sample_id                          race 
#               5207                          5207                          7301 


tcgabatabase_feature_python<-tcgabatabase_feature
head(tcgabatabase_feature_python)
nrow(tcgabatabase_feature_python)
tcgabatabase_feature_python<-unique(tcgabatabase_feature_python)

#age grouping
str(tcgabatabase_feature_python)
tcgabatabase_feature_python$donor_age_at_diagnosis_rename<-as.numeric(tcgabatabase_feature_python$donor_age_at_diagnosis_rename)
tcgabatabase_feature_python$donor_age_at_diagnosis<-tcgabatabase_feature_python$donor_age_at_diagnosis_rename
table(tcgabatabase_feature_python$donor_age_at_diagnosis)
tcgabatabase_feature_python[tcgabatabase_feature_python$donor_age_at_diagnosis > 0 & tcgabatabase_feature_python$donor_age_at_diagnosis < 10, "donor_age_at_diagnosis"] = 1
tcgabatabase_feature_python[tcgabatabase_feature_python$donor_age_at_diagnosis > 9 & tcgabatabase_feature_python$donor_age_at_diagnosis < 20, "donor_age_at_diagnosis"] = 10
tcgabatabase_feature_python[tcgabatabase_feature_python$donor_age_at_diagnosis > 19 & tcgabatabase_feature_python$donor_age_at_diagnosis < 30, "donor_age_at_diagnosis"] = 20
tcgabatabase_feature_python[tcgabatabase_feature_python$donor_age_at_diagnosis > 29 & tcgabatabase_feature_python$donor_age_at_diagnosis < 40, "donor_age_at_diagnosis"] = 30
tcgabatabase_feature_python[tcgabatabase_feature_python$donor_age_at_diagnosis > 39 & tcgabatabase_feature_python$donor_age_at_diagnosis < 50, "donor_age_at_diagnosis"] = 40
tcgabatabase_feature_python[tcgabatabase_feature_python$donor_age_at_diagnosis > 49 & tcgabatabase_feature_python$donor_age_at_diagnosis < 60, "donor_age_at_diagnosis"] = 50
tcgabatabase_feature_python[tcgabatabase_feature_python$donor_age_at_diagnosis > 59 & tcgabatabase_feature_python$donor_age_at_diagnosis < 70, "donor_age_at_diagnosis"] = 60
tcgabatabase_feature_python[tcgabatabase_feature_python$donor_age_at_diagnosis > 69 & tcgabatabase_feature_python$donor_age_at_diagnosis < 80, "donor_age_at_diagnosis"] = 70
tcgabatabase_feature_python[tcgabatabase_feature_python$donor_age_at_diagnosis > 79 & tcgabatabase_feature_python$donor_age_at_diagnosis < 90, "donor_age_at_diagnosis"] = 80
tcgabatabase_feature_python[tcgabatabase_feature_python$donor_age_at_diagnosis > 89 & tcgabatabase_feature_python$donor_age_at_diagnosis < 100, "donor_age_at_diagnosis"] = 90
table(tcgabatabase_feature_python$donor_age_at_diagnosis)
#  0    1   10   20   30   40   50   60   70   80   90 
#418  556  160  249  667 1670 3029 3499 2566  688   39

table(tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("BLCA-US", 0, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("BRCA-US", 1, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("CESC-US", 2, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("COAD-US", 3, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("GBM-US", 4, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("HNSC-US", 5, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("KIRC-US", 6, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("KIRP-US", 7, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("LGG-US", 8, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("LIHC-US", 9, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("LUAD-US", 10, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("LUSC-US", 11, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("OV-US", 12, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("PAAD-US", 13, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("PRAD-US", 14, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("READ-US", 15, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("SKCM-US", 16, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("STAD-US", 17, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("THCA-US", 18, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("UCEC-US", 19, tcgabatabase_feature_python$project_code_rename)

tcgabatabase_feature_python$project_code_rename<-gsub("BOCA-UK", 20, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("BRCA-UK", 21, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("CLLE-ES", 22, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("CMDI-UK", 23, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("EOPC-DE", 24, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("ESAD-UK", 25, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("LAML-US", 26, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("LICA-FR", 27, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("LINC-JP", 28, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("LIRI-JP", 29, tcgabatabase_feature_python$project_code_rename)

tcgabatabase_feature_python$project_code_rename<-gsub("MALY-DE", 30, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("NBL-US", 31, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("ORCA-IN", 32, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("OV-AU", 33, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("PACA-AU", 34, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("PACA-CA", 35, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("PAEN-AU", 36, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("PBCA-DE", 37, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("PRAD-CA", 38, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("RECA-CN", 39, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("RECA-EU", 40, tcgabatabase_feature_python$project_code_rename)
tcgabatabase_feature_python$project_code_rename<-gsub("THCA-SA", 41, tcgabatabase_feature_python$project_code_rename)
table(tcgabatabase_feature_python$project_code_rename)
#   0    1   10   11   12   13   14   15   16   17   18   19    2   20   21   22   23   24   25   26   27   28   29    3   30   31   32   33   34   35   36   37   38   39    4   40   41    5    6    7 
# 202 1098  517  697 2013   59  213  177  321  325  548  515   85   69  141  271  129   17   16  200  153  244  158  481   55  389   50   93  373  142   22  306   10   10 1664  167   15  407  579  166 
#   8    9 
# 271  173 

table(tcgabatabase_feature_python$donor_sex_rename)
tcgabatabase_feature_python$donor_sex_rename<-gsub("female", "0", tcgabatabase_feature_python$donor_sex_rename)
tcgabatabase_feature_python$donor_sex_rename<-gsub("male", "1", tcgabatabase_feature_python$donor_sex_rename)
table(tcgabatabase_feature_python$donor_sex_rename)
# na   0    1 
#161 7629 5751 
tcgabatabase_feature_python$race_rename<-tcgabatabase_feature_python$race
head(tcgabatabase_feature_python)
table(tcgabatabase_feature_python$race_rename)
tcgabatabase_feature_python$race_rename<-gsub("AMERICAN INDIAN OR ALASKA NATIVE", "0", tcgabatabase_feature_python$race_rename)
tcgabatabase_feature_python$race_rename<-gsub("ASIAN", "1", tcgabatabase_feature_python$race_rename)
tcgabatabase_feature_python$race_rename<-gsub("BLACK OR AFRICAN AMERICAN", "2", tcgabatabase_feature_python$race_rename)
tcgabatabase_feature_python$race_rename<-gsub("NATIVE HAWAIIAN OR OTHER PACIFIC", "3", tcgabatabase_feature_python$race_rename)
tcgabatabase_feature_python$race_rename<-gsub("WHITE", "4", tcgabatabase_feature_python$race_rename)
table(tcgabatabase_feature_python$race_rename)

table(tcgabatabase_feature_python$donor_age_at_diagnosis)
tcgabatabase_feature_python$donor_age_at_diagnosis<-gsub("10", "0", tcgabatabase_feature_python$donor_age_at_diagnosis)
tcgabatabase_feature_python$donor_age_at_diagnosis<-gsub("20", "1", tcgabatabase_feature_python$donor_age_at_diagnosis)
tcgabatabase_feature_python$donor_age_at_diagnosis<-gsub("30", "2", tcgabatabase_feature_python$donor_age_at_diagnosis)
tcgabatabase_feature_python$donor_age_at_diagnosis<-gsub("40", "3", tcgabatabase_feature_python$donor_age_at_diagnosis)
tcgabatabase_feature_python$donor_age_at_diagnosis<-gsub("50", "4", tcgabatabase_feature_python$donor_age_at_diagnosis)
tcgabatabase_feature_python$donor_age_at_diagnosis<-gsub("60", "5", tcgabatabase_feature_python$donor_age_at_diagnosis)
tcgabatabase_feature_python$donor_age_at_diagnosis<-gsub("70", "6", tcgabatabase_feature_python$donor_age_at_diagnosis)
tcgabatabase_feature_python$donor_age_at_diagnosis<-gsub("80", "7", tcgabatabase_feature_python$donor_age_at_diagnosis)
tcgabatabase_feature_python$donor_age_at_diagnosis<-gsub("90", "8", tcgabatabase_feature_python$donor_age_at_diagnosis)
table(tcgabatabase_feature_python$donor_age_at_diagnosis)
tcgabatabase_feature_python <- tcgabatabase_feature_python[,-c(4,6)]
head(tcgabatabase_feature_python)



table(tcgabatabase_feature_python$Stage_rename)
tcgabatabase_feature_python$Stage_rename<-gsub("Stage IV", "3", tcgabatabase_feature_python$Stage_rename)
tcgabatabase_feature_python$Stage_rename<-gsub("Stage III", "2", tcgabatabase_feature_python$Stage_rename)
tcgabatabase_feature_python$Stage_rename<-gsub("Stage II", "1", tcgabatabase_feature_python$Stage_rename)
tcgabatabase_feature_python$Stage_rename<-gsub("Stage I", "0", tcgabatabase_feature_python$Stage_rename)
table(tcgabatabase_feature_python$Stage_rename)


head(tcgabatabase_feature_python)

tcgabatabase_feature_python_rename<-tcgabatabase_feature_python[-c(2, 3, 7, 8, 10)]
head(tcgabatabase_feature_python_rename)
head(tcgabatabase_feature)
colnames(tcgabatabase_feature) = c("icgc_donor_id" ,"submitted_donor_id", "project_code", "donor_sex", "donor_age_at_diagnosis", "bcr_patient_barcode", "submitted_sample_id", "race")

tcgabatabase_feature2<-merge(x=tcgabatabase_feature, y=tcgabatabase_feature_python_rename, by = c("submitted_sample_id"), all.x = TRUE)
head(tcgabatabase_feature2)
write.csv(tcgabatabase_feature2, file = "/home/solbi/20190118/tcgadatabase_patient_feature.csv", row.names = FALSE)


##### normal samplels table  ####
#all
all_4838<-read.csv("/home/solbi/20190118/tcga_primary_dataset_4838_geneexpres_dataset_4featuredata.csv")

all_feature<-all_4838[,c(1:14)]
head(all_feature)
sapply(all_feature, function(x) sum(is.na(x))) # count each columns 

table(all_feature$donor_age_at_diagnosis)
all_feature[all_feature$donor_age_at_diagnosis > 9 & all_feature$donor_age_at_diagnosis < 20, "donor_age_at_diagnosis"] = 10
all_feature[all_feature$donor_age_at_diagnosis > 19 & all_feature$donor_age_at_diagnosis < 30, "donor_age_at_diagnosis"] = 20
all_feature[all_feature$donor_age_at_diagnosis > 29 & all_feature$donor_age_at_diagnosis < 40, "donor_age_at_diagnosis"] = 30
all_feature[all_feature$donor_age_at_diagnosis > 39 & all_feature$donor_age_at_diagnosis < 50, "donor_age_at_diagnosis"] = 40
all_feature[all_feature$donor_age_at_diagnosis > 49 & all_feature$donor_age_at_diagnosis < 60, "donor_age_at_diagnosis"] = 50
all_feature[all_feature$donor_age_at_diagnosis > 59 & all_feature$donor_age_at_diagnosis < 70, "donor_age_at_diagnosis"] = 60
all_feature[all_feature$donor_age_at_diagnosis > 69 & all_feature$donor_age_at_diagnosis < 80, "donor_age_at_diagnosis"] = 70
all_feature[all_feature$donor_age_at_diagnosis > 79 & all_feature$donor_age_at_diagnosis < 90, "donor_age_at_diagnosis"] = 80
all_feature[all_feature$donor_age_at_diagnosis > 89 & all_feature$donor_age_at_diagnosis < 100, "donor_age_at_diagnosis"] = 90
table(all_feature$donor_age_at_diagnosis)

a1<-table(all_feature$donor_age_at_diagnosis, all_feature$donor_sex)
a2<-table(all_feature$race, all_feature$donor_sex)
a3<-table(all_feature$project_code, all_feature$donor_sex)
a4<-table(all_feature$Stage_rename, all_feature$donor_sex)
a<-rbind(a1, a2, a3, a4)
a

b1<-prop.table(ftable(donor_sex ~ donor_age_at_diagnosis, all_feature), margin = 2)*100
b2<-prop.table(ftable(donor_sex ~ race, all_feature), margin = 2)*100
b3<-prop.table(ftable(donor_sex ~ project_code, all_feature), margin = 2)*100
b4<-prop.table(ftable(donor_sex ~ Stage_rename, all_feature), margin = 2)*100
b<-rbind(b1, b2, b3, b4)
head(b)
colnames(b) = c("female", "male")
head(b)
b<-round(b,  digits =2)
b
ab<-cbind(a, b)
ab

write.csv(ab, file = "/home/solbi/20190118/cancersample/table1.csv", row.names = FALSE)


#BRCA
BRCA<-read.csv("/home/solbi/20190118/cancersample/tcga_primary_dataset_3782_geneexpres_dataset_stage_featuredata_BRCA-US_python.csv")

BRCA_feature<-BRCA[,c(1:14)]
head(BRCA_feature)
sapply(BRCA_feature, function(x) sum(is.na(x))) # count each columns 

table(BRCA_feature$donor_age_at_diagnosis)
BRCA_feature[BRCA_feature$donor_age_at_diagnosis > 0 & BRCA_feature$donor_age_at_diagnosis < 10, "donor_age_at_diagnosis"] = 1
BRCA_feature[BRCA_feature$donor_age_at_diagnosis > 9 & BRCA_feature$donor_age_at_diagnosis < 20, "donor_age_at_diagnosis"] = 10
BRCA_feature[BRCA_feature$donor_age_at_diagnosis > 19 & BRCA_feature$donor_age_at_diagnosis < 30, "donor_age_at_diagnosis"] = 20
BRCA_feature[BRCA_feature$donor_age_at_diagnosis > 29 & BRCA_feature$donor_age_at_diagnosis < 40, "donor_age_at_diagnosis"] = 30
BRCA_feature[BRCA_feature$donor_age_at_diagnosis > 39 & BRCA_feature$donor_age_at_diagnosis < 50, "donor_age_at_diagnosis"] = 40
BRCA_feature[BRCA_feature$donor_age_at_diagnosis > 49 & BRCA_feature$donor_age_at_diagnosis < 60, "donor_age_at_diagnosis"] = 50
BRCA_feature[BRCA_feature$donor_age_at_diagnosis > 59 & BRCA_feature$donor_age_at_diagnosis < 70, "donor_age_at_diagnosis"] = 60
BRCA_feature[BRCA_feature$donor_age_at_diagnosis > 69 & BRCA_feature$donor_age_at_diagnosis < 80, "donor_age_at_diagnosis"] = 70
BRCA_feature[BRCA_feature$donor_age_at_diagnosis > 79 & BRCA_feature$donor_age_at_diagnosis < 90, "donor_age_at_diagnosis"] = 80
BRCA_feature[BRCA_feature$donor_age_at_diagnosis > 89 & BRCA_feature$donor_age_at_diagnosis < 100, "donor_age_at_diagnosis"] = 90
table(BRCA_feature$donor_age_at_diagnosis)



a1<-as.data.table(table(BRCA_feature$donor_sex))
a2<-as.data.table(table(BRCA_feature$donor_age_at_diagnosis))
a3<-as.data.table(table(BRCA_feature$race))
a4<-as.data.table(table(BRCA_feature$Stage))
BRCA_feature_table<-rbind(a1, a2, a3, a4)
BRCA_feature_table

#KIRC
KIRC<-read.csv("/home/solbi/20190118/cancersample/tcga_primary_dataset_3782_geneexpres_dataset_stage_featuredata_KIRC-US_python.csv")

KIRC_feature<-KIRC[,c(1:14)]
head(KIRC_feature)
sapply(KIRC_feature, function(x) sum(is.na(x))) # count each columns 

table(KIRC_feature$donor_age_at_diagnosis)
KIRC_feature[KIRC_feature$donor_age_at_diagnosis > 0 & KIRC_feature$donor_age_at_diagnosis < 10, "donor_age_at_diagnosis"] = 1
KIRC_feature[KIRC_feature$donor_age_at_diagnosis > 9 & KIRC_feature$donor_age_at_diagnosis < 20, "donor_age_at_diagnosis"] = 10
KIRC_feature[KIRC_feature$donor_age_at_diagnosis > 19 & KIRC_feature$donor_age_at_diagnosis < 30, "donor_age_at_diagnosis"] = 20
KIRC_feature[KIRC_feature$donor_age_at_diagnosis > 29 & KIRC_feature$donor_age_at_diagnosis < 40, "donor_age_at_diagnosis"] = 30
KIRC_feature[KIRC_feature$donor_age_at_diagnosis > 39 & KIRC_feature$donor_age_at_diagnosis < 50, "donor_age_at_diagnosis"] = 40
KIRC_feature[KIRC_feature$donor_age_at_diagnosis > 49 & KIRC_feature$donor_age_at_diagnosis < 60, "donor_age_at_diagnosis"] = 50
KIRC_feature[KIRC_feature$donor_age_at_diagnosis > 59 & KIRC_feature$donor_age_at_diagnosis < 70, "donor_age_at_diagnosis"] = 60
KIRC_feature[KIRC_feature$donor_age_at_diagnosis > 69 & KIRC_feature$donor_age_at_diagnosis < 80, "donor_age_at_diagnosis"] = 70
KIRC_feature[KIRC_feature$donor_age_at_diagnosis > 79 & KIRC_feature$donor_age_at_diagnosis < 90, "donor_age_at_diagnosis"] = 80
KIRC_feature[KIRC_feature$donor_age_at_diagnosis > 89 & KIRC_feature$donor_age_at_diagnosis < 100, "donor_age_at_diagnosis"] = 90
table(KIRC_feature$donor_age_at_diagnosis)



a1<-as.data.table(table(KIRC_feature$donor_sex))
a2<-as.data.table(table(KIRC_feature$donor_age_at_diagnosis))
a3<-as.data.table(table(KIRC_feature$race))
a4<-as.data.table(table(KIRC_feature$Stage))
KIRC_feature_table<-rbind(a1, a2, a3, a4)
KIRC_feature_table

#UCEC
UCEC<-read.csv("/home/solbi/20190118/cancersample/tcga_primary_dataset_3782_geneexpres_dataset_stage_featuredata_UCEC-US_python.csv")

UCEC_feature<-UCEC[,c(1:14)]
head(UCEC_feature)
sapply(UCEC_feature, function(x) sum(is.na(x))) # count each columns 

table(UCEC_feature$donor_age_at_diagnosis)
UCEC_feature[UCEC_feature$donor_age_at_diagnosis > 0 & UCEC_feature$donor_age_at_diagnosis < 10, "donor_age_at_diagnosis"] = 1
UCEC_feature[UCEC_feature$donor_age_at_diagnosis > 9 & UCEC_feature$donor_age_at_diagnosis < 20, "donor_age_at_diagnosis"] = 10
UCEC_feature[UCEC_feature$donor_age_at_diagnosis > 19 & UCEC_feature$donor_age_at_diagnosis < 30, "donor_age_at_diagnosis"] = 20
UCEC_feature[UCEC_feature$donor_age_at_diagnosis > 29 & UCEC_feature$donor_age_at_diagnosis < 40, "donor_age_at_diagnosis"] = 30
UCEC_feature[UCEC_feature$donor_age_at_diagnosis > 39 & UCEC_feature$donor_age_at_diagnosis < 50, "donor_age_at_diagnosis"] = 40
UCEC_feature[UCEC_feature$donor_age_at_diagnosis > 49 & UCEC_feature$donor_age_at_diagnosis < 60, "donor_age_at_diagnosis"] = 50
UCEC_feature[UCEC_feature$donor_age_at_diagnosis > 59 & UCEC_feature$donor_age_at_diagnosis < 70, "donor_age_at_diagnosis"] = 60
UCEC_feature[UCEC_feature$donor_age_at_diagnosis > 69 & UCEC_feature$donor_age_at_diagnosis < 80, "donor_age_at_diagnosis"] = 70
UCEC_feature[UCEC_feature$donor_age_at_diagnosis > 79 & UCEC_feature$donor_age_at_diagnosis < 90, "donor_age_at_diagnosis"] = 80
UCEC_feature[UCEC_feature$donor_age_at_diagnosis > 89 & UCEC_feature$donor_age_at_diagnosis < 100, "donor_age_at_diagnosis"] = 90
table(UCEC_feature$donor_age_at_diagnosis)



a1<-as.data.table(table(UCEC_feature$donor_sex))
a2<-as.data.table(table(UCEC_feature$donor_age_at_diagnosis))
a3<-as.data.table(table(UCEC_feature$race))
a4<-as.data.table(table(UCEC_feature$Stage))
UCEC_feature_table<-rbind(a1, a2, a3, a4)
UCEC_feature_table

#LUAD
LUAD<-read.csv("/home/solbi/20190118/cancersample/tcga_primary_dataset_3782_geneexpres_dataset_stage_featuredata_LUAD-US_python.csv")

LUAD_feature<-LUAD[,c(1:14)]
head(LUAD_feature)
sapply(LUAD_feature, function(x) sum(is.na(x))) # count each columns 

table(LUAD_feature$donor_age_at_diagnosis)
LUAD_feature[LUAD_feature$donor_age_at_diagnosis > 0 & LUAD_feature$donor_age_at_diagnosis < 10, "donor_age_at_diagnosis"] = 1
LUAD_feature[LUAD_feature$donor_age_at_diagnosis > 9 & LUAD_feature$donor_age_at_diagnosis < 20, "donor_age_at_diagnosis"] = 10
LUAD_feature[LUAD_feature$donor_age_at_diagnosis > 19 & LUAD_feature$donor_age_at_diagnosis < 30, "donor_age_at_diagnosis"] = 20
LUAD_feature[LUAD_feature$donor_age_at_diagnosis > 29 & LUAD_feature$donor_age_at_diagnosis < 40, "donor_age_at_diagnosis"] = 30
LUAD_feature[LUAD_feature$donor_age_at_diagnosis > 39 & LUAD_feature$donor_age_at_diagnosis < 50, "donor_age_at_diagnosis"] = 40
LUAD_feature[LUAD_feature$donor_age_at_diagnosis > 49 & LUAD_feature$donor_age_at_diagnosis < 60, "donor_age_at_diagnosis"] = 50
LUAD_feature[LUAD_feature$donor_age_at_diagnosis > 59 & LUAD_feature$donor_age_at_diagnosis < 70, "donor_age_at_diagnosis"] = 60
LUAD_feature[LUAD_feature$donor_age_at_diagnosis > 69 & LUAD_feature$donor_age_at_diagnosis < 80, "donor_age_at_diagnosis"] = 70
LUAD_feature[LUAD_feature$donor_age_at_diagnosis > 79 & LUAD_feature$donor_age_at_diagnosis < 90, "donor_age_at_diagnosis"] = 80
LUAD_feature[LUAD_feature$donor_age_at_diagnosis > 89 & LUAD_feature$donor_age_at_diagnosis < 100, "donor_age_at_diagnosis"] = 90
table(LUAD_feature$donor_age_at_diagnosis)



a1<-as.data.table(table(LUAD_feature$donor_sex))
a2<-as.data.table(table(LUAD_feature$donor_age_at_diagnosis))
a3<-as.data.table(table(LUAD_feature$race))
a4<-as.data.table(table(LUAD_feature$Stage))
LUAD_feature_table<-rbind(a1, a2, a3, a4)
LUAD_feature_table

#THCA
THCA<-read.csv("/home/solbi/20190118/cancersample/tcga_primary_dataset_3782_geneexpres_dataset_stage_featuredata_THCA-US_python.csv")

THCA_feature<-THCA[,c(1:14)]
head(THCA_feature)
sapply(THCA_feature, function(x) sum(is.na(x))) # count each columns 

table(THCA_feature$donor_age_at_diagnosis)
THCA_feature[THCA_feature$donor_age_at_diagnosis > 0 & THCA_feature$donor_age_at_diagnosis < 10, "donor_age_at_diagnosis"] = 1
THCA_feature[THCA_feature$donor_age_at_diagnosis > 9 & THCA_feature$donor_age_at_diagnosis < 20, "donor_age_at_diagnosis"] = 10
THCA_feature[THCA_feature$donor_age_at_diagnosis > 19 & THCA_feature$donor_age_at_diagnosis < 30, "donor_age_at_diagnosis"] = 20
THCA_feature[THCA_feature$donor_age_at_diagnosis > 29 & THCA_feature$donor_age_at_diagnosis < 40, "donor_age_at_diagnosis"] = 30
THCA_feature[THCA_feature$donor_age_at_diagnosis > 39 & THCA_feature$donor_age_at_diagnosis < 50, "donor_age_at_diagnosis"] = 40
THCA_feature[THCA_feature$donor_age_at_diagnosis > 49 & THCA_feature$donor_age_at_diagnosis < 60, "donor_age_at_diagnosis"] = 50
THCA_feature[THCA_feature$donor_age_at_diagnosis > 59 & THCA_feature$donor_age_at_diagnosis < 70, "donor_age_at_diagnosis"] = 60
THCA_feature[THCA_feature$donor_age_at_diagnosis > 69 & THCA_feature$donor_age_at_diagnosis < 80, "donor_age_at_diagnosis"] = 70
THCA_feature[THCA_feature$donor_age_at_diagnosis > 79 & THCA_feature$donor_age_at_diagnosis < 90, "donor_age_at_diagnosis"] = 80
THCA_feature[THCA_feature$donor_age_at_diagnosis > 89 & THCA_feature$donor_age_at_diagnosis < 100, "donor_age_at_diagnosis"] = 90
table(THCA_feature$donor_age_at_diagnosis)



a1<-as.data.table(table(THCA_feature$donor_sex))
a2<-as.data.table(table(THCA_feature$donor_age_at_diagnosis))
a3<-as.data.table(table(THCA_feature$race))
a4<-as.data.table(table(THCA_feature$Stage))
THCA_feature_table<-rbind(a1, a2, a3, a4)
THCA_feature_table

colnames(BRCA_feature_table) = c("feature", "BRCA")
colnames(KIRC_feature_table) = c("feature", "KIRC")
colnames(UCEC_feature_table) = c("feature", "UCEC")
colnames(LUAD_feature_table) = c("feature", "LUAD")
colnames(THCA_feature_table) = c("feature", "THCA")

cancer_table1<-merge(x = BRCA_feature_table, y = KIRC_feature_table, by = 'feature', all = T)
cancer_table2<-merge(x = cancer_table1, y = UCEC_feature_table, by = 'feature', all = T)
cancer_table3<-merge(x = cancer_table2, y = LUAD_feature_table, by = 'feature', all = T)
cancer_table4<-merge(x = cancer_table3, y = THCA_feature_table, by = 'feature', all = T)

write.csv(cancer_table4, file = "/home/solbi/20190118/cancersample/table_cancerTOP5.csv", row.names = FALSE)



