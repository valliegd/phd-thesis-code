library(dplyr)
library(plyr)
library(tidyr)
library(stringr)
library(data.table)
library(tidyverse)

setwd("/home/vgalassideforie/re_gecip/neurology/Valentina/Complete_quad_analysis")
raw_quad <- read.table(file = 'rare_diseases_family_2022-04-26_10-22-52.tsv', sep = '\t', header = TRUE)
total_analysis <- read.table(file = 'rare_disease_analysis_2022-04-26_10-48-28.tsv', sep = '\t', header = TRUE)
complete_rawquad_table <- left_join(raw_quad, total_analysis, by = "Rare.Diseases.Family.Id")
# selecting families with proband, mother, father and siblings of any combination
complete_rawquad_table <- complete_rawquad_table[complete.cases(complete_rawquad_table), ]


amended_complete_rawquad_table <- complete_rawquad_table
amended_complete_rawquad_table <- within(amended_complete_rawquad_table, Biological.Relationship.To.Proband[Biological.Relationship.To.Proband == "Twins Unknown" & Rare.Diseases.Family.Id == "a"] <- "Twins Monozygous")
amended_complete_rawquad_table <- within(amended_complete_rawquad_table, Biological.Relationship.To.Proband[Biological.Relationship.To.Proband == "Twins Unknown" & Rare.Diseases.Family.Id == "b"] <- "Twins Monozygous")
amended_complete_rawquad_table <- within(amended_complete_rawquad_table, Biological.Relationship.To.Proband[Biological.Relationship.To.Proband == "Twins Unknown" & Rare.Diseases.Family.Id == "c"] <- "Twins Monozygous")
amended_complete_rawquad_table <- within(amended_complete_rawquad_table, Biological.Relationship.To.Proband[Biological.Relationship.To.Proband == "Twins Unknown" & Rare.Diseases.Family.Id == "d"] <- "Twins Dizygous")
amended_complete_rawquad_table <- within(amended_complete_rawquad_table, Biological.Relationship.To.Proband[Biological.Relationship.To.Proband == "Twins Unknown" & Rare.Diseases.Family.Id == "e"] <- "Twins Dizygous")
amended_complete_rawquad_table <- within(amended_complete_rawquad_table, Biological.Relationship.To.Proband[Biological.Relationship.To.Proband == "Twins Unknown" & Rare.Diseases.Family.Id == "f"] <- "Twins Dizygous")
amended_complete_rawquad_table <- within(amended_complete_rawquad_table, Biological.Relationship.To.Proband[Biological.Relationship.To.Proband == "Twins Unknown" & Rare.Diseases.Family.Id == "g"] <- "Twins Dizygous")
amended_complete_rawquad_table <- within(amended_complete_rawquad_table, Biological.Relationship.To.Proband[Biological.Relationship.To.Proband == "Twins Unknown" & Rare.Diseases.Family.Id == "h"] <- "Twins Dizygous")

amended_complete_rawquad_table_proband <- amended_complete_rawquad_table %>% filter(Biological.Relationship.To.Proband == "N/A")
amended_complete_rawquad_table_unique <- amended_complete_rawquad_table_proband %>% select(Rare.Diseases.Family.Id, Family.Group.Type)
amended_complete_rawquad_table_unique <- amended_complete_rawquad_table_unique %>% unique()
amended_complete_rawquad_table_proband <- amended_complete_rawquad_table_proband %>% select(Rare.Diseases.Family.Id, Participant.Id, Plate.Key, Genome.Build, Genetic.Vs.Reported.Results)
amended_complete_rawquad_table_proband$ProbandID = amended_complete_rawquad_table_proband$Participant.Id
amended_complete_rawquad_table_proband$ProbandKey = amended_complete_rawquad_table_proband$Plate.Key
amended_complete_rawquad_table_proband$ProbandGenome = amended_complete_rawquad_table_proband$Genome.Build
amended_complete_rawquad_table_proband$ProbandGeneticsvsReported = amended_complete_rawquad_table_proband$Genetic.Vs.Reported.Results

# filtering out mother and father from original table and joining to table with unique family_ids
amended_complete_rawquad_table_mother <- amended_complete_rawquad_table %>% filter(Biological.Relationship.To.Proband == "Mother")
amended_complete_rawquad_table_mother <- amended_complete_rawquad_table_mother %>% select(Rare.Diseases.Family.Id, Participant.Id, Plate.Key, Genome.Build, Genetic.Vs.Reported.Results)
amended_complete_rawquad_table_mother$MothersID = amended_complete_rawquad_table_mother$Participant.Id
amended_complete_rawquad_table_mother$MothersKey = amended_complete_rawquad_table_mother$Plate.Key
amended_complete_rawquad_table_mother$MothersGenome = amended_complete_rawquad_table_mother$Genome.Build
amended_complete_rawquad_table_mother$MothersGeneticsvsReported = amended_complete_rawquad_table_mother$Genetic.Vs.Reported.Results
amended_complete_rawquad_table_father <- amended_complete_rawquad_table %>% filter(Biological.Relationship.To.Proband == "Father")
amended_complete_rawquad_table_father <- amended_complete_rawquad_table_father %>% select(Rare.Diseases.Family.Id, Participant.Id, Plate.Key, Genome.Build, Genetic.Vs.Reported.Results)
amended_complete_rawquad_table_father$FathersID = amended_complete_rawquad_table_father$Participant.Id
amended_complete_rawquad_table_father$FathersKey = amended_complete_rawquad_table_father$Plate.Key
amended_complete_rawquad_table_father$FathersGenome = amended_complete_rawquad_table_father$Genome.Build
amended_complete_rawquad_table_father$FathersGeneticsvsReported = amended_complete_rawquad_table_father$Genetic.Vs.Reported.Results
amended_complete_rawquad_table_fullsib <- amended_complete_rawquad_table %>% filter(Biological.Relationship.To.Proband == "Full Sibling")
amended_complete_rawquad_table_fullsib <- amended_complete_rawquad_table_fullsib %>% select(Rare.Diseases.Family.Id, Participant.Id, Plate.Key, Genome.Build, Genetic.Vs.Reported.Results)
amended_complete_rawquad_table_fullsib$FullsiblingID = amended_complete_rawquad_table_fullsib$Participant.Id
amended_complete_rawquad_table_fullsib$FullsiblingKey = amended_complete_rawquad_table_fullsib$Plate.Key
amended_complete_rawquad_table_fullsib$FullsiblingGenome = amended_complete_rawquad_table_fullsib$Genome.Build
amended_complete_rawquad_table_fullsib$FullsiblingGeneticsvsReported = amended_complete_rawquad_table_fullsib$Genetic.Vs.Reported.Results

n_occur <- amended_complete_rawquad_table_fullsib %>% group_by(Rare.Diseases.Family.Id) %>% filter(n()>1) %>% dplyr::summarize(n=n())

amended_complete_rawquad_table_twinsdz <- amended_complete_rawquad_table %>% filter(Biological.Relationship.To.Proband == "Twins Dizygous")
amended_complete_rawquad_table_twinsdz <- amended_complete_rawquad_table_twinsdz %>% select(Rare.Diseases.Family.Id, Participant.Id, Plate.Key, Genome.Build, Genetic.Vs.Reported.Results)
amended_complete_rawquad_table_twinsdz$TwinDZID = amended_complete_rawquad_table_twinsdz$Participant.Id
amended_complete_rawquad_table_twinsdz$TwinDZKey = amended_complete_rawquad_table_twinsdz$Plate.Key
amended_complete_rawquad_table_twinsdz$TwinDZGenome = amended_complete_rawquad_table_twinsdz$Genome.Build
amended_complete_rawquad_table_twinsdz$TwinDZGeneticsvsReported = amended_complete_rawquad_table_twinsdz$Genetic.Vs.Reported.Results

amended_complete_rawquad_table_mz <- amended_complete_rawquad_table %>% filter(Biological.Relationship.To.Proband == "Twins Monozygous")
amended_complete_rawquad_table_mz <- amended_complete_rawquad_table_mz %>% select(Rare.Diseases.Family.Id, Participant.Id, Plate.Key, Genome.Build, Genetic.Vs.Reported.Results)
amended_complete_rawquad_table_mz$TwinMZID = amended_complete_rawquad_table_mz$Participant.Id
amended_complete_rawquad_table_mz$TwinMZKey = amended_complete_rawquad_table_mz$Plate.Key
amended_complete_rawquad_table_mz$TwinMZGenome = amended_complete_rawquad_table_mz$Genome.Build
amended_complete_rawquad_table_mz$TwinMZGeneticsvsReported = amended_complete_rawquad_table_mz$Genetic.Vs.Reported.Results

Joinmothers = left_join(amended_complete_rawquad_table_unique, amended_complete_rawquad_table_mother, by = "Rare.Diseases.Family.Id")
Joinparents = left_join(Joinmothers, amended_complete_rawquad_table_father, by = "Rare.Diseases.Family.Id")
Joinallproband = left_join(Joinparents, amended_complete_rawquad_table_proband, by = "Rare.Diseases.Family.Id")
Joinall = left_join(Joinallproband, amended_complete_rawquad_table_fullsib, by = "Rare.Diseases.Family.Id")
Joinall <- Joinall %>% unique() %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, FullsiblingID, FullsiblingKey, FullsiblingGenome, FullsiblingGeneticsvsReported)
Joinall <- Joinall[complete.cases(Joinall), ]
Joinall_genome38 <- Joinall %>% filter(ProbandGenome == "GRCh38", MothersGenome == "GRCh38", FathersGenome == "GRCh38", FullsiblingGenome == "GRCh38")
#572 families with full siblings as well that have been sequenced against GRCh38
Joinall_genome38_geneticsvsreported <- Joinall_genome38 %>% filter(ProbandGeneticsvsReported == "familyPassesGvsRChecks", MothersGeneticsvsReported == "familyPassesGvsRChecks", FathersGeneticsvsReported == "familyPassesGvsRChecks", FullsiblingGeneticsvsReported == "familyPassesGvsRChecks")
#414 families with full siblings as well that have been sequenced against GRCh38 and passed genetics vs reported

Joinalltwindz = left_join(Joinallproband, amended_complete_rawquad_table_twinsdz, by = "Rare.Diseases.Family.Id")
Joinalltwindz <- Joinalltwindz %>% unique() %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, TwinDZID, TwinDZKey, TwinDZGenome, TwinDZGeneticsvsReported)
Joinalltwindz <- Joinalltwindz[complete.cases(Joinalltwindz), ]
Joinall_genome38_twinsdz <- Joinalltwindz %>% filter(ProbandGenome == "GRCh38", MothersGenome == "GRCh38", FathersGenome == "GRCh38", TwinDZGenome == "GRCh38")
#12 families with dizygous twins as well that have been sequenced against GRCh38
Joinall_genome38_geneticsvsreported_twindz <- Joinall_genome38_twinsdz %>% filter(ProbandGeneticsvsReported == "familyPassesGvsRChecks", MothersGeneticsvsReported == "familyPassesGvsRChecks", FathersGeneticsvsReported == "familyPassesGvsRChecks", TwinDZGeneticsvsReported == "familyPassesGvsRChecks")
#11 families with dizygous twins as well that have been sequenced against GRCh38 and passed genetics vs reported - INCLUDING ORIGINAL UNKNOWN

Joinalltwinmz = left_join(Joinallproband, amended_complete_rawquad_table_mz, by = "Rare.Diseases.Family.Id")
Joinalltwinmz <- Joinalltwinmz %>% unique() %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, TwinMZID, TwinMZKey, TwinMZGenome, TwinMZGeneticsvsReported)
Joinalltwinmz <- Joinalltwinmz[complete.cases(Joinalltwinmz), ]
Joinall_genome38_twinsmz <- Joinalltwinmz %>% filter(ProbandGenome == "GRCh38", MothersGenome == "GRCh38", FathersGenome == "GRCh38", TwinMZGenome == "GRCh38")
#7 families with monozygous twins as well that have been sequenced against GRCh38
Joinall_genome38_geneticsvsreported_twinmz <- Joinall_genome38_twinsmz %>% filter(ProbandGeneticsvsReported == "familyPassesGvsRChecks", MothersGeneticsvsReported == "familyPassesGvsRChecks", FathersGeneticsvsReported == "familyPassesGvsRChecks", TwinMZGeneticsvsReported == "familyPassesGvsRChecks")
#5 families with monozygous twins as well that have been sequenced against GRCh38 and passed genetics vs reported

Joinalltwinunknown = left_join(Joinallproband, complete_rawquad_table_twinunknown, by = "Rare.Diseases.Family.Id")
Joinalltwinunknown <- Joinalltwinunknown %>% unique() %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, TwinUnknownID, TwinUnknownKey, TwinUnknownGenome, TwinUnknownGeneticsvsReported)
Joinalltwinunknown <- Joinalltwinunknown[complete.cases(Joinalltwinunknown), ]
Joinall_genome38_twinsunknown <- Joinalltwinunknown %>% filter(ProbandGenome == "GRCh38", MothersGenome == "GRCh38", FathersGenome == "GRCh38", TwinUnknownGenome == "GRCh38")
#8 families with unknown twins as well that have been sequenced against GRCh38
Joinall_genome38_geneticsvsreported_twinunknown <- Joinall_genome38_twinsunknown %>% filter(ProbandGeneticsvsReported == "familyPassesGvsRChecks", MothersGeneticsvsReported == "familyPassesGvsRChecks", FathersGeneticsvsReported == "familyPassesGvsRChecks", TwinUnknownGeneticsvsReported == "familyPassesGvsRChecks")
#4 families with unknown twins as well that have been sequenced against GRCh38 and passed genetics vs reported

#GRCh37
Joinall_genome37 <- Joinall %>% filter(ProbandGenome == "GRCh37", MothersGenome == "GRCh37", FathersGenome == "GRCh37", FullsiblingGenome == "GRCh37")
#108 families with full siblings as well that have been sequenced against GRCh37
Joinall_genome37_twinsdz <- Joinalltwindz %>% filter(ProbandGenome == "GRCh37", MothersGenome == "GRCh37", FathersGenome == "GRCh37", TwinDZGenome == "GRCh37")
#5 families with dizygous twins as well that have been sequenced against GRCh37 - INCLUDING ORIGINALLY UNKNOWN
Joinall_genome37_twinsmz <- Joinalltwinmz %>% filter(ProbandGenome == "GRCh37", MothersGenome == "GRCh37", FathersGenome == "GRCh37", TwinMZGenome == "GRCh37")
#3 families with monozygous twins as well that have been sequenced against GRCh37 - INCLUDING ORIGINALLY UNKNOWN
#Joinall_genome37_twinsunknown <- Joinalltwinunknown %>% filter(ProbandGenome == "GRCh37", MothersGenome == "GRCh37", FathersGenome == "GRCh37", TwinUnknownGenome == "GRCh37")
#4 families with unknown twins as well that have been sequenced against GRCh37


interpreted <- read.delim("rare_disease_interpreted_2021-06-08_16-03-49.tsv")
quad_GRCh37 <- amended_complete_rawquad_table %>% filter(Genome.Build == "GRCh37" & Family.Group.Type == "Families with more than three Participants") %>% unique()
quad_GRCh37 <- quad_GRCh37 %>% unique()
quad_GRCh37_unique <- quad_GRCh37$Rare.Diseases.Family.Id %>% unique() %>% length() 
# 170 genomes
quad_GRCh37_interpreted <- interpreted %>% filter(Assembly == "GRCh37") %>% select(Plate.Key) %>% unique() %>% pull()
# 7754 genomes
quad_GRCh37_interpreted_total <- quad_GRCh37 %>% filter(Plate.Key %in% quad_GRCh37_interpreted)
# 356 genomes
quad_GRCh37_interpreted_total_unique <- quad_GRCh37_interpreted_total$Rare.Diseases.Family.Id %>% unique() %>% length() 
# 94 families
Proband37 <- quad_GRCh37_interpreted_total %>% filter(Biological.Relationship.To.Proband == "N/A")
Proband_37 <- Proband37 %>% select(Rare.Diseases.Family.Id, Participant.Id)
Proband_37$ProbandID = Proband_37$Participant.Id

# filtering out mother and father from original table and joining to table with unique family_ids
Mother37 <- quad_GRCh37_interpreted_total %>% filter(Biological.Relationship.To.Proband == "Mother")
Mother_37 <- Mother37 %>% select(Rare.Diseases.Family.Id, Participant.Id)
Mother_37$MothersID = Mother_37$Participant.Id
Father37 <- quad_GRCh37_interpreted_total %>% filter(Biological.Relationship.To.Proband == "Father")
Father_37 <- Father37 %>% select(Rare.Diseases.Family.Id, Participant.Id)
Father_37$FathersID = Father_37$Participant.Id

Fullsibling37 <- quad_GRCh37_interpreted_total %>% filter(Biological.Relationship.To.Proband == "Full Sibling")
Fullsibling_37 <- Fullsibling37 %>% select(Rare.Diseases.Family.Id, Participant.Id)
Fullsibling_37$FullsiblingID = Fullsibling_37$Participant.Id

TwinsDZ37 <- quad_GRCh37_interpreted_total %>% filter(Biological.Relationship.To.Proband == "Twins Dizygous")
TwinsDZ_37 <- TwinsDZ37 %>% select(Rare.Diseases.Family.Id, Participant.Id)
TwinsDZ_37$TwinDZID = TwinsDZ_37$Participant.Id

TwinsMZ37 <- quad_GRCh37_interpreted_total %>% filter(Biological.Relationship.To.Proband == "Twins Monozygous")
TwinsMZ_37 <- TwinsMZ37 %>% select(Rare.Diseases.Family.Id, Participant.Id)
TwinsMZ_37$TwinMZID = TwinsMZ_37$Participant.Id

Joinmothers37 = left_join(quad_GRCh37_interpreted_total, Mother_37, by = "Rare.Diseases.Family.Id")
Joinparents37 = left_join(Joinmothers37, Father_37, by = "Rare.Diseases.Family.Id")
Joinall37 = left_join(Joinparents37, Proband_37, by = "Rare.Diseases.Family.Id")
Joinallfullsib37 = left_join(Joinall37, Fullsibling_37, by = "Rare.Diseases.Family.Id")
JoinallDZ37 = left_join(Joinall37, TwinsDZ_37, by = "Rare.Diseases.Family.Id")
JoinallMZ37 = left_join(Joinall37, TwinsMZ_37, by = "Rare.Diseases.Family.Id")

# tidy up
Jointidyfullsib37 <- Joinallfullsib37 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, MothersID, FathersID, FullsiblingID)
Jointidyfullsib37 <- Jointidyfullsib37 %>% unique()
Jointidyfullsib37 <- Jointidyfullsib37[complete.cases(Jointidyfullsib37), ]
# 59 full sib families with IDs without NAs
JointidyDZ37 <- JoinallDZ37 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, MothersID, FathersID, TwinDZID)
JointidyDZ37 <- JointidyDZ37 %>% unique()
JointidyDZ37 <- JointidyDZ37[complete.cases(JointidyDZ37), ]
# 4 DZ families with IDs without NAs
JointidyMZ37 <- JoinallMZ37 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, MothersID, FathersID, TwinMZID)
JointidyMZ37 <- JointidyMZ37 %>% unique()
JointidyMZ37 <- JointidyMZ37[complete.cases(JointidyMZ37), ]
# 2 MZ families with IDs without NAs

ProbandKey37 <- Proband37 %>% select(Rare.Diseases.Family.Id, Plate.Key)
ProbandKey37$ProbandKey = ProbandKey37$Plate.Key
MotherKey37 <- Mother37 %>% select(Rare.Diseases.Family.Id, Plate.Key)
MotherKey37$MothersKey = MotherKey37$Plate.Key
FatherKey37 <- Father37 %>% select(Rare.Diseases.Family.Id, Plate.Key)
FatherKey37$FathersKey = FatherKey37$Plate.Key
FullsiblingKey37 <- Fullsibling37 %>% select(Rare.Diseases.Family.Id, Plate.Key)
FullsiblingKey37$FullsiblingKey = FullsiblingKey37$Plate.Key
TwinsDZKey37 <- TwinsDZ37 %>% select(Rare.Diseases.Family.Id, Plate.Key)
TwinsDZKey37$TwinDZKey = TwinsDZKey37$Plate.Key
TwinsMZKey37 <- TwinsMZ37 %>% select(Rare.Diseases.Family.Id, Plate.Key)
TwinsMZKey37$TwinMZKey = TwinsMZKey37$Plate.Key

Joinmotherkey37 = left_join(Jointidyfullsib37, MotherKey37, by = "Rare.Diseases.Family.Id")
Joinparentkey37 = left_join(Joinmotherkey37, FatherKey37, by = "Rare.Diseases.Family.Id")
Joinallkey37 = left_join(Joinparentkey37, ProbandKey37, by = "Rare.Diseases.Family.Id")
Joinallfullsibkey37 = left_join(Joinallkey37, FullsiblingKey37, by = "Rare.Diseases.Family.Id")
Jointidyfullsibkey37 <- Joinallfullsibkey37 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, MothersID, MothersKey, FathersID, FathersKey, FullsiblingID, FullsiblingKey)
Jointidyfullsibkey37 <- Jointidyfullsibkey37[!duplicated(Jointidyfullsibkey37$Rare.Diseases.Family.Id), ]
Jointidyfullsibkey37 <- Jointidyfullsibkey37[complete.cases(Jointidyfullsibkey37), ]
# 50 full sib families
JoinallkeyDZ37 = left_join(JointidyDZ37, MotherKey37, by = "Rare.Diseases.Family.Id")
JoinallkeyDZ37 = left_join(JoinallkeyDZ37, FatherKey37, by = "Rare.Diseases.Family.Id")
JoinallkeyDZ37 = left_join(JoinallkeyDZ37, ProbandKey37, by = "Rare.Diseases.Family.Id")
JoinallkeyDZ37 = left_join(JoinallkeyDZ37, TwinsDZKey37, by = "Rare.Diseases.Family.Id")
JoinallkeyDZ37 <- JoinallkeyDZ37 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, MothersID, MothersKey, FathersID, FathersKey, TwinDZID, TwinDZKey)
JoinallkeyDZ37 <- JoinallkeyDZ37[!duplicated(JoinallkeyDZ37$Rare.Diseases.Family.Id), ]
JoinallkeyDZ37 <- JoinallkeyDZ37[complete.cases(JoinallkeyDZ37), ]
# 4 DZ families 
JoinallkeyMZ37 = left_join(JointidyMZ37, MotherKey37, by = "Rare.Diseases.Family.Id")
JoinallkeyMZ37 = left_join(JoinallkeyMZ37, FatherKey37, by = "Rare.Diseases.Family.Id")
JoinallkeyMZ37 = left_join(JoinallkeyMZ37, ProbandKey37, by = "Rare.Diseases.Family.Id")
JoinallkeyMZ37 = left_join(JoinallkeyMZ37, TwinsMZKey37, by = "Rare.Diseases.Family.Id")
JoinallkeyMZ37 <- JoinallkeyMZ37 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, MothersID, MothersKey, FathersID, FathersKey, TwinMZID, TwinMZKey)
JoinallkeyMZ37 <- JoinallkeyMZ37[!duplicated(JoinallkeyMZ37$Rare.Diseases.Family.Id), ]
JoinallkeyMZ37 <- JoinallkeyMZ37[complete.cases(JoinallkeyMZ37), ]
# 2 MZ families 

#Add genome build and genetics vs reported
ProbandBuild37 <- Proband37 %>% select(Rare.Diseases.Family.Id, Genome.Build)
ProbandBuild37$ProbandGenome = ProbandBuild37$Genome.Build
MotherBuild37 <- Mother37 %>% select(Rare.Diseases.Family.Id, Genome.Build)
MotherBuild37$MothersGenome = MotherBuild37$Genome.Build
FatherBuild37 <- Father37 %>% select(Rare.Diseases.Family.Id, Genome.Build)
FatherBuild37$FathersGenome = FatherBuild37$Genome.Build
FullsiblingBuild37 <- Fullsibling37 %>% select(Rare.Diseases.Family.Id, Genome.Build)
FullsiblingBuild37$FullsiblingGenome = FullsiblingBuild37$Genome.Build
TwinsDZBuild37 <- TwinsDZ37 %>% select(Rare.Diseases.Family.Id, Genome.Build)
TwinsDZBuild37$TwinDZGenome = TwinsDZBuild37$Genome.Build
TwinsMZBuild37 <- TwinsMZ37 %>% select(Rare.Diseases.Family.Id, Genome.Build)
TwinsMZBuild37$TwinMZGenome = TwinsMZBuild37$Genome.Build

Joinallfullsibbuild37 = left_join(Jointidyfullsibkey37, MotherBuild37, by = "Rare.Diseases.Family.Id")
Joinallfullsibbuild37 = left_join(Joinallfullsibbuild37, FatherBuild37, by = "Rare.Diseases.Family.Id")
Joinallfullsibbuild37 = left_join(Joinallfullsibbuild37, ProbandBuild37, by = "Rare.Diseases.Family.Id")
Joinallfullsibbuild37 = left_join(Joinallfullsibbuild37, FullsiblingBuild37, by = "Rare.Diseases.Family.Id")
Joinallfullsibbuild37 <- Joinallfullsibbuild37 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, MothersID, MothersKey, MothersGenome, FathersID, FathersKey, FathersGenome, FullsiblingID, FullsiblingKey, FullsiblingGenome)
Joinallfullsibbuild37 <- Joinallfullsibbuild37[!duplicated(Joinallfullsibbuild37$Rare.Diseases.Family.Id), ]
Joinallfullsibbuild37 <- Joinallfullsibbuild37[complete.cases(Joinallfullsibbuild37), ]
# 50 full sib families
JoinallbuildDZ37 = left_join(JoinallkeyDZ37, MotherBuild37, by = "Rare.Diseases.Family.Id")
JoinallbuildDZ37 = left_join(JoinallbuildDZ37, FatherBuild37, by = "Rare.Diseases.Family.Id")
JoinallbuildDZ37 = left_join(JoinallbuildDZ37, ProbandBuild37, by = "Rare.Diseases.Family.Id")
JoinallbuildDZ37 = left_join(JoinallbuildDZ37, TwinsDZBuild37, by = "Rare.Diseases.Family.Id")
JoinallbuildDZ37 <- JoinallbuildDZ37 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, MothersID, MothersKey, MothersGenome, FathersID, FathersKey, FathersGenome, TwinDZID, TwinDZKey, TwinDZGenome)
JoinallbuildDZ37 <- JoinallbuildDZ37[!duplicated(JoinallbuildDZ37$Rare.Diseases.Family.Id), ]
JoinallbuildDZ37 <- JoinallbuildDZ37[complete.cases(JoinallbuildDZ37), ]
# 4 DZ families 
JoinallbuildMZ37 = left_join(JoinallkeyMZ37, MotherBuild37, by = "Rare.Diseases.Family.Id")
JoinallbuildMZ37 = left_join(JoinallbuildMZ37, FatherBuild37, by = "Rare.Diseases.Family.Id")
JoinallbuildMZ37 = left_join(JoinallbuildMZ37, ProbandBuild37, by = "Rare.Diseases.Family.Id")
JoinallbuildMZ37 = left_join(JoinallbuildMZ37, TwinsMZBuild37, by = "Rare.Diseases.Family.Id")
JoinallbuildMZ37 <- JoinallbuildMZ37 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, MothersID, MothersKey, MothersGenome, FathersID, FathersKey, FathersGenome, TwinMZID, TwinMZKey, TwinMZGenome)
JoinallbuildMZ37 <- JoinallbuildMZ37[!duplicated(JoinallbuildMZ37$Rare.Diseases.Family.Id), ]
JoinallbuildMZ37 <- JoinallbuildMZ37[complete.cases(JoinallbuildMZ37), ]
# 2 MZ families 

ProbandCheck37 <- Proband37 %>% select(Rare.Diseases.Family.Id, Genetic.Vs.Reported.Results)
ProbandCheck37$ProbandGeneticsvsReported = ProbandCheck37$Genetic.Vs.Reported.Results
MotherCheck37 <- Mother37 %>% select(Rare.Diseases.Family.Id, Genetic.Vs.Reported.Results)
MotherCheck37$MothersGeneticsvsReported = MotherCheck37$Genetic.Vs.Reported.Results
FatherCheck37 <- Father37 %>% select(Rare.Diseases.Family.Id, Genetic.Vs.Reported.Results)
FatherCheck37$FathersGeneticsvsReported = FatherCheck37$Genetic.Vs.Reported.Results
FullsiblingCheck37 <- Fullsibling37 %>% select(Rare.Diseases.Family.Id, Genetic.Vs.Reported.Results)
FullsiblingCheck37$FullsiblingGeneticsvsReported = FullsiblingCheck37$Genetic.Vs.Reported.Results
TwinsDZCheck37 <- TwinsDZ37 %>% select(Rare.Diseases.Family.Id, Genetic.Vs.Reported.Results)
TwinsDZCheck37$TwinDZGeneticsvsReported = TwinsDZCheck37$Genetic.Vs.Reported.Results
TwinsMZCheck37 <- TwinsMZ37 %>% select(Rare.Diseases.Family.Id, Genetic.Vs.Reported.Results)
TwinsMZCheck37$TwinMZGeneticsvsReported = TwinsMZCheck37$Genetic.Vs.Reported.Results

Joinallfullsibcheck37 = left_join(Joinallfullsibbuild37, MotherCheck37, by = "Rare.Diseases.Family.Id")
Joinallfullsibcheck37 = left_join(Joinallfullsibcheck37, FatherCheck37, by = "Rare.Diseases.Family.Id")
Joinallfullsibcheck37 = left_join(Joinallfullsibcheck37, ProbandCheck37, by = "Rare.Diseases.Family.Id")
Joinallfullsibcheck37 = left_join(Joinallfullsibcheck37, FullsiblingCheck37, by = "Rare.Diseases.Family.Id")
Joinallfullsibcheck37 <- Joinallfullsibcheck37 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, FullsiblingID, FullsiblingKey, FullsiblingGenome, FullsiblingGeneticsvsReported)
Joinallfullsibcheck37 <- Joinallfullsibcheck37[!duplicated(Joinallfullsibcheck37$Rare.Diseases.Family.Id), ]
Joinallfullsibcheck37 <- Joinallfullsibcheck37[complete.cases(Joinallfullsibcheck37), ]
# 50 full sib families
JoinallcheckDZ37 = left_join(JoinallbuildDZ37, MotherCheck37, by = "Rare.Diseases.Family.Id")
JoinallcheckDZ37 = left_join(JoinallcheckDZ37, FatherCheck37, by = "Rare.Diseases.Family.Id")
JoinallcheckDZ37 = left_join(JoinallcheckDZ37, ProbandCheck37, by = "Rare.Diseases.Family.Id")
JoinallcheckDZ37 = left_join(JoinallcheckDZ37, TwinsDZCheck37, by = "Rare.Diseases.Family.Id")
JoinallcheckDZ37 <- JoinallcheckDZ37 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, TwinDZID, TwinDZKey, TwinDZGenome, TwinDZGeneticsvsReported)
JoinallcheckDZ37 <- JoinallcheckDZ37[!duplicated(JoinallcheckDZ37$Rare.Diseases.Family.Id), ]
JoinallcheckDZ37 <- JoinallcheckDZ37[complete.cases(JoinallcheckDZ37), ]
# 4 DZ families 
JoinallcheckMZ37 = left_join(JoinallbuildMZ37, MotherCheck37, by = "Rare.Diseases.Family.Id")
JoinallcheckMZ37 = left_join(JoinallcheckMZ37, FatherCheck37, by = "Rare.Diseases.Family.Id")
JoinallcheckMZ37 = left_join(JoinallcheckMZ37, ProbandCheck37, by = "Rare.Diseases.Family.Id")
JoinallcheckMZ37 = left_join(JoinallcheckMZ37, TwinsMZCheck37, by = "Rare.Diseases.Family.Id")
JoinallcheckMZ37 <- JoinallcheckMZ37 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, TwinMZID, TwinMZKey, TwinMZGenome, TwinMZGeneticsvsReported)
JoinallcheckMZ37 <- JoinallcheckMZ37[!duplicated(JoinallcheckMZ37$Rare.Diseases.Family.Id), ]
JoinallcheckMZ37 <- JoinallcheckMZ37[complete.cases(JoinallcheckMZ37), ]
# 2 MZ families 

Quad_fullsibling_metastart <- rbind(Joinall_genome38_geneticsvsreported, Joinallfullsibcheck37)
Quad_twinsdizygous_metastart <- rbind(Joinall_genome38_geneticsvsreported_twindz, JoinallcheckDZ37)
Quad_twinsmonozygous_metastart <- rbind(Joinall_genome38_geneticsvsreported_twinmz, JoinallcheckMZ37)

write.table(Quad_twinsdizygous_metastart, file = "/home/vgalassideforie/re_gecip/neurology/Valentina/Complete_quad_analysis/Quad_twinsdizygous_metastart.tsv", row.names=FALSE, sep="\t", quote = FALSE)
write.table(Quad_twinsmonozygous_metastart, file = "/home/vgalassideforie/re_gecip/neurology/Valentina/Complete_quad_analysis/Quad_twinsmonozygous_metastart.tsv", row.names=FALSE, sep="\t", quote = FALSE)

n_occur <- Quad_fullsibling_metastart %>% group_by(Rare.Diseases.Family.Id) %>% filter(n()>1) %>% dplyr::summarize(n=n())
n_occur <- Joinall_genome38_geneticsvsreported %>% group_by(Rare.Diseases.Family.Id) %>% filter(n()>1) %>% dplyr::summarize(n=n())


# amending for multiple siblings
n_occur_fullsib <- as.data.frame(Quad_fullsibling_metastart[which(Quad_fullsibling_metastart$Rare.Diseases.Family.Id %in% n_occur$Rare.Diseases.Family.Id),])
n_occur_fullsib <- n_occur_fullsib %>% select(Rare.Diseases.Family.Id, FullsiblingID)
n_occur_fullsib <- n_occur_fullsib %>% group_by(Rare.Diseases.Family.Id)

#add EH3 expansion sizes
raw_eh3repeat <- read.delim("table_STR_repeat_size_each_row_allele_EHv3.2.2_HTT_simplified.tsv")
raw_eh3repeat_grouped <- raw_eh3repeat %>% group_by(platekey) %>% dplyr::summarise(repeat_size_list = toString(repeat_size)) %>% as.data.frame()
amended_complete_rawquad_table$platekey = amended_complete_rawquad_table$Plate.Key
plate = left_join(raw_eh3repeat_grouped, amended_complete_rawquad_table, by = "platekey")
plate <- plate %>% select(Rare.Diseases.Family.Id, Family.Group.Type, Biological.Relationship.To.Proband, platekey, repeat_size_list, Participant.Id, Genome.Build, Genetic.Vs.Reported.Results)
plate <- plate %>% filter(Family.Group.Type == "Families with more than three Participants")
plate <- plate[!duplicated(plate$Rare.Diseases.Family.Id), ] 
plate <- plate[complete.cases(plate), ]

plate_proband <- plate %>% filter(Biological.Relationship.To.Proband == "N/A")
plate_proband <- plate_proband %>% select(Rare.Diseases.Family.Id, repeat_size_list)
plate_proband$Probandrepeatsize = plate_proband$repeat_size_list

plate_mother <- plate %>% filter(Biological.Relationship.To.Proband == "Mother")
plate_mother <- plate_mother %>% select(Rare.Diseases.Family.Id, repeat_size_list)
plate_mother$Mothersrepeatsize = plate_mother$repeat_size_list

plate_father <- plate %>% filter(Biological.Relationship.To.Proband == "Father")
plate_father <- plate_father %>% select(Rare.Diseases.Family.Id, repeat_size_list)
plate_father$Fathersrepeatsize = plate_father$repeat_size_list

plate_fullsib <- plate %>% filter(Biological.Relationship.To.Proband == "Full Sibling")
plate_fullsib <- plate_fullsib %>% select(Rare.Diseases.Family.Id, repeat_size_list)
plate_fullsib$Fullsiblingrepeatsize = plate_fullsib$repeat_size_list

plate_twinsdz <- plate %>% filter(Biological.Relationship.To.Proband == "Twins Dizygous")
plate_twinsdz <- plate_twinsdz %>% select(Rare.Diseases.Family.Id, repeat_size_list)
plate_twinsdz$TwinsDZrepeatsize = plate_twinsdz$repeat_size_list

plate_twinsmz <- plate %>% filter(Biological.Relationship.To.Proband == "Twins Monozygous")
plate_twinsmz <- plate_twinsmz %>% select(Rare.Diseases.Family.Id, repeat_size_list)
plate_twinsmz$TwinsMZrepeatsize = plate_twinsmz$repeat_size_list

plate_twinsunknown <- plate %>% filter(Biological.Relationship.To.Proband == "Twins Unknown")
plate_twinsunknown <- plate_twinsunknown %>% select(Rare.Diseases.Family.Id, repeat_size_list)
plate_twinsunknown$TwinsUnknownrepeatsize = plate_twinsunknown$repeat_size_list

Joinmothersrepeat = left_join(Joinall_genome38_geneticsvsreported, plate_mother, by = "Rare.Diseases.Family.Id")
Joinfathersrepeat = left_join(Joinmothersrepeat, plate_father, by = "Rare.Diseases.Family.Id")
Joinprobandrepeat = left_join(Joinfathersrepeat, plate_proband, by = "Rare.Diseases.Family.Id")
Joinfullrepeat = left_join(Joinprobandrepeat, plate_fullsib, by = "Rare.Diseases.Family.Id")
Joinfullrepeat <- Joinfullrepeat %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, Probandrepeatsize, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, Mothersrepeatsize, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, Fathersrepeatsize, FullsiblingID, FullsiblingKey, FullsiblingGenome, FullsiblingGeneticsvsReported, Fullsiblingrepeatsize)

Joinmothersrepeat = left_join(Joinall_genome38_geneticsvsreported_twindz, plate_mother, by = "Rare.Diseases.Family.Id")
Joinfathersrepeat = left_join(Joinmothersrepeat, plate_father, by = "Rare.Diseases.Family.Id")
Joinprobandrepeat = left_join(Joinfathersrepeat, plate_proband, by = "Rare.Diseases.Family.Id")
Jointwindzrepeat = left_join(Joinprobandrepeat, plate_twinsdz, by = "Rare.Diseases.Family.Id")
Jointwindzrepeat <- Jointwindzrepeat %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, Probandrepeatsize, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, Mothersrepeatsize, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, Fathersrepeatsize, TwinDZID, TwinDZKey, TwinDZGenome, TwinDZGeneticsvsReported, TwinsDZrepeatsize)

Joinmothersrepeat = left_join(Joinall_genome38_geneticsvsreported_twinmz, plate_mother, by = "Rare.Diseases.Family.Id")
Joinfathersrepeat = left_join(Joinmothersrepeat, plate_father, by = "Rare.Diseases.Family.Id")
Joinprobandrepeat = left_join(Joinfathersrepeat, plate_proband, by = "Rare.Diseases.Family.Id")
Jointwinmzrepeat = left_join(Joinprobandrepeat, plate_twinsmz, by = "Rare.Diseases.Family.Id")
Jointwinmzrepeat <- Jointwinmzrepeat %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, Probandrepeatsize, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, Mothersrepeatsize, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, Fathersrepeatsize, TwinMZID, TwinMZKey, TwinMZGenome, TwinMZGeneticsvsReported, TwinsMZrepeatsize)

Joinmothersrepeat = left_join(Joinall_genome38_geneticsvsreported_twinunknown, plate_mother, by = "Rare.Diseases.Family.Id")
Joinfathersrepeat = left_join(Joinmothersrepeat, plate_father, by = "Rare.Diseases.Family.Id")
Joinprobandrepeat = left_join(Joinfathersrepeat, plate_proband, by = "Rare.Diseases.Family.Id")
Jointwinunknownrepeat = left_join(Joinprobandrepeat, plate_twinsunknown, by = "Rare.Diseases.Family.Id")
Jointwinunknownrepeat <- Jointwinunknownrepeat %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, Probandrepeatsize, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, Mothersrepeatsize, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, Fathersrepeatsize, TwinUnknownID, TwinUnknownKey, TwinUnknownGenome, TwinUnknownGeneticsvsReported, TwinsUnknownrepeatsize)
Jointwinunknownrepeat <- Jointwinunknownrepeat[!duplicated(Jointwinunknownrepeat$Rare.Diseases.Family.Id), ]

Joinprobandrepeat_findna <- Joinfullrepeat %>% select(ProbandKey, Probandrepeatsize) %>% filter(is.na(Probandrepeatsize)) %>% select(ProbandKey)
names(Joinprobandrepeat_findna)[names(Joinprobandrepeat_findna) == 'ProbandKey'] <- 'Platekey'
Joinmotherrepeat_findna <- Joinfullrepeat %>% select(MothersKey, Mothersrepeatsize) %>% filter(is.na(Mothersrepeatsize)) %>% select(MothersKey)
names(Joinmotherrepeat_findna)[names(Joinmotherrepeat_findna) == 'MothersKey'] <- 'Platekey'
Joinfatherrepeat_findna <- Joinfullrepeat %>% select(FathersKey, Fathersrepeatsize) %>% filter(is.na(Fathersrepeatsize)) %>% select(FathersKey)
names(Joinfatherrepeat_findna)[names(Joinfatherrepeat_findna) == 'FathersKey'] <- 'Platekey'
Joinfullrepeat_findna <- Joinfullrepeat %>% select(FullsiblingKey, Fullsiblingrepeatsize) %>% filter(is.na(Fullsiblingrepeatsize)) %>% select(FullsiblingKey)
names(Joinfullrepeat_findna)[names(Joinfullrepeat_findna) == 'FullsiblingKey'] <- 'Platekey'
Jointwindzrepeat_findna <- Jointwindzrepeat %>% select(TwinDZKey, TwinsDZrepeatsize) %>% filter(is.na(TwinsDZrepeatsize)) %>% select(TwinDZKey)
names(Jointwindzrepeat_findna)[names(Jointwindzrepeat_findna) == 'TwinDZKey'] <- 'Platekey'
Jointwinmzrepeat_findna <- Jointwinmzrepeat %>% select(TwinMZKey, TwinsMZrepeatsize) %>% filter(is.na(TwinsMZrepeatsize)) %>% select(TwinMZKey)
names(Jointwinmzrepeat_findna)[names(Jointwinmzrepeat_findna) == 'TwinMZKey'] <- 'Platekey'
Jointwinunknownrepeat_findna <- Jointwinunknownrepeat %>% select(TwinUnknownKey, TwinsUnknownrepeatsize) %>% filter(is.na(TwinsUnknownrepeatsize)) %>% select(TwinUnknownKey)
names(Jointwinunknownrepeat_findna)[names(Jointwinunknownrepeat_findna) == 'TwinUnknownKey'] <- 'Platekey'
quadrepeatsnotfound <- do.call("rbind", list(Joinprobandrepeat_findna, Joinmotherrepeat_findna, Joinfatherrepeat_findna, Joinfullrepeat_findna, Jointwindzrepeat_findna, Jointwinmzrepeat_findna, Jointwinunknownrepeat_findna))
write.table(quadrepeatsnotfound, file = "/home/vgalassideforie/re_gecip/shared_allGeCIPs/AD_VGD/quadrepeatsnotfound.tsv", row.names=FALSE, sep="\t")
quadrepeatsnotfound <- read.delim("quadrepeatsnotfound.tsv")

platekeyfilepaths <- read.delim("genome_file_paths_and_types_2021-10-13_10-34-19.tsv")
platekeyfilepaths <- platekeyfilepaths %>% filter(grepl('.bam', File.Path))
platekeyfilepaths <- platekeyfilepaths %>% select(Platekey, File.Path)
platekeyfilepaths_common <- as.data.frame(platekeyfilepaths[which(platekeyfilepaths$Platekey %in% quadrepeatsnotfound$Platekey),])
platekeyfilepaths_common <- platekeyfilepaths_common %>% select(Platekey)
platekeyfilepaths_common$Platekey = paste(platekeyfilepaths_common$Platekey,"_HTT.vcf", sep = "")
platekeyfilepaths_common <- platekeyfilepaths_common %>% unique()
write.table(platekeyfilepaths_common, file = "/home/vgalassideforie/re_gecip/neurology/Valentina/Running_EHv3/EH3_quad_output/platekey_common.txt", row.names=FALSE, sep="\t", quote = FALSE)

# extract repeat size values from EH3
quad_vcfextraction <- read.csv("Quad_vcfextraction.csv")
quad_vcfextraction <- mutate(quad_vcfextraction, Genotype = gsub(":.*", "", quad_vcfextraction$X.1))
quad_vcfextraction <- mutate(quad_vcfextraction, Type = gsub(":\\d+.*", "", quad_vcfextraction$X.1))
quad_vcfextraction <- mutate(quad_vcfextraction, Type = gsub(".*:", "", quad_vcfextraction$Type))

quad_vcfextraction <- mutate(quad_vcfextraction, A = gsub("[0-9]/[0-9]:[A-Z]*/[A-Z]*:", "", quad_vcfextraction$X.1))
quad_vcfextraction <- mutate(quad_vcfextraction, Allele_1 = gsub("/.*", "", quad_vcfextraction$A))
quad_vcfextraction <- mutate(quad_vcfextraction, A = sub("\\d*?/", "", quad_vcfextraction$A))
quad_vcfextraction <- mutate(quad_vcfextraction, Allele_2 = gsub(":.*", "", quad_vcfextraction$A))
quad_vcfextraction <- mutate(quad_vcfextraction, A = sub("\\d*?:", "", quad_vcfextraction$A))
quad_vcfextraction <- mutate(quad_vcfextraction, Range_Allele_1 = gsub("/.*", "", quad_vcfextraction$A))
quad_vcfextraction <- mutate(quad_vcfextraction, A = sub("\\d*?-\\d*?/", "", quad_vcfextraction$A))
quad_vcfextraction <- mutate(quad_vcfextraction, Range_Allele_2 = gsub(":.*", "", quad_vcfextraction$A))
quad_vcfextraction <- mutate(quad_vcfextraction, A = sub("\\d*?-\\d*?:", "", quad_vcfextraction$A))
quad_vcfextraction <- mutate(quad_vcfextraction, Spanning = gsub(":.*", "", quad_vcfextraction$A))
quad_vcfextraction <- mutate(quad_vcfextraction, A = sub("\\d*?/\\d*?:", "", quad_vcfextraction$A))
quad_vcfextraction <- mutate(quad_vcfextraction, Flanking = gsub(":.*", "", quad_vcfextraction$A))
quad_vcfextraction <- mutate(quad_vcfextraction, A = sub("\\d*?/\\d*?:", "", quad_vcfextraction$A))
quad_vcfextraction <- mutate(quad_vcfextraction, Inrepeat = gsub(":.*", "", quad_vcfextraction$A))
quad_vcfextraction <- mutate(quad_vcfextraction, Locus_coverage = gsub(".*:", "", quad_vcfextraction$A))
names(quad_vcfextraction)[1] <- "platekey"
quad_vcfextraction <- quad_vcfextraction %>% select(platekey, Genotype, Type, Allele_1, Allele_2, Range_Allele_1, Range_Allele_2, Spanning, Flanking, Inrepeat, Locus_coverage)
quad_vcfextraction <- mutate(quad_vcfextraction, Type_allele_1 = gsub("/.*", "", quad_vcfextraction$Type))
quad_vcfextraction <- mutate(quad_vcfextraction, Type_allele_2 = gsub(".*/", "", quad_vcfextraction$Type))
quad_vcfextraction <- mutate(quad_vcfextraction, Spanning_allele_1 = gsub("/.*", "", quad_vcfextraction$Spanning))
quad_vcfextraction <- mutate(quad_vcfextraction, Spanning_allele_2 = gsub(".*/", "", quad_vcfextraction$Spanning))
quad_vcfextraction <- mutate(quad_vcfextraction, Flanking_allele_1 = gsub("/.*", "", quad_vcfextraction$Flanking))
quad_vcfextraction <- mutate(quad_vcfextraction, Flanking_allele_2 = gsub(".*/", "", quad_vcfextraction$Flanking))
quad_vcfextraction <- mutate(quad_vcfextraction, Inrepeat_allele_1 = gsub("/.*", "", quad_vcfextraction$Inrepeat))
quad_vcfextraction <- mutate(quad_vcfextraction, Inrepeat_allele_2 = gsub(".*/", "", quad_vcfextraction$Inrepeat))
quad_vcfextraction <- quad_vcfextraction %>% select(platekey, Allele_1, Allele_2, Spanning_allele_1, Spanning_allele_2, Flanking_allele_1, Flanking_allele_2, Inrepeat_allele_1, Inrepeat_allele_2)
Alleleperrow <- quad_vcfextraction %>% select(platekey, Allele_1, Allele_2)
Alleleperrow$repeat_size_list <- paste(Alleleperrow$Allele_1, Alleleperrow$Allele_2, sep = ",")
Alleleperrow <- Alleleperrow %>% select(platekey, repeat_size_list)
full_repeatsizes <- rbind(raw_eh3repeat_grouped, Alleleperrow, by = "platekey")

amended_complete_rawquad_table$platekey = amended_complete_rawquad_table$Plate.Key
plate = left_join(full_repeatsizes, amended_complete_rawquad_table, by = "platekey")

plate_proband <- plate %>% filter(Biological.Relationship.To.Proband == "N/A")
plate_proband <- plate_proband %>% select(Rare.Diseases.Family.Id, repeat_size_list)
plate_proband$Probandrepeatsize = plate_proband$repeat_size_list

plate_mother <- plate %>% filter(Biological.Relationship.To.Proband == "Mother")
plate_mother <- plate_mother %>% select(Rare.Diseases.Family.Id, repeat_size_list)
plate_mother$Mothersrepeatsize = plate_mother$repeat_size_list

plate_father <- plate %>% filter(Biological.Relationship.To.Proband == "Father")
plate_father <- plate_father %>% select(Rare.Diseases.Family.Id, repeat_size_list)
plate_father$Fathersrepeatsize = plate_father$repeat_size_list

plate_fullsib <- plate %>% filter(Biological.Relationship.To.Proband == "Full Sibling")
plate_fullsib <- plate_fullsib %>% select(Rare.Diseases.Family.Id, repeat_size_list)
plate_fullsib$Fullsiblingrepeatsize = plate_fullsib$repeat_size_list

plate_twinsdz <- plate %>% filter(Biological.Relationship.To.Proband == "Twins Dizygous")
plate_twinsdz <- plate_twinsdz %>% select(Rare.Diseases.Family.Id, repeat_size_list)
plate_twinsdz$TwinsDZrepeatsize = plate_twinsdz$repeat_size_list

plate_twinsmz <- plate %>% filter(Biological.Relationship.To.Proband == "Twins Monozygous")
plate_twinsmz <- plate_twinsmz %>% select(Rare.Diseases.Family.Id, repeat_size_list)
plate_twinsmz$TwinsMZrepeatsize = plate_twinsmz$repeat_size_list

plate_twinsunknown <- plate %>% filter(Biological.Relationship.To.Proband == "Twins Unknown")
plate_twinsunknown <- plate_twinsunknown %>% select(Rare.Diseases.Family.Id, repeat_size_list)
plate_twinsunknown$TwinsUnknownrepeatsize = plate_twinsunknown$repeat_size_list

Joinmothersrepeat = left_join(Joinall_genome38_geneticsvsreported, plate_mother, by = "Rare.Diseases.Family.Id")
Joinfathersrepeat = left_join(Joinmothersrepeat, plate_father, by = "Rare.Diseases.Family.Id")
Joinprobandrepeat = left_join(Joinfathersrepeat, plate_proband, by = "Rare.Diseases.Family.Id")
Joinfullrepeat = left_join(Joinprobandrepeat, plate_fullsib, by = "Rare.Diseases.Family.Id")
Joinfullrepeat <- Joinfullrepeat %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, Probandrepeatsize, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, Mothersrepeatsize, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, Fathersrepeatsize, FullsiblingID, FullsiblingKey, FullsiblingGenome, FullsiblingGeneticsvsReported, Fullsiblingrepeatsize)

Joinmothersrepeat = left_join(Joinall_genome38_geneticsvsreported_twindz, plate_mother, by = "Rare.Diseases.Family.Id")
Joinfathersrepeat = left_join(Joinmothersrepeat, plate_father, by = "Rare.Diseases.Family.Id")
Joinprobandrepeat = left_join(Joinfathersrepeat, plate_proband, by = "Rare.Diseases.Family.Id")
Jointwindzrepeat = left_join(Joinprobandrepeat, plate_twinsdz, by = "Rare.Diseases.Family.Id")
Jointwindzrepeat <- Jointwindzrepeat %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, Probandrepeatsize, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, Mothersrepeatsize, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, Fathersrepeatsize, TwinDZID, TwinDZKey, TwinDZGenome, TwinDZGeneticsvsReported, TwinsDZrepeatsize)

Joinmothersrepeat = left_join(Joinall_genome38_geneticsvsreported_twinmz, plate_mother, by = "Rare.Diseases.Family.Id")
Joinfathersrepeat = left_join(Joinmothersrepeat, plate_father, by = "Rare.Diseases.Family.Id")
Joinprobandrepeat = left_join(Joinfathersrepeat, plate_proband, by = "Rare.Diseases.Family.Id")
Jointwinmzrepeat = left_join(Joinprobandrepeat, plate_twinsmz, by = "Rare.Diseases.Family.Id")
Jointwinmzrepeat <- Jointwinmzrepeat %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, Probandrepeatsize, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, Mothersrepeatsize, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, Fathersrepeatsize, TwinMZID, TwinMZKey, TwinMZGenome, TwinMZGeneticsvsReported, TwinsMZrepeatsize)

Joinmothersrepeat = left_join(Joinall_genome38_geneticsvsreported_twinunknown, plate_mother, by = "Rare.Diseases.Family.Id")
Joinfathersrepeat = left_join(Joinmothersrepeat, plate_father, by = "Rare.Diseases.Family.Id")
Joinprobandrepeat = left_join(Joinfathersrepeat, plate_proband, by = "Rare.Diseases.Family.Id")
Jointwinunknownrepeat = left_join(Joinprobandrepeat, plate_twinsunknown, by = "Rare.Diseases.Family.Id")
Jointwinunknownrepeat <- Jointwinunknownrepeat %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, Probandrepeatsize, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, Mothersrepeatsize, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, Fathersrepeatsize, TwinUnknownID, TwinUnknownKey, TwinUnknownGenome, TwinUnknownGeneticsvsReported, TwinsUnknownrepeatsize)

Joinmothersrepeat = left_join(Quad_fullsibling_metastart, plate_mother, by = "Rare.Diseases.Family.Id")
Joinfathersrepeat = left_join(Joinmothersrepeat, plate_father, by = "Rare.Diseases.Family.Id")
Joinprobandrepeat = left_join(Joinfathersrepeat, plate_proband, by = "Rare.Diseases.Family.Id")
Joinfullrepeat = left_join(Joinprobandrepeat, plate_fullsib, by = "Rare.Diseases.Family.Id")
Joinfullrepeat <- Joinfullrepeat %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, Probandrepeatsize, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, Mothersrepeatsize, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, Fathersrepeatsize, FullsiblingID, FullsiblingKey, FullsiblingGenome, FullsiblingGeneticsvsReported, Fullsiblingrepeatsize)
Joinfullrepeat <- Joinfullrepeat[!duplicated(Joinfullrepeat$FullsiblingID), ]

# replace entries with 4 repeat sizes with correct 2 - due to sequenced against both 37 and 38
Joinfullrepeat <- within(Joinfullrepeat, Probandrepeatsize[Probandrepeatsize == "17, 17, 16, 16" & Rare.Diseases.Family.Id == "i"] <- "17, 16")
Joinfullrepeat <- within(Joinfullrepeat, Probandrepeatsize[Probandrepeatsize == "17, 17, 18, 18" & Rare.Diseases.Family.Id == "j"] <- "17, 18")
Joinfullrepeat <- within(Joinfullrepeat, Mothersrepeatsize[Mothersrepeatsize == "17, 17, 13, 13" & Rare.Diseases.Family.Id == "k"] <- "17, 13")
Joinfullrepeat <- within(Joinfullrepeat, Mothersrepeatsize[Mothersrepeatsize == "17, 17, 17, 17" & Rare.Diseases.Family.Id == "l"] <- "17, 17")
Joinfullrepeat <- within(Joinfullrepeat, Mothersrepeatsize[Mothersrepeatsize == "15, 15, 15, 15" & Rare.Diseases.Family.Id == "m"] <- "15, 15")
Joinfullrepeat <- within(Joinfullrepeat, Mothersrepeatsize[Mothersrepeatsize == "18, 18, 22, 22" & Rare.Diseases.Family.Id == "n"] <- "18, 22")
Joinfullrepeat <- within(Joinfullrepeat, Mothersrepeatsize[Mothersrepeatsize == "17, 17, 20, 20" & Rare.Diseases.Family.Id == "o"] <- "17, 20")
Joinfullrepeat <- within(Joinfullrepeat, Mothersrepeatsize[Mothersrepeatsize == "18, 18, 19, 19" & Rare.Diseases.Family.Id == "p"] <- "18, 19")
Joinfullrepeat <- within(Joinfullrepeat, Fathersrepeatsize[Fathersrepeatsize == "24, 24, 25, 25" & Rare.Diseases.Family.Id == "q"] <- "24, 25")
Joinfullrepeat <- within(Joinfullrepeat, Fathersrepeatsize[Fathersrepeatsize == "17, 17, 20, 20" & Rare.Diseases.Family.Id == "r"] <- "17, 20")
Joinfullrepeat <- within(Joinfullrepeat, Fathersrepeatsize[Fathersrepeatsize == "19, 19, 15, 15" & Rare.Diseases.Family.Id == "s"] <- "19, 15")
Joinfullrepeat <- within(Joinfullrepeat, Fathersrepeatsize[Fathersrepeatsize == "19, 19, 30, 30" & Rare.Diseases.Family.Id == "t"] <- "19, 30")
Joinfullrepeat <- within(Joinfullrepeat, Fullsiblingrepeatsize[Fullsiblingrepeatsize == "17, 17, 13, 13" & Rare.Diseases.Family.Id == "u"] <- "17, 13")
Joinfullrepeat <- within(Joinfullrepeat, Fullsiblingrepeatsize[Fullsiblingrepeatsize == "18, 18, 16, 16" & Rare.Diseases.Family.Id == "v"] <- "18, 16")
Joinfullrepeat <- within(Joinfullrepeat, Fullsiblingrepeatsize[Fullsiblingrepeatsize == "15, 15, 15, 15" & Rare.Diseases.Family.Id == "w"] <- "15, 15")

#split so each allele in each column
probandsplitrepeats <- Joinfullrepeat
probandsplitrepeats <- probandsplitrepeats %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, Probandrepeatsize)
splitrepeatsproband <- probandsplitrepeats %>% separate(Probandrepeatsize, c("Proband_allele_1", "Proband_allele_2"), ",")
mothersplitrepeats <- Joinfullrepeat
mothersplitrepeats <- mothersplitrepeats %>% select(Rare.Diseases.Family.Id, Family.Group.Type, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, Mothersrepeatsize)
splitrepeatsmother <- mothersplitrepeats %>% separate(Mothersrepeatsize, c("Mothers_allele_1", "Mothers_allele_2"), ",")
fathersplitrepeats <- Joinfullrepeat
fathersplitrepeats <- fathersplitrepeats %>% select(Rare.Diseases.Family.Id, Family.Group.Type, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, Fathersrepeatsize)
splitrepeatsfather <- fathersplitrepeats %>% separate(Fathersrepeatsize, c("Fathers_allele_1", "Fathers_allele_2"), ",")
fullsibsplitrepeats <- Joinfullrepeat
fullsibsplitrepeats <- fullsibsplitrepeats %>% select(Rare.Diseases.Family.Id, Family.Group.Type, FullsiblingID, FullsiblingKey, FullsiblingGenome, FullsiblingGeneticsvsReported, Fullsiblingrepeatsize)
splitrepeatsfullsib <- fullsibsplitrepeats %>% separate(Fullsiblingrepeatsize, c("Fullsibling_allele_1", "Fullsibling_allele_2"), ",")
Joinrepeatsizeparents <- left_join(splitrepeatsmother, splitrepeatsfather, by = "Rare.Diseases.Family.Id")
Joinrepeatsizetrio <- left_join(Joinrepeatsizeparents, splitrepeatsproband, by = "Rare.Diseases.Family.Id")
Joinrepeatsizeall <- left_join(Joinrepeatsizetrio, splitrepeatsfullsib, by = "Rare.Diseases.Family.Id")
Joinrepeatsizetidy <- Joinrepeatsizeall[!duplicated(Joinrepeatsizeall$FullsiblingID), ]
Joinrepeatsizetidy <- Joinrepeatsizetidy %>% select(Rare.Diseases.Family.Id, Family.Group.Type.x, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, Proband_allele_1, Proband_allele_2, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, Mothers_allele_1, Mothers_allele_2, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, Fathers_allele_1, Fathers_allele_2, FullsiblingID, FullsiblingKey, FullsiblingGenome, FullsiblingGeneticsvsReported, Fullsibling_allele_1, Fullsibling_allele_2)

Joinmothersrepeat = left_join(Quad_twinsdizygous_metastart, plate_mother, by = "Rare.Diseases.Family.Id")
Joinfathersrepeat = left_join(Joinmothersrepeat, plate_father, by = "Rare.Diseases.Family.Id")
Joinprobandrepeat = left_join(Joinfathersrepeat, plate_proband, by = "Rare.Diseases.Family.Id")
Jointwindzrepeat = left_join(Joinprobandrepeat, plate_twinsdz, by = "Rare.Diseases.Family.Id")
Jointwindzrepeat <- Jointwindzrepeat %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, Probandrepeatsize, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, Mothersrepeatsize, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, Fathersrepeatsize, TwinDZID, TwinDZKey, TwinDZGenome, TwinDZGeneticsvsReported, TwinsDZrepeatsize)

#split so each allele in each column
probandsplitrepeats <- Jointwindzrepeat
probandsplitrepeats <- probandsplitrepeats %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, Probandrepeatsize)
splitrepeatsproband <- probandsplitrepeats %>% separate(Probandrepeatsize, c("Proband_allele_1", "Proband_allele_2"), ",")
mothersplitrepeats <- Jointwindzrepeat
mothersplitrepeats <- mothersplitrepeats %>% select(Rare.Diseases.Family.Id, Family.Group.Type, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, Mothersrepeatsize)
splitrepeatsmother <- mothersplitrepeats %>% separate(Mothersrepeatsize, c("Mothers_allele_1", "Mothers_allele_2"), ",")
fathersplitrepeats <- Jointwindzrepeat
fathersplitrepeats <- fathersplitrepeats %>% select(Rare.Diseases.Family.Id, Family.Group.Type, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, Fathersrepeatsize)
splitrepeatsfather <- fathersplitrepeats %>% separate(Fathersrepeatsize, c("Fathers_allele_1", "Fathers_allele_2"), ",")
twindzsplitrepeats <- Jointwindzrepeat
twindzsplitrepeats <- twindzsplitrepeats %>% select(Rare.Diseases.Family.Id, Family.Group.Type, TwinDZID, TwinDZKey, TwinDZGenome, TwinDZGeneticsvsReported, TwinsDZrepeatsize)
splitrepeatstwindz <- twindzsplitrepeats %>% separate(TwinsDZrepeatsize, c("TwinsDZ_allele_1", "TwinsDZ_allele_2"), ",")
Joinrepeatsizeparents <- left_join(splitrepeatsmother, splitrepeatsfather, by = "Rare.Diseases.Family.Id")
Joinrepeatsizetrio <- left_join(Joinrepeatsizeparents, splitrepeatsproband, by = "Rare.Diseases.Family.Id")
Joinrepeatsizeall <- left_join(Joinrepeatsizetrio, splitrepeatstwindz, by = "Rare.Diseases.Family.Id")
Joinrepeatsizetidy_twindz <- Joinrepeatsizeall[!duplicated(Joinrepeatsizeall$Rare.Diseases.Family.Id), ]
Joinrepeatsizetidy_twindz <- Joinrepeatsizetidy_twindz %>% select(Rare.Diseases.Family.Id, Family.Group.Type.x, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, Proband_allele_1, Proband_allele_2, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, Mothers_allele_1, Mothers_allele_2, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, Fathers_allele_1, Fathers_allele_2, TwinDZID, TwinDZKey, TwinDZGenome, TwinDZGeneticsvsReported, TwinsDZ_allele_1, TwinsDZ_allele_2)

Fullsib_families_and_IDs <- Quad_fullsibling_metastart %>% select(Rare.Diseases.Family.Id, ProbandID, MothersID, FathersID, FullsiblingID)

Joinmothersrepeat = left_join(Quad_twinsmonozygous_metastart, plate_mother, by = "Rare.Diseases.Family.Id")
Joinfathersrepeat = left_join(Joinmothersrepeat, plate_father, by = "Rare.Diseases.Family.Id")
Joinprobandrepeat = left_join(Joinfathersrepeat, plate_proband, by = "Rare.Diseases.Family.Id")
Jointwinmzrepeat = left_join(Joinprobandrepeat, plate_twinsmz, by = "Rare.Diseases.Family.Id")
Jointwinmzrepeat <- Jointwinmzrepeat %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, Probandrepeatsize, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, Mothersrepeatsize, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, Fathersrepeatsize, TwinMZID, TwinMZKey, TwinMZGenome, TwinMZGeneticsvsReported, TwinsMZrepeatsize)

#split so each allele in each column
probandsplitrepeats <- Jointwinmzrepeat
probandsplitrepeats <- probandsplitrepeats %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, Probandrepeatsize)
splitrepeatsproband <- probandsplitrepeats %>% separate(Probandrepeatsize, c("Proband_allele_1", "Proband_allele_2"), ",")
mothersplitrepeats <- Jointwinmzrepeat
mothersplitrepeats <- mothersplitrepeats %>% select(Rare.Diseases.Family.Id, Family.Group.Type, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, Mothersrepeatsize)
splitrepeatsmother <- mothersplitrepeats %>% separate(Mothersrepeatsize, c("Mothers_allele_1", "Mothers_allele_2"), ",")
fathersplitrepeats <- Jointwinmzrepeat
fathersplitrepeats <- fathersplitrepeats %>% select(Rare.Diseases.Family.Id, Family.Group.Type, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, Fathersrepeatsize)
splitrepeatsfather <- fathersplitrepeats %>% separate(Fathersrepeatsize, c("Fathers_allele_1", "Fathers_allele_2"), ",")
twinmzsplitrepeats <- Jointwinmzrepeat
twinmzsplitrepeats <- twinmzsplitrepeats %>% select(Rare.Diseases.Family.Id, Family.Group.Type, TwinMZID, TwinMZKey, TwinMZGenome, TwinMZGeneticsvsReported, TwinsMZrepeatsize)
splitrepeatstwinmz <- twinmzsplitrepeats %>% separate(TwinsMZrepeatsize, c("TwinsMZ_allele_1", "TwinsMZ_allele_2"), ",")
Joinrepeatsizeparents <- left_join(splitrepeatsmother, splitrepeatsfather, by = "Rare.Diseases.Family.Id")
Joinrepeatsizetrio <- left_join(Joinrepeatsizeparents, splitrepeatsproband, by = "Rare.Diseases.Family.Id")
Joinrepeatsizeall <- left_join(Joinrepeatsizetrio, splitrepeatstwinmz, by = "Rare.Diseases.Family.Id")
Joinrepeatsizetidy_twinmz <- Joinrepeatsizeall[!duplicated(Joinrepeatsizeall$Rare.Diseases.Family.Id), ]
Joinrepeatsizetidy_twinmz <- Joinrepeatsizetidy_twinmz %>% select(Rare.Diseases.Family.Id, Family.Group.Type.x, ProbandID, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, Proband_allele_1, Proband_allele_2, MothersID, MothersKey, MothersGenome, MothersGeneticsvsReported, Mothers_allele_1, Mothers_allele_2, FathersID, FathersKey, FathersGenome, FathersGeneticsvsReported, Fathers_allele_1, Fathers_allele_2, TwinMZID, TwinMZKey, TwinMZGenome, TwinMZGeneticsvsReported, TwinsMZ_allele_1, TwinsMZ_allele_2)

# merging 3 repeat size tables into one with extra column stating the sibling type
Joinrepeatsizetidy$Sibling <- "Full sibling"
Joinrepeatsizetidy_twindz$Sibling <- "Twins Dizygous"
Joinrepeatsizetidy_twinmz$Sibling <- "Twins Monozygous"
names(Joinrepeatsizetidy)[names(Joinrepeatsizetidy) == 'FullsiblingID'] <- 'SiblingID'
names(Joinrepeatsizetidy)[names(Joinrepeatsizetidy) == 'FullsiblingKey'] <- 'SiblingKey'
names(Joinrepeatsizetidy)[names(Joinrepeatsizetidy) == 'FullsiblingGenome'] <- 'SiblingGenome'
names(Joinrepeatsizetidy)[names(Joinrepeatsizetidy) == 'FullsiblingGeneticsvsReported'] <- 'SiblingGeneticsvsReported'
names(Joinrepeatsizetidy)[names(Joinrepeatsizetidy) == 'Fullsibling_allele_1'] <- 'Sibling_allele_1'
names(Joinrepeatsizetidy)[names(Joinrepeatsizetidy) == 'Fullsibling_allele_2'] <- 'Sibling_allele_2'
names(Joinrepeatsizetidy_twindz)[names(Joinrepeatsizetidy_twindz) == 'TwinsDZ_allele_1'] <- 'Sibling_allele_1'
names(Joinrepeatsizetidy_twindz)[names(Joinrepeatsizetidy_twindz) == 'TwinsDZ_allele_2'] <- 'Sibling_allele_2'
names(Joinrepeatsizetidy_twindz)[names(Joinrepeatsizetidy_twindz) == 'TwinDZID'] <- 'SiblingID'
names(Joinrepeatsizetidy_twindz)[names(Joinrepeatsizetidy_twindz) == 'TwinDZKey'] <- 'SiblingKey'
names(Joinrepeatsizetidy_twindz)[names(Joinrepeatsizetidy_twindz) == 'TwinDZGenome'] <- 'SiblingGenome'
names(Joinrepeatsizetidy_twindz)[names(Joinrepeatsizetidy_twindz) == 'TwinDZGeneticsvsReported'] <- 'SiblingGeneticsvsReported'
names(Joinrepeatsizetidy_twinmz)[names(Joinrepeatsizetidy_twinmz) == 'TwinsMZ_allele_1'] <- 'Sibling_allele_1'
names(Joinrepeatsizetidy_twinmz)[names(Joinrepeatsizetidy_twinmz) == 'TwinsMZ_allele_2'] <- 'Sibling_allele_2'
names(Joinrepeatsizetidy_twinmz)[names(Joinrepeatsizetidy_twinmz) == 'TwinMZID'] <- 'SiblingID'
names(Joinrepeatsizetidy_twinmz)[names(Joinrepeatsizetidy_twinmz) == 'TwinMZKey'] <- 'SiblingKey'
names(Joinrepeatsizetidy_twinmz)[names(Joinrepeatsizetidy_twinmz) == 'TwinMZGenome'] <- 'SiblingGenome'
names(Joinrepeatsizetidy_twinmz)[names(Joinrepeatsizetidy_twinmz) == 'TwinMZGeneticsvsReported'] <- 'SiblingGeneticsvsReported'

repeatsizelist_allsiblingtypes <- rbind(Joinrepeatsizetidy, Joinrepeatsizetidy_twindz)
repeatsizelist_allsiblingtypes <- rbind(repeatsizelist_allsiblingtypes, Joinrepeatsizetidy_twinmz)

write.table(repeatsizelist_allsiblingtypes, file = "/home/vgalassideforie/re_gecip/neurology/Valentina/Complete_quad_analysis/completequadlist_allsiblingtypes.tsv", row.names=FALSE, sep="\t", quote = FALSE)

#adding metadata

# One row per participant - full siblings family
Fullsib_families_and_IDs <- Quad_fullsibling_metastart %>% select(Rare.Diseases.Family.Id, ProbandID, MothersID, FathersID, FullsiblingID)
Fullsib_families_and_IDs <- data.frame(Fullsib_families_and_IDs[1], stack(Fullsib_families_and_IDs[2:ncol(Fullsib_families_and_IDs)]))
Fullsib_families_and_IDs$Participant.Id <- Fullsib_families_and_IDs$values
Fullsib_families_and_IDs <- Fullsib_families_and_IDs %>% select(Rare.Diseases.Family.Id, Participant.Id)
# One row per participant - twins dizygous family
Twindz_families_and_IDs <- Quad_twinsdizygous_metastart %>% select(Rare.Diseases.Family.Id, ProbandID, MothersID, FathersID, TwinDZID)
Twindz_families_and_IDs <- data.frame(Twindz_families_and_IDs[1], stack(Twindz_families_and_IDs[2:ncol(Twindz_families_and_IDs)]))
Twindz_families_and_IDs$Participant.Id <- Twindz_families_and_IDs$values
Twindz_families_and_IDs <- Twindz_families_and_IDs %>% select(Rare.Diseases.Family.Id, Participant.Id)
# One row per participant - twins monozygous family
Twinmz_families_and_IDs <- Quad_twinsmonozygous_metastart %>% select(Rare.Diseases.Family.Id, ProbandID, MothersID, FathersID, TwinMZID)
Twinmz_families_and_IDs <- data.frame(Twinmz_families_and_IDs[1], stack(Twinmz_families_and_IDs[2:ncol(Twinmz_families_and_IDs)]))
Twinmz_families_and_IDs$Participant.Id <- Twinmz_families_and_IDs$values
Twinmz_families_and_IDs <- Twinmz_families_and_IDs %>% select(Rare.Diseases.Family.Id, Participant.Id)

# Consent - full siblings family
metatable_extraction <- read.delim("participant_2022-06-14_15-37-31.tsv")
Consent_extraction <- metatable_extraction %>% select(Participant.Id, Date.Of.Consent)
Fullsib_families_and_IDs_Consent <- left_join(Fullsib_families_and_IDs, Consent_extraction, by = "Participant.Id")
Fullsib_families_and_IDs_Consent <- Fullsib_families_and_IDs_Consent %>% unique()
Fullsib_families_and_IDs_Consent <- Fullsib_families_and_IDs_Consent %>% filter(!grepl("N/A", Date.Of.Consent))
Fullsib_families_and_IDs_Consent <- Fullsib_families_and_IDs_Consent[complete.cases(Fullsib_families_and_IDs_Consent), ]
# Consent - twins dizygous family
Twindz_families_and_IDs_Consent <- left_join(Twindz_families_and_IDs, Consent_extraction, by = "Participant.Id")
Twindz_families_and_IDs_Consent <- Twindz_families_and_IDs_Consent %>% unique()
Twindz_families_and_IDs_Consent <- Twindz_families_and_IDs_Consent %>% filter(!grepl("N/A", Date.Of.Consent))
Twindz_families_and_IDs_Consent <- Twindz_families_and_IDs_Consent[complete.cases(Twindz_families_and_IDs_Consent), ]
# Consent - twins monozygous family
Twinmz_families_and_IDs_Consent <- left_join(Twinmz_families_and_IDs, Consent_extraction, by = "Participant.Id")
Twinmz_families_and_IDs_Consent <- Twinmz_families_and_IDs_Consent %>% unique()
Twinmz_families_and_IDs_Consent <- Twinmz_families_and_IDs_Consent %>% filter(!grepl("N/A", Date.Of.Consent))
Twinmz_families_and_IDs_Consent <- Twinmz_families_and_IDs_Consent[complete.cases(Twinmz_families_and_IDs_Consent), ]

# YOB - full siblings family
YOB_extraction <- metatable_extraction %>% select(Participant.Id, Year.Of.Birth)
Fullsib_families_and_IDs_Consent_YOB <- left_join(Fullsib_families_and_IDs_Consent, YOB_extraction, by = "Participant.Id")
Fullsib_families_and_IDs_Consent_YOB <- Fullsib_families_and_IDs_Consent_YOB %>% unique()
Fullsib_families_and_IDs_Consent_YOB <- Fullsib_families_and_IDs_Consent_YOB %>% filter(!grepl("N/A", Year.Of.Birth))
Fullsib_families_and_IDs_Consent_YOB <- Fullsib_families_and_IDs_Consent_YOB[complete.cases(Fullsib_families_and_IDs_Consent_YOB), ]
# YOB - twin dizygous family
Twindz_families_and_IDs_Consent_YOB <- left_join(Twindz_families_and_IDs_Consent, YOB_extraction, by = "Participant.Id")
Twindz_families_and_IDs_Consent_YOB <- Twindz_families_and_IDs_Consent_YOB %>% unique()
Twindz_families_and_IDs_Consent_YOB <- Twindz_families_and_IDs_Consent_YOB %>% filter(!grepl("N/A", Year.Of.Birth))
Twindz_families_and_IDs_Consent_YOB <- Twindz_families_and_IDs_Consent_YOB[complete.cases(Twindz_families_and_IDs_Consent_YOB), ]
# YOB - twin monozygous family
Twinmz_families_and_IDs_Consent_YOB <- left_join(Twinmz_families_and_IDs_Consent, YOB_extraction, by = "Participant.Id")
Twinmz_families_and_IDs_Consent_YOB <- Twinmz_families_and_IDs_Consent_YOB %>% unique()
Twinmz_families_and_IDs_Consent_YOB <- Twinmz_families_and_IDs_Consent_YOB %>% filter(!grepl("N/A", Year.Of.Birth))
Twinmz_families_and_IDs_Consent_YOB <- Twinmz_families_and_IDs_Consent_YOB[complete.cases(Twinmz_families_and_IDs_Consent_YOB), ]

# Disease group - full siblings family
disease_metatable_extraction <- read.delim("rare_diseases_participant_dise_2021-10-12_10-39-17.tsv")
Disease_extraction <- disease_metatable_extraction %>% select(Participant.Id, Normalised.Disease.Sub.Group)
Fullsib_families_and_IDs_Consent_YOB_Disease <- left_join(Fullsib_families_and_IDs_Consent_YOB, Disease_extraction, by = "Participant.Id")
Fullsib_families_and_IDs_Consent_YOB_Disease <- Fullsib_families_and_IDs_Consent_YOB_Disease %>% unique()
Findduplicates <- Fullsib_families_and_IDs_Consent_YOB_Disease$Participant.Id[duplicated(Fullsib_families_and_IDs_Consent_YOB_Disease$Participant.Id)]
Findduplicates <- as.data.frame(Findduplicates)
Findduplicates$Participant.Id <- Findduplicates$Findduplicates
Findduplicates <- Findduplicates %>% select(Participant.Id)
Duplicated_info <- left_join(Findduplicates, Fullsib_families_and_IDs_Consent_YOB_Disease, by = "Participant.Id")
# More than one disease associated with an individual
Individual_diseases_merge <- Duplicated_info %>% 
  dplyr::group_by(Participant.Id) %>%
  dplyr::summarise(Normalised.Disease.Sub.Group = paste(Normalised.Disease.Sub.Group, collapse = ","))
Fullsib_families_and_IDs_Consent_YOB_Disease_x <- left_join(Fullsib_families_and_IDs_Consent_YOB, Individual_diseases_merge, by = "Participant.Id")
Fullsib_families_and_IDs_Consent_YOB_Disease_x$Normalised.Disease.Sub.Group <- as.character(Fullsib_families_and_IDs_Consent_YOB_Disease_x$Normalised.Disease.Sub.Group)
Fullsib_families_and_IDs_Consent_YOB_Disease <- rbind(Fullsib_families_and_IDs_Consent_YOB_Disease_x, Fullsib_families_and_IDs_Consent_YOB_Disease, by = "Rare.Diseases.Family.Id")
Fullsib_families_and_IDs_Consent_YOB_Disease <- Fullsib_families_and_IDs_Consent_YOB_Disease %>%
  group_by(Participant.Id) %>%
  slice(which.max(!is.na(Normalised.Disease.Sub.Group)))
# extra row checking why
check<-as.data.frame(Fullsib_families_and_IDs_Consent_YOB_Disease[which(!Fullsib_families_and_IDs_Consent_YOB_Disease$Participant.Id %in% Fullsib_families_and_IDs_Consent_YOB$Participant.Id),])
Fullsib_families_and_IDs_Consent_YOB_Disease <- Fullsib_families_and_IDs_Consent_YOB_Disease[!(Fullsib_families_and_IDs_Consent_YOB_Disease$Rare.Diseases.Family.Id=="Rare.Diseases.Family.Id"),]
Fullsib_families_and_IDs_Consent_YOB_Disease$Normalised.Disease.Sub.Group[is.na(Fullsib_families_and_IDs_Consent_YOB_Disease$Normalised.Disease.Sub.Group)] <- "Not reported"
Fullsib_families_and_IDs_Consent_YOB_Disease <- Fullsib_families_and_IDs_Consent_YOB_Disease %>% ungroup(Participant.Id)
# Disease group - twin dizygous family
Twindz_families_and_IDs_Consent_YOB_Disease <- left_join(Twindz_families_and_IDs_Consent_YOB, Disease_extraction, by = "Participant.Id")
Twindz_families_and_IDs_Consent_YOB_Disease <- Twindz_families_and_IDs_Consent_YOB_Disease %>% unique()
Findduplicates <- Twindz_families_and_IDs_Consent_YOB_Disease$Participant.Id[duplicated(Twindz_families_and_IDs_Consent_YOB_Disease$Participant.Id)]
Findduplicates <- as.data.frame(Findduplicates)
Findduplicates$Participant.Id <- Findduplicates$Findduplicates
Findduplicates <- Findduplicates %>% select(Participant.Id)
Duplicated_info <- left_join(Findduplicates, Twindz_families_and_IDs_Consent_YOB_Disease, by = "Participant.Id")
# More than one disease associated with an individual
Individual_diseases_merge <- Duplicated_info %>% 
  dplyr::group_by(Participant.Id) %>%
  dplyr::summarise(Normalised.Disease.Sub.Group = paste(Normalised.Disease.Sub.Group, collapse = ","))
Twindz_families_and_IDs_Consent_YOB_Disease_x <- left_join(Twindz_families_and_IDs_Consent_YOB, Individual_diseases_merge, by = "Participant.Id")
Twindz_families_and_IDs_Consent_YOB_Disease_x$Normalised.Disease.Sub.Group <- as.character(Twindz_families_and_IDs_Consent_YOB_Disease_x$Normalised.Disease.Sub.Group)
Twindz_families_and_IDs_Consent_YOB_Disease <- rbind(Twindz_families_and_IDs_Consent_YOB_Disease_x, Twindz_families_and_IDs_Consent_YOB_Disease, by = "Rare.Diseases.Family.Id")
Twindz_families_and_IDs_Consent_YOB_Disease <- Twindz_families_and_IDs_Consent_YOB_Disease %>%
  group_by(Participant.Id) %>%
  slice(which.max(!is.na(Normalised.Disease.Sub.Group)))
# extra row checking why
check<-as.data.frame(Twindz_families_and_IDs_Consent_YOB_Disease[which(!Twindz_families_and_IDs_Consent_YOB_Disease$Participant.Id %in% Twindz_families_and_IDs_Consent_YOB$Participant.Id),])
Twindz_families_and_IDs_Consent_YOB_Disease <- Twindz_families_and_IDs_Consent_YOB_Disease[!(Twindz_families_and_IDs_Consent_YOB_Disease$Rare.Diseases.Family.Id=="Rare.Diseases.Family.Id"),]
Twindz_families_and_IDs_Consent_YOB_Disease$Normalised.Disease.Sub.Group[is.na(Twindz_families_and_IDs_Consent_YOB_Disease$Normalised.Disease.Sub.Group)] <- "Not reported"
Twindz_families_and_IDs_Consent_YOB_Disease <- Twindz_families_and_IDs_Consent_YOB_Disease %>% ungroup(Participant.Id)
# Disease group - twin monozygous family
Twinmz_families_and_IDs_Consent_YOB_Disease <- left_join(Twinmz_families_and_IDs_Consent_YOB, Disease_extraction, by = "Participant.Id")
Twinmz_families_and_IDs_Consent_YOB_Disease <- Twinmz_families_and_IDs_Consent_YOB_Disease %>% unique()
Findduplicates <- Twinmz_families_and_IDs_Consent_YOB_Disease$Participant.Id[duplicated(Twinmz_families_and_IDs_Consent_YOB_Disease$Participant.Id)]
Findduplicates <- as.data.frame(Findduplicates)
Findduplicates$Participant.Id <- Findduplicates$Findduplicates
Findduplicates <- Findduplicates %>% select(Participant.Id)
Duplicated_info <- left_join(Findduplicates, Twinmz_families_and_IDs_Consent_YOB_Disease, by = "Participant.Id")
# More than one disease associated with an individual
Individual_diseases_merge <- Duplicated_info %>% 
  dplyr::group_by(Participant.Id) %>%
  dplyr::summarise(Normalised.Disease.Sub.Group = paste(Normalised.Disease.Sub.Group, collapse = ","))
Twinmz_families_and_IDs_Consent_YOB_Disease_x <- left_join(Twinmz_families_and_IDs_Consent_YOB_Disease, Individual_diseases_merge, by = "Participant.Id")
Twinmz_families_and_IDs_Consent_YOB_Disease_x$Normalised.Disease.Sub.Group <- as.character(Twinmz_families_and_IDs_Consent_YOB_Disease_x$Normalised.Disease.Sub.Group)
Twinmz_families_and_IDs_Consent_YOB_Disease <- rbind(Twinmz_families_and_IDs_Consent_YOB_Disease_x, Twinmz_families_and_IDs_Consent_YOB_Disease, by = "Rare.Diseases.Family.Id")
Twinmz_families_and_IDs_Consent_YOB_Disease <- Twinmz_families_and_IDs_Consent_YOB_Disease %>%
  group_by(Participant.Id) %>%
  slice(which.max(!is.na(Normalised.Disease.Sub.Group)))
# extra row checking why
check<-as.data.frame(Twinmz_families_and_IDs_Consent_YOB_Disease[which(!Twinmz_families_and_IDs_Consent_YOB_Disease$Participant.Id %in% Twinmz_families_and_IDs_Consent_YOB$Participant.Id),])
Twinmz_families_and_IDs_Consent_YOB_Disease <- Twinmz_families_and_IDs_Consent_YOB_Disease[!(Twinmz_families_and_IDs_Consent_YOB_Disease$Rare.Diseases.Family.Id=="Rare.Diseases.Family.Id"),]
Twinmz_families_and_IDs_Consent_YOB_Disease$Normalised.Disease.Sub.Group[is.na(Twinmz_families_and_IDs_Consent_YOB_Disease$Normalised.Disease.Sub.Group)] <- "Not reported"
Twinmz_families_and_IDs_Consent_YOB_Disease <- Twinmz_families_and_IDs_Consent_YOB_Disease %>% ungroup(Participant.Id)

Fullsib_families_and_IDs_1 <- Quad_fullsibling_metastart %>% select(Rare.Diseases.Family.Id, ProbandID, MothersID, FathersID, FullsiblingID)
Fullsib_families_and_IDs_1 <- data.frame(Fullsib_families_and_IDs_1[1], stack(Fullsib_families_and_IDs_1[2:ncol(Fullsib_families_and_IDs_1)]))
Fullsib_families_and_IDs_1$Participant.Id <- Fullsib_families_and_IDs_1$values
Fullsib_families_and_IDs_1 <- Fullsib_families_and_IDs_1 %>% select(Participant.Id, ind)
Unstack_families <- unstack(Fullsib_families_and_IDs_1)
Unstack_families$ProbandID <- as.character(Unstack_families$ProbandID)
Unstack_families$MothersID <- as.character(Unstack_families$MothersID)
Unstack_families$FathersID <- as.character(Unstack_families$FathersID)
Unstack_families$FullsiblingID <- as.character(Unstack_families$FullsiblingID)
Fullsib_families_and_IDs_Consent_YOB_Disease$ProbandID <- Fullsib_families_and_IDs_Consent_YOB_Disease$Participant.Id
Join_Proband <- left_join(Fullsib_families_and_IDs_Consent_YOB_Disease, Unstack_families, by = "ProbandID")
Join_Proband <- Join_Proband[complete.cases(Join_Proband), ]
Join_Proband$ProbandDate.Of.Consent <- Join_Proband$Date.Of.Consent
Join_Proband$ProbandYear.Of.Birth <- Join_Proband$Year.Of.Birth
Join_Proband$ProbandDisease.Group <- Join_Proband$Normalised.Disease.Sub.Group
Join_Proband <- Join_Proband %>% select(Rare.Diseases.Family.Id, ProbandID, ProbandDate.Of.Consent, ProbandYear.Of.Birth, ProbandDisease.Group)
Fullsib_families_and_IDs_Consent_YOB_Disease$MothersID <- Fullsib_families_and_IDs_Consent_YOB_Disease$Participant.Id
Join_Mother <- left_join(Fullsib_families_and_IDs_Consent_YOB_Disease, Unstack_families, by = "MothersID")
Join_Mother <- Join_Mother[complete.cases(Join_Mother), ]
Join_Mother$MothersDate.Of.Consent <- Join_Mother$Date.Of.Consent
Join_Mother$MothersYear.Of.Birth <- Join_Mother$Year.Of.Birth
Join_Mother$MothersDisease.Group <- Join_Mother$Normalised.Disease.Sub.Group
Join_Mother <- Join_Mother %>% select(Rare.Diseases.Family.Id, MothersID, MothersDate.Of.Consent, MothersYear.Of.Birth, MothersDisease.Group)
Fullsib_families_and_IDs_Consent_YOB_Disease$FathersID <- Fullsib_families_and_IDs_Consent_YOB_Disease$Participant.Id
Join_Father <- left_join(Fullsib_families_and_IDs_Consent_YOB_Disease, Unstack_families, by = "FathersID")
Join_Father <- Join_Father[complete.cases(Join_Father), ]
Join_Father$FathersDate.Of.Consent <- Join_Father$Date.Of.Consent
Join_Father$FathersYear.Of.Birth <- Join_Father$Year.Of.Birth
Join_Father$FathersDisease.Group <- Join_Father$Normalised.Disease.Sub.Group
Join_Father <- Join_Father %>% select(Rare.Diseases.Family.Id, FathersID, FathersDate.Of.Consent, FathersYear.Of.Birth, FathersDisease.Group)
Fullsib_families_and_IDs_Consent_YOB_Disease$FullsiblingID <- Fullsib_families_and_IDs_Consent_YOB_Disease$Participant.Id
Join_Fullsib <- left_join(Fullsib_families_and_IDs_Consent_YOB_Disease, Unstack_families, by = "FullsiblingID")
Join_Fullsib <- Join_Fullsib[complete.cases(Join_Fullsib), ]
Join_Fullsib$FullsiblingDate.Of.Consent <- Join_Fullsib$Date.Of.Consent
Join_Fullsib$FullsiblingYear.Of.Birth <- Join_Fullsib$Year.Of.Birth
Join_Fullsib$FullsiblingDisease.Group <- Join_Fullsib$Normalised.Disease.Sub.Group
Join_Fullsib <- Join_Fullsib %>% select(Rare.Diseases.Family.Id, FullsiblingID, FullsiblingDate.Of.Consent, FullsiblingYear.Of.Birth, FullsiblingDisease.Group)
Fullsib_overall_additional_meta <- left_join(Join_Proband, Join_Mother, by = "Rare.Diseases.Family.Id")
Fullsib_overall_additional_meta <- left_join(Fullsib_overall_additional_meta, Join_Father, by = "Rare.Diseases.Family.Id")
Fullsib_overall_additional_meta <- left_join(Fullsib_overall_additional_meta, Join_Fullsib, by = "Rare.Diseases.Family.Id")
Fullsib_overall_additional_meta <- left_join(Fullsib_overall_additional_meta, Quad_fullsibling_metastart, by = "Rare.Diseases.Family.Id")
Fullsib_overall_additional_meta <- Fullsib_overall_additional_meta %>% select(Rare.Diseases.Family.Id, ProbandID.x, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, ProbandDate.Of.Consent, ProbandYear.Of.Birth, ProbandDisease.Group, MothersID.x, MothersKey, MothersGenome, MothersGeneticsvsReported, MothersDate.Of.Consent, MothersYear.Of.Birth, MothersDisease.Group, FathersID.x, FathersKey, FathersGenome, FathersGeneticsvsReported, FathersDate.Of.Consent, FathersYear.Of.Birth, FathersDisease.Group, FullsiblingID.x, FullsiblingKey, FullsiblingGenome, FullsiblingGeneticsvsReported, FullsiblingDate.Of.Consent, FullsiblingYear.Of.Birth, FullsiblingDisease.Group)
# Generate columns for parent age at birth of child
Fullsib_overall_additional_meta$Mother_AgeatBirth <- with(Fullsib_overall_additional_meta, ProbandYear.Of.Birth - MothersYear.Of.Birth)
Fullsib_overall_additional_meta$Father_AgeatBirth <- with(Fullsib_overall_additional_meta, ProbandYear.Of.Birth - FathersYear.Of.Birth)
Fullsib_overall_additional_meta <- Fullsib_overall_additional_meta %>% select(Rare.Diseases.Family.Id, ProbandID.x, ProbandKey, ProbandGenome, ProbandGeneticsvsReported, ProbandDate.Of.Consent, ProbandYear.Of.Birth, ProbandDisease.Group, MothersID.x, MothersKey, MothersGenome, MothersGeneticsvsReported, MothersDate.Of.Consent, MothersYear.Of.Birth, Mother_AgeatBirth, MothersDisease.Group, FathersID.x, FathersKey, FathersGenome, FathersGeneticsvsReported, FathersDate.Of.Consent, FathersYear.Of.Birth, Father_AgeatBirth, FathersDisease.Group, FullsiblingID.x, FullsiblingKey, FullsiblingGenome, FullsiblingGeneticsvsReported, FullsiblingDate.Of.Consent, FullsiblingYear.Of.Birth, FullsiblingDisease.Group)
# Generate age at Sequencing
Fullsib_overall_additional_meta <- mutate(Fullsib_overall_additional_meta, ProbandSequencing.Year = gsub("-.*", "", Fullsib_overall_additional_meta$ProbandDate.Of.Consent))
Fullsib_overall_additional_meta <- mutate(Fullsib_overall_additional_meta, MothersSequencing.Year = gsub("-.*", "", Fullsib_overall_additional_meta$MothersDate.Of.Consent))
Fullsib_overall_additional_meta <- mutate(Fullsib_overall_additional_meta, FathersSequencing.Year = gsub("-.*", "", Fullsib_overall_additional_meta$FathersDate.Of.Consent))
Fullsib_overall_additional_meta <- mutate(Fullsib_overall_additional_meta, FullsiblingSequencing.Year = gsub("-.*", "", Fullsib_overall_additional_meta$FullsiblingDate.Of.Consent))
Fullsib_overall_additional_meta$ProbandSequencing.Year <- as.integer(Fullsib_overall_additional_meta$ProbandSequencing.Year)
Fullsib_overall_additional_meta$MothersSequencing.Year <- as.integer(Fullsib_overall_additional_meta$MothersSequencing.Year)
Fullsib_overall_additional_meta$FathersSequencing.Year <- as.integer(Fullsib_overall_additional_meta$FathersSequencing.Year)
Fullsib_overall_additional_meta$FullsiblingSequencing.Year <- as.integer(Fullsib_overall_additional_meta$FullsiblingSequencing.Year)
Fullsib_overall_additional_meta$ProbandAge.At.Sequencing <- with(Fullsib_overall_additional_meta, ProbandSequencing.Year - ProbandYear.Of.Birth)
Fullsib_overall_additional_meta$MothersAge.At.Sequencing <- with(Fullsib_overall_additional_meta, MothersSequencing.Year - MothersYear.Of.Birth)
Fullsib_overall_additional_meta$FathersAge.At.Sequencing <- with(Fullsib_overall_additional_meta, FathersSequencing.Year - FathersYear.Of.Birth)
Fullsib_overall_additional_meta$FullsiblingthersAge.At.Sequencing <- with(Fullsib_overall_additional_meta, FullsiblingSequencing.Year - FullsiblingsYear.Of.Birth)
Fullsib_overall_additional_meta <- Fullsib_overall_additional_meta %>% select(Rare.Diseases.Family.Id, ProbandID.x, ProbandPlatekey, ProbandGenomeBuild, ProbandGeneticsvsReported, ProbandDate.Of.Consent, ProbandAge.At.Sequencing, ProbandYear.Of.Birth, ProbandDisease.Group, MothersID.x, MothersPlatekey, MothersGenomeBuild, MothersGeneticsvsReported, MothersDate.Of.Consent, MothersAge.At.Sequencing, MothersYear.Of.Birth, Mother_AgeatBirth, MothersDisease.Group, FathersID.x, FathersPlatekey, FathersGenomeBuild, FathersGeneticsvsReported, FathersDate.Of.Consent, FathersAge.At.Sequencing, FathersYear.Of.Birth, Father_AgeatBirth, FathersDisease.Group, FullsiblingID.x, FullsiblingKey, FullsiblingGenome, FullsiblingGeneticsvsReported, FullsiblingDate.Of.Consent, FullsiblingthersAge.At.Sequencing, FullsiblingYear.Of.Birth, FullsiblingDisease.Group)


x <- repeatsizelist_allsiblingtypes %>% filter(ProbandGenome == "GRCh38", MothersGenome == "GRCh38", FathersGenome == "GRCh38", SiblingGenome == "GRCh38")
a <- x %>% select(MothersKey)
names(a)[names(a) == "MothersKey"] <- "Key"
b <- x %>% select(FathersKey)
names(b)[names(b) == "FathersKey"] <- "Key"
c <- x %>% select(ProbandKey)
names(c)[names(c) == "ProbandKey"] <- "Key"
d <- x %>% select(SiblingKey)
names(d)[names(d) == "SiblingKey"] <- "Key"
e <- do.call("rbind", list(a, b, c, d))
write.table(e, file = "/home/vgalassideforie/re_gecip/neurology/Valentina/Running_EHv3/platekeyfinal_38.txt", row.names=FALSE, sep="\t", quote = FALSE)

platekeyfilepaths <- read.delim("genome_file_paths_and_types_2021-10-13_10-34-19.tsv")
platekeyfilepaths <- platekeyfilepaths %>% filter(grepl('.bam', File.Path))
platekeyfilepaths <- platekeyfilepaths %>% select(Platekey, File.Path)
platekeyfilepaths_common <- as.data.frame(platekeyfilepaths[which(platekeyfilepaths$Platekey %in% e$Key),])
platekeyfilepaths_common <- platekeyfilepaths_common %>% select(File.Path)
write.table(platekeyfilepaths_common, file = "/home/vgalassideforie/re_gecip/neurology/Valentina/Running_EHv3/platekeyfilepaths_38.txt", row.names=FALSE, sep="\t", quote = FALSE)

x2 <- repeatsizelist_allsiblingtypes %>% filter(ProbandGenome == "GRCh37", MothersGenome == "GRCh37", FathersGenome == "GRCh37", SiblingGenome == "GRCh37")
a2 <- x2 %>% select(MothersKey)
names(a2)[names(a2) == "MothersKey"] <- "Key"
b2 <- x2 %>% select(FathersKey)
names(b2)[names(b2) == "FathersKey"] <- "Key"
c2 <- x2 %>% select(ProbandKey)
names(c2)[names(c2) == "ProbandKey"] <- "Key"
d2 <- x2 %>% select(SiblingKey)
names(d2)[names(d2) == "SiblingKey"] <- "Key"
e2 <- do.call("rbind", list(a2, b2, c2, d2))
write.table(e2, file = "/home/vgalassideforie/re_gecip/neurology/Valentina/Running_EHv3/platekeyfinal_37.txt", row.names=FALSE, sep="\t", quote = FALSE)

platekeyfilepaths <- read.delim("genome_file_paths_and_types_2021-10-13_10-34-19.tsv")
platekeyfilepaths <- platekeyfilepaths %>% filter(grepl('.bam', File.Path))
platekeyfilepaths <- platekeyfilepaths %>% select(Platekey, File.Path)
platekeyfilepaths_common2 <- as.data.frame(platekeyfilepaths[which(platekeyfilepaths$Platekey %in% e2$Key),])
platekeyfilepaths_common2 <- platekeyfilepaths_common2 %>% select(File.Path)
write.table(platekeyfilepaths_common2, file = "/home/vgalassideforie/re_gecip/neurology/Valentina/Running_EHv3/platekeyfilepaths_37.txt", row.names=FALSE, sep="\t", quote = FALSE)


e <- e %>% mutate(File_Name = paste(Key, "_HTT.vcf", sep = ""))
e <- e %>% select(File_Name)
write.table(e, file = "/home/vgalassideforie/re_gecip/neurology/Valentina/Running_EHv3/vcfnames_38.txt", row.names=FALSE, sep="\t", quote = FALSE)
vcfextraction_38 <- vcfextraction_38 %>% mutate(File_Name = paste(X, "_HTT.vcf", sep = ""))
vcfextraction_38 <- vcfextraction_38 %>% select(File_Name)
write.table(e, file = "/home/vgalassideforie/re_gecip/neurology/Valentina/Running_EHv3/platekey38_ransuccessfully.txt", row.names=FALSE, sep="\t", quote = FALSE)

bams <- bams %>% select(File.Path)
write.table(bams, file = "/home/vgalassideforie/re_gecip/neurology/Valentina/Running_EHv3/EH3_finalquad_output/bams.txt", row.names=FALSE, sep="\t", quote = FALSE)
