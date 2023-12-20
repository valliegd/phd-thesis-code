library(dplyr)
library(plyr)
library(tidyr)
library(stringr)
library(data.table)


# use masking to make individual columns for family members
complete_v7_table <- left_join(rare_disease_analysis_2021.04.16_11.49.15, rare_diseases_family_2021.04.16_11.50.39, by = "Rare.Diseases.Family.Id")
test <- complete_v7_table
test$is_proband = test$Biological.Relationship.To.Proband == "N/A"
test$is_mother = test$Biological.Relationship.To.Proband == "Mother"
test$is_father = test$Biological.Relationship.To.Proband == "Father"
test$is_proband <- as.character(test$is_proband)
test$is_proband[test$is_proband == "FALSE"] <- ""
test$is_mother <- as.character(test$is_mother)
test$is_mother[test$is_mother == "FALSE"] <- ""
test$is_father <- as.character(test$is_father)
test$is_father[test$is_father == "FALSE"] <- ""
test$is_proband[test$is_proband == "TRUE"] <- "Proband"
test$is_mother[test$is_mother == "TRUE"] <- "Mother"
test$is_father[test$is_father == "TRUE"] <- "Father"

test1 <- filter(test, Family.Group.Type == "Trio with Mother and Father" | is_proband == "Proband")
# test1 <- arrange(test, desc(rare_diseases_family_id))

Filterv7 <- test1 %>% select(Rare.Diseases.Family.Id, Family.Group.Type)
Filtertriov7 <- Filterv7 %>% filter(Family.Group.Type == "Trio with Mother and Father")
Filtertriouniquev7 <- Filtertriov7 %>% unique()

# filter probands in main and pilot programme
Mainprobandv7 <- complete_v7_table %>% filter(Biological.Relationship.To.Proband == "N/A")
Mainfilterv7 <- Mainprobandv7 %>% select(Rare.Diseases.Family.Id, Participant.Id)
Mainfilterv7$ProbandID = Mainfilterv7$Participant.Id

# filtering out mother and father from original table and joining to table with unique family_ids
Mothersv7 <- complete_v7_table %>% filter(Biological.Relationship.To.Proband == "Mother")
Mothersfilterv7 <- Mothersv7 %>% select(Rare.Diseases.Family.Id, Participant.Id)
Mothersfilterv7$MothersID = Mothersfilterv7$Participant.Id
Fathersv7 <- complete_v7_table %>% filter(Biological.Relationship.To.Proband == "Father")
Fathersfilterv7 <- Fathersv7 %>% select(Rare.Diseases.Family.Id, Participant.Id)
Fathersfilterv7$FathersID = Fathersfilterv7$Participant.Id

Joinmothersv7 = left_join(Filtertriouniquev7, Mothersfilterv7, by = "Rare.Diseases.Family.Id")
Joinparentsv7 = left_join(Joinmothersv7, Fathersfilterv7, by = "Rare.Diseases.Family.Id")
Joinallv7 = left_join(Joinparentsv7, Mainfilterv7, by = "Rare.Diseases.Family.Id")

# tidy up
Jointidyv7 <- Joinallv7 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, MothersID, FathersID)
Jointidy1v7 <- Jointidyv7 %>% unique()
# remove NAs from Mothers, Fathers and Probands
Joinblankv7 <- Jointidy1v7
Joinblankv7 <- sapply(Joinblankv7, as.character)
Joinblankv7[is.na(Joinblankv7)] <- " "
JoinNAv7 <- na.omit(Jointidy1v7)

# check platekeys for NA and left join
Platekeysv7 <- complete_v7_table %>% filter(Plate.Key == "N/A")
# 381 NAs in platekeys
test$platekey = test$Plate.Key
Main_Pv7 <- test %>% filter(Biological.Relationship.To.Proband == "N/A")
Main_P_filterv7 <- Main_Pv7 %>% select(Rare.Diseases.Family.Id, Plate.Key)
Main_P_filterv7$ProbandPlatekey = Main_P_filterv7$Plate.Key

Mothersfilter_platekeyv7 <- Mothersv7 %>% select(Rare.Diseases.Family.Id, Plate.Key)
Mothersfilter_platekeyv7$MothersPlatekey = Mothersfilter_platekeyv7$Plate.Key
Fathersfilter_platekeyv7 <- Fathersv7 %>% select(Rare.Diseases.Family.Id, Plate.Key)
Fathersfilter_platekeyv7$FathersPlatekey = Fathersfilter_platekeyv7$Plate.Key

Joinmothersplatekeyv7 = left_join(Jointidy1v7, Mothersfilter_platekeyv7, by = "Rare.Diseases.Family.Id")
Joinparentsplatekeyv7 = left_join(Joinmothersplatekeyv7, Fathersfilter_platekeyv7, by = "Rare.Diseases.Family.Id")
Joinallplatekeyv7 = left_join(Joinparentsplatekeyv7, Main_P_filterv7, by = "Rare.Diseases.Family.Id")

Jointidyplatekeyv7 <- Joinallplatekeyv7 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandPlatekey, MothersID, MothersPlatekey, FathersID, FathersPlatekey)
Jointidyplatekey1v7 <- Jointidyplatekeyv7[!duplicated(Jointidyplatekeyv7$Rare.Diseases.Family.Id), ]

# Find the rare disease family id which have na in any IDs or platekeys, create vector for participant id

filternav7 <- Jointidyplatekey1v7[rowSums(is.na(Jointidyplatekey1v7)) > 0, ]
listfilternav7 <- filternav7$Rare.Diseases.Family.Id
filternaprobandv7 <- Jointidyplatekey1v7[is.na(Jointidyplatekey1v7$ProbandPlatekey) , ]
listfilternaprobandv7 <- filternaprobandv7$Rare.Diseases.Family.Id
filternamotherv7 <- Jointidyplatekey1v7[is.na(Jointidyplatekey1v7$MothersPlatekey) , ]
listfilternamotherv7 <- filternamotherv7$Rare.Diseases.Family.Id
filternafatherv7 <- Jointidyplatekey1v7[is.na(Jointidyplatekey1v7$FathersPlatekey) , ]
listfilternafatherv7 <- filternafatherv7$Rare.Diseases.Family.Id

#Add genome build and genetics vs reported
Main_buildsv7 <- complete_v7_table %>% filter(Biological.Relationship.To.Proband == "N/A")
Main_builds_filterv7 <- Main_buildsv7 %>% select(Rare.Diseases.Family.Id, Genome.Build)
Main_builds_filterv7$ProbandGenomeBuild = Main_builds_filterv7$Genome.Build

Mothersfilter_buildsv7 <- Mothersv7 %>% select(Rare.Diseases.Family.Id, Genome.Build)
Mothersfilter_buildsv7$MothersGenomeBuild = Mothersfilter_buildsv7$Genome.Build
Fathersfilter_buildsv7 <- Fathersv7 %>% select(Rare.Diseases.Family.Id, Genome.Build)
Fathersfilter_buildsv7$FathersGenomeBuild = Fathersfilter_buildsv7$Genome.Build

Joinmothersbuildv7 = left_join(Jointidyplatekey1v7, Mothersfilter_buildsv7, by = "Rare.Diseases.Family.Id")
Joinparentsbuildv7 = left_join(Joinmothersbuildv7, Fathersfilter_buildsv7, by = "Rare.Diseases.Family.Id")
Joinallbuildv7 = left_join(Joinparentsbuildv7, Main_builds_filterv7, by = "Rare.Diseases.Family.Id")

Jointidybuildv7 <- Joinallbuildv7 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandPlatekey, ProbandGenomeBuild, MothersID, MothersPlatekey, MothersGenomeBuild, FathersID, FathersPlatekey, FathersGenomeBuild)
Jointidybuild1v7 <- Jointidybuildv7[!duplicated(Jointidybuildv7$Rare.Diseases.Family.Id), ]

Main_checksv7 <- complete_v7_table %>% filter(Biological.Relationship.To.Proband == "N/A")
Main_checks_filterv7 <- Main_checksv7 %>% select(Rare.Diseases.Family.Id, Genetic.Vs.Reported.Results)
Main_checks_filterv7$ProbandGeneticsvsReported = Main_checks_filterv7$Genetic.Vs.Reported.Results

Mothersfilter_checksv7 <- Mothersv7 %>% select(Rare.Diseases.Family.Id, Genetic.Vs.Reported.Results)
Mothersfilter_checksv7$MothersGeneticsvsReported = Mothersfilter_checksv7$Genetic.Vs.Reported.Results
Fathersfilter_checksv7 <- Fathersv7 %>% select(Rare.Diseases.Family.Id, Genetic.Vs.Reported.Results)
Fathersfilter_checksv7$FathersGeneticsvsReported = Fathersfilter_checksv7$Genetic.Vs.Reported.Results

Joinmotherscheckv7 = left_join(Jointidybuild1v7, Mothersfilter_checksv7, by = "Rare.Diseases.Family.Id")
Joinparentscheckv7 = left_join(Joinmotherscheckv7, Fathersfilter_checksv7, by = "Rare.Diseases.Family.Id")
Joinallcheckv7 = left_join(Joinparentscheckv7, Main_checks_filterv7, by = "Rare.Diseases.Family.Id")

Jointidycheckv7 <- Joinallcheckv7 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandPlatekey, ProbandGenomeBuild, ProbandGeneticsvsReported, MothersID, MothersPlatekey, MothersGenomeBuild, MothersGeneticsvsReported, FathersID, FathersPlatekey, FathersGenomeBuild, FathersGeneticsvsReported)
Jointidycheck1v7 <- Jointidycheckv7[!duplicated(Jointidycheckv7$Rare.Diseases.Family.Id), ]

#Proband genetics failed = N/A = 2034 (Ch37+N/A), NA's = 10, Failed = 156, Awaiting = 728, Passed = 9143, Missing = 2
#Mothers genetics failed = N/A = 2002, NA's = 16 Failed = 159, Awaiting = 753, Passed = 9141, Missing = 2
#Fathers genetics failed = N/A = 2018, NA's = 43, Failed = 158, Awaiting = 745, Passed = 9107, Missing = 2

#Filter for rows with all geneticsvsreported passed ADD PATHS AS WELL.

GRCh38_geneticspassedv7 <- subset(Jointidycheck1v7, Jointidycheck1v7$ProbandGeneticsvsReported == "familyPassesGvsRChecks" & Jointidycheck1v7$MothersGeneticsvsReported == "familyPassesGvsRChecks" & Jointidycheck1v7$FathersGeneticsvsReported == "familyPassesGvsRChecks")

#Extract repeat_size per platekey
x <- table_STR_repeat_size_each_row_allele_EHv3.2.2_HTT_simplified
x2 <- x %>% group_by(platekey) %>% dplyr::summarise(repeat_size_list = toString(repeat_size)) %>% as.data.frame()

#Want all repeat sizes for all trios, want all repeat sizes for just duos and want list of families with no repeat sizes for probands
test$platekey = test$Plate.Key
plate = left_join(x2, test, by = "platekey")
plate2 <- plate %>% select(Rare.Diseases.Family.Id, Family.Group.Type, Biological.Relationship.To.Proband, platekey, repeat_size_list, Participant.Id, Genome.Build, Genetic.Vs.Reported.Results)
plate2trio <- plate2 %>% filter(Genome.Build == "GRCh38", Family.Group.Type == "Trio with Mother and Father")
plate2triounique <- plate2trio[!duplicated(plate2trio$Rare.Diseases.Family.Id), ] 

Probandv7 <- plate2 %>% filter(Biological.Relationship.To.Proband == "N/A", Genome.Build == "GRCh38", Family.Group.Type == "Trio with Mother and Father", Genetic.Vs.Reported.Results == "familyPassesGvsRChecks")
Probandfilterv7 <- Probandv7 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, Biological.Relationship.To.Proband, platekey, repeat_size_list, Participant.Id, Genome.Build, Genetic.Vs.Reported.Results)
Probandfilterv7$ProbandID = Probandfilterv7$Participant.Id

Motherv7 <- plate2 %>% filter(Biological.Relationship.To.Proband == "Mother", Genome.Build == "GRCh38", Family.Group.Type == "Trio with Mother and Father", Genetic.Vs.Reported.Results == "familyPassesGvsRChecks")
Motherfilterv7 <- Motherv7 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, Biological.Relationship.To.Proband, platekey, repeat_size_list, Participant.Id, Genome.Build, Genetic.Vs.Reported.Results)
Motherfilterv7$MothersID = Motherfilterv7$Participant.Id
Fatherv7 <- plate2 %>% filter(Biological.Relationship.To.Proband == "Father", Genome.Build == "GRCh38", Family.Group.Type == "Trio with Mother and Father", Genetic.Vs.Reported.Results == "familyPassesGvsRChecks")
Fatherfilterv7 <- Fatherv7 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, Biological.Relationship.To.Proband, platekey, repeat_size_list, Participant.Id, Genome.Build, Genetic.Vs.Reported.Results)
Fatherfilterv7$FathersID = Fatherfilterv7$Participant.Id

# repeat trio with repeat size but add one by one
Main_repeatv7 <- plate2 %>% filter(Biological.Relationship.To.Proband == "N/A")
Main_repeat_filterv7 <- Main_repeatv7 %>% select(Rare.Diseases.Family.Id, repeat_size_list)
Main_repeat_filterv7$Probandrepeatsize = Main_repeat_filterv7$repeat_size_list

Mothers_repeatv7 <- plate2 %>% filter(Biological.Relationship.To.Proband == "Mother")
Mothersfilter_repeatv7 <- Mothers_repeatv7 %>% select(Rare.Diseases.Family.Id, repeat_size_list)
Mothersfilter_repeatv7$Mothersrepeatsize = Mothersfilter_repeatv7$repeat_size_list

Fathers_repeatv7 <- plate2 %>% filter(Biological.Relationship.To.Proband == "Father")
Fathersfilter_repeatv7 <- Fathers_repeatv7 %>% select(Rare.Diseases.Family.Id, repeat_size_list)
Fathersfilter_repeatv7$Fathersrepeatsize = Fathersfilter_repeatv7$repeat_size_list

Joinmothersrepeatv7 = left_join(GRCh38_geneticspassedv7, Mothersfilter_repeatv7, by = "Rare.Diseases.Family.Id")
Joinfathersrepeatv7 = left_join(GRCh38_geneticspassedv7, Fathersfilter_repeatv7, by = "Rare.Diseases.Family.Id")
Joinprobandrepeatv7 = left_join(GRCh38_geneticspassedv7, Main_repeat_filterv7, by = "Rare.Diseases.Family.Id")
Joinparentsrepeatv7 = left_join(Joinmothersrepeatv7, Fathersfilter_repeatv7, by = "Rare.Diseases.Family.Id")
Joinallrepeatv7 = left_join(Joinparentsrepeatv7, Main_repeat_filterv7, by = "Rare.Diseases.Family.Id")

Jointidyrepeatv7 <- Joinallrepeatv7 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandPlatekey, ProbandGenomeBuild, ProbandGeneticsvsReported, Probandrepeatsize, MothersID, MothersPlatekey, MothersGenomeBuild, MothersGeneticsvsReported, Mothersrepeatsize, FathersID, FathersPlatekey, FathersGenomeBuild, FathersGeneticsvsReported, Fathersrepeatsize)
Jointidyrepeat1v7 <- Jointidyrepeatv7[!duplicated(Jointidyrepeatv7$Rare.Diseases.Family.Id), ]
#8,976 entries

Joinmothersselectv7 <- Joinmothersrepeatv7 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, MothersID, MothersPlatekey, MothersGenomeBuild, MothersGeneticsvsReported, Mothersrepeatsize)
Joinmothersuniquev7 <- Joinmothersselectv7[!duplicated(Joinmothersselectv7$Rare.Diseases.Family.Id), ]
Joinfathersselectv7 <- Joinfathersrepeatv7 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, FathersID, FathersPlatekey, FathersGenomeBuild, FathersGeneticsvsReported, Fathersrepeatsize)
Joinfathersuniquev7 <- Joinfathersselectv7[!duplicated(Joinfathersselectv7$Rare.Diseases.Family.Id), ]
Joinprobandsselectv7 <- Joinprobandrepeatv7 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandPlatekey, ProbandGenomeBuild, ProbandGeneticsvsReported, Probandrepeatsize)
Joinprobandsuniquev7 <- Joinprobandsselectv7[!duplicated(Joinprobandsselectv7$Rare.Diseases.Family.Id), ]

#split so each allele in each column
probandsplitrepeats <- Joinprobandsuniquev7
splitrepeatsproband <- probandsplitrepeats %>% separate(Probandrepeatsize, c("Proband_allele_1", "Proband_allele_2"), ",")
mothersplitrepeats <- Joinmothersuniquev7
splitrepeatsmother <- mothersplitrepeats %>% separate(Mothersrepeatsize, c("Mother_allele_1", "Mother_allele_2"), ",")
fathersplitrepeats <- Joinfathersuniquev7
splitrepeatsfather <- fathersplitrepeats %>% separate(Fathersrepeatsize, c("Father_allele_1", "Father_allele_2"), ",")
Joinrepeatsizeparents <- left_join(splitrepeatsmother, splitrepeatsfather, by = "Rare.Diseases.Family.Id")
Joinrepeatsizeall <- left_join(Joinrepeatsizeparents, splitrepeatsproband, by = "Rare.Diseases.Family.Id")
Joinrepeatsizetidy <- Joinrepeatsizeall[!duplicated(Joinrepeatsizeall$Rare.Diseases.Family.Id), ]


# Interpreted build37 trios


complete_b37_table <- left_join(rare_disease_analysis_2021.06.08_10.42.53, rare_diseases_family_2021.06.08_10.42.09, by = "Rare.Diseases.Family.Id")

AriGRCh37 <- complete_b37_table %>% filter(Genome.Build == "GRCh37" & Family.Group.Type == "Trio with Mother and Father") %>% unique()
AriGRCh37platekeys <- AriGRCh37 %>% unique()
# 5696 genomes

Interpretetedplatekeysb37 <- rare_disease_interpreted_2021.06.08_16.03.49 %>% filter(Assembly == "GRCh37") %>% select(Plate.Key) %>% unique() %>% pull()
# 7754 genomes

TrioGRCh37Interpreted <- AriGRCh37platekeys %>% filter(Plate.Key %in% Interpretetedplatekeysb37)
# 4675 genomes



Proband37 <- TrioGRCh37Interpreted %>% filter(Biological.Relationship.To.Proband == "N/A")
Probandfilter37 <- Proband37 %>% select(Rare.Diseases.Family.Id, Participant.Id)
Probandfilter37$ProbandID = Probandfilter37$Participant.Id

# filtering out mother and father from original table and joining to table with unique family_ids
Mother37 <- TrioGRCh37Interpreted %>% filter(Biological.Relationship.To.Proband == "Mother")
Motherfilter37 <- Mother37 %>% select(Rare.Diseases.Family.Id, Participant.Id)
Motherfilter37$MothersID = Motherfilter37$Participant.Id
Father37 <- TrioGRCh37Interpreted %>% filter(Biological.Relationship.To.Proband == "Father")
Fatherfilter37 <- Father37 %>% select(Rare.Diseases.Family.Id, Participant.Id)
Fatherfilter37$FathersID = Fatherfilter37$Participant.Id

Joinmothers37 = left_join(TrioGRCh37Interpreted, Motherfilter37, by = "Rare.Diseases.Family.Id")
Joinparents37 = left_join(Joinmothers37, Fatherfilter37, by = "Rare.Diseases.Family.Id")
Joinall37 = left_join(Joinparents37, Probandfilter37, by = "Rare.Diseases.Family.Id")

# tidy up
Jointidy37 <- Joinall37 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, MothersID, FathersID)
Jointidy137 <- Jointidy37 %>% unique()

# 1580 genomes

ProbandKey37 <- Proband37 %>% select(Rare.Diseases.Family.Id, Plate.Key)
ProbandKey37$ProbandPlatekey = ProbandKey37$Plate.Key
MotherKey37 <- Mother37 %>% select(Rare.Diseases.Family.Id, Plate.Key)
MotherKey37$MothersPlatekey = MotherKey37$Plate.Key
FatherKey37 <- Father37 %>% select(Rare.Diseases.Family.Id, Plate.Key)
FatherKey37$FathersPlatekey = FatherKey37$Plate.Key

Joinmotherkey37 = left_join(Jointidy137, MotherKey37, by = "Rare.Diseases.Family.Id")
Joinparentkey37 = left_join(Joinmotherkey37, FatherKey37, by = "Rare.Diseases.Family.Id")
Joinallkey37 = left_join(Joinparentkey37, ProbandKey37, by = "Rare.Diseases.Family.Id")
Jointidykey37 <- Joinallkey37 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandPlatekey, MothersID, MothersPlatekey, FathersID, FathersPlatekey)
Jointidykey137 <- Jointidykey37[!duplicated(Jointidykey37$Rare.Diseases.Family.Id), ]

#Add genome build and genetics vs reported
ProbandBuild37 <- Proband37 %>% select(Rare.Diseases.Family.Id, Genome.Build)
ProbandBuild37$ProbandGenomeBuild = ProbandBuild37$Genome.Build
MotherBuild37 <- Mother37 %>% select(Rare.Diseases.Family.Id, Genome.Build)
MotherBuild37$MothersGenomeBuild = MotherBuild37$Genome.Build
FatherBuild37 <- Father37 %>% select(Rare.Diseases.Family.Id, Genome.Build)
FatherBuild37$FathersGenomeBuild = FatherBuild37$Genome.Build

Joinmothersbuild37 = left_join(Jointidykey137, MotherBuild37, by = "Rare.Diseases.Family.Id")
Joinparentsbuild37 = left_join(Joinmothersbuild37, FatherBuild37, by = "Rare.Diseases.Family.Id")
Joinallbuild37 = left_join(Joinparentsbuild37, ProbandBuild37, by = "Rare.Diseases.Family.Id")
Jointidybuild37 <- Joinallbuild37 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandPlatekey, ProbandGenomeBuild, MothersID, MothersPlatekey, MothersGenomeBuild, FathersID, FathersPlatekey, FathersGenomeBuild)
Jointidybuild137 <- Jointidybuild37[!duplicated(Jointidybuild37$Rare.Diseases.Family.Id), ]

ProbandCheck37 <- Proband37 %>% select(Rare.Diseases.Family.Id, Genetic.Vs.Reported.Results)
ProbandCheck37$ProbandGeneticsvsReported = ProbandCheck37$Genetic.Vs.Reported.Results
MotherCheck37 <- Mother37 %>% select(Rare.Diseases.Family.Id, Genetic.Vs.Reported.Results)
MotherCheck37$MothersGeneticsvsReported = MotherCheck37$Genetic.Vs.Reported.Results
FatherCheck37 <- Father37 %>% select(Rare.Diseases.Family.Id, Genetic.Vs.Reported.Results)
FatherCheck37$FathersGeneticsvsReported = FatherCheck37$Genetic.Vs.Reported.Results

Joinmotherscheck37 = left_join(Jointidybuild137, MotherCheck37, by = "Rare.Diseases.Family.Id")
Joinparentscheck37 = left_join(Joinmotherscheck37, FatherCheck37, by = "Rare.Diseases.Family.Id")
Joinallcheck37 = left_join(Joinparentscheck37, ProbandCheck37, by = "Rare.Diseases.Family.Id")

Jointidycheck37 <- Joinallcheck37 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandPlatekey, ProbandGenomeBuild, ProbandGeneticsvsReported, MothersID, MothersPlatekey, MothersGenomeBuild, MothersGeneticsvsReported, FathersID, FathersPlatekey, FathersGenomeBuild, FathersGeneticsvsReported)
Jointidycheck137 <- Jointidycheck37[!duplicated(Jointidycheck37$Rare.Diseases.Family.Id), ]

# 1580 GRCh37 trio are interpreted
# 13 Proband ID N/As, 40 Mother ID N/As, 38 Father ID N/As

TrioGRCh37completecases <- Jointidycheck137[complete.cases(Jointidycheck137), ]

# 1513 complete GRCh37 trio are interpreted (all rows with at least 1 NA removed)

x37 <- table_STR_repeat_size_each_row_allele_EHv3.2.2_HTT_simplified
y37 <- x37 %>% group_by(platekey) %>% dplyr::summarise(repeat_size_list = toString(repeat_size)) %>% as.data.frame()
complete_b37_table$platekey = complete_b37_table$Plate.Key
plate37 = left_join(y37, complete_b37_table, by = "platekey")
platey37 <- plate37 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, Biological.Relationship.To.Proband, platekey, repeat_size_list, Participant.Id, Genome.Build, Genetic.Vs.Reported.Results)

Proband_repeat37 <- platey37 %>% filter(Biological.Relationship.To.Proband == "N/A")
Proband_repeat_filter37 <- Proband_repeat37 %>% select(Rare.Diseases.Family.Id, repeat_size_list)
Proband_repeat_filter37$Probandrepeatsize = Proband_repeat_filter37$repeat_size_list

Mothers_repeat37 <- plate2 %>% filter(Biological.Relationship.To.Proband == "Mother")
Mothersfilter_repeat37 <- Mothers_repeat37 %>% select(Rare.Diseases.Family.Id, repeat_size_list)
Mothersfilter_repeat37$Mothersrepeatsize = Mothersfilter_repeat37$repeat_size_list

Fathers_repeat37 <- plate2 %>% filter(Biological.Relationship.To.Proband == "Father")
Fathersfilter_repeat37 <- Fathers_repeat37 %>% select(Rare.Diseases.Family.Id, repeat_size_list)
Fathersfilter_repeat37$Fathersrepeatsize = Fathersfilter_repeat37$repeat_size_list

Joinmothersrepeat37 = left_join(TrioGRCh37completecases, Mothersfilter_repeat37, by = "Rare.Diseases.Family.Id")
Joinfathersrepeat37 = left_join(TrioGRCh37completecases, Fathersfilter_repeat37, by = "Rare.Diseases.Family.Id")
Joinprobandrepeat37 = left_join(TrioGRCh37completecases, Proband_repeat_filter37, by = "Rare.Diseases.Family.Id")

Joinmothersselect37 <- Joinmothersrepeat37 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, MothersID, MothersPlatekey, MothersGenomeBuild, MothersGeneticsvsReported, Mothersrepeatsize)
Joinmothersunique37 <- Joinmothersselect37[!duplicated(Joinmothersselect37$Rare.Diseases.Family.Id), ]
Joinfathersselect37 <- Joinfathersrepeat37 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, FathersID, FathersPlatekey, FathersGenomeBuild, FathersGeneticsvsReported, Fathersrepeatsize)
Joinfathersunique37 <- Joinfathersselect37[!duplicated(Joinfathersselect37$Rare.Diseases.Family.Id), ]
Joinprobandsselect37 <- Joinprobandrepeat37 %>% select(Rare.Diseases.Family.Id, Family.Group.Type, ProbandID, ProbandPlatekey, ProbandGenomeBuild, ProbandGeneticsvsReported, Probandrepeatsize)
Joinprobandsunique37 <- Joinprobandsselect37[!duplicated(Joinprobandsselect37$Rare.Diseases.Family.Id), ]

#split so each allele in each column
probandsplitrepeats37 <- Joinprobandsunique37
splitrepeatsproband37 <- probandsplitrepeats37 %>% separate(Probandrepeatsize, c("Proband_allele_1", "Proband_allele_2"), ",")
mothersplitrepeats37 <- Joinmothersunique37
splitrepeatsmother37 <- mothersplitrepeats37 %>% separate(Mothersrepeatsize, c("Mother_allele_1", "Mother_allele_2"), ",")
fathersplitrepeats37 <- Joinfathersunique37
splitrepeatsfather37 <- fathersplitrepeats37 %>% separate(Fathersrepeatsize, c("Father_allele_1", "Father_allele_2"), ",")
Joinrepeatsizeparents37 <- left_join(splitrepeatsmother37, splitrepeatsfather37, by = "Rare.Diseases.Family.Id")
Joinrepeatsizeall37 <- left_join(Joinrepeatsizeparents37, splitrepeatsproband37, by = "Rare.Diseases.Family.Id")
Joinrepeatsizetidy37 <- Joinrepeatsizeall37[!duplicated(Joinrepeatsizeall37$Rare.Diseases.Family.Id), ]


# Join b38 and b37 trio data

Completetriodata <- dplyr::bind_rows(Joinrepeatsizetidy37, Joinrepeatsizetidy)
write.csv(Completetriodata, "~/re_gecip/shared_allGeCIPs/VG_DH_AT/Complete_trio_analysis/Completetriodata.csv")
write.csv(Completetriodata, "/mnt/ovd/slaveserver/o1630486472PUTmi/sharedFolder_5e06ef7d85ae976dbb6b59b83c35d909/re_gecip/Valentina/Complete_trio_analysis/Completetriodata.csv")

# 10,489 total genomes - COMPLETE DATA FOR B37 AND B38

# each allele own row
Alleleperrow <- Completetriodata %>% select(MothersPlatekey, FathersPlatekey, ProbandPlatekey, Proband_allele_1, Proband_allele_2, Mother_allele_1, Mother_allele_2, Father_allele_1, Father_allele_2)
Alleleperrow_analysis <- data.frame(Alleleperrow[1:3], stack(Alleleperrow[4:5]))
Alleleperrow_analysis1 <- data.frame(Alleleperrow[1:3], stack(Alleleperrow[6:7]))
Alleleperrow_analysis2 <- data.frame(Alleleperrow[1:3], stack(Alleleperrow[8:9]))
names(Alleleperrow_analysis)[4] <- "Proband"
names(Alleleperrow_analysis1)[4] <- "Mother"
names(Alleleperrow_analysis2)[4] <- "Father"
Alleleperrow_analysis$index <- seq.int(nrow(Alleleperrow_analysis))
Alleleperrow_analysis1$index <- seq.int(nrow(Alleleperrow_analysis1))
Alleleperrow_analysis2$index <- seq.int(nrow(Alleleperrow_analysis2))
Alleleperrow_analysis_join <- left_join(Alleleperrow_analysis, Alleleperrow_analysis1, by = "index")
Alleleperrow_analysis_join_1 <- left_join(Alleleperrow_analysis_join, Alleleperrow_analysis2, by = "index")
Alleleperrow_analysis_complete = select(Alleleperrow_analysis_join_1, MothersPlatekey, FathersPlatekey, ProbandPlatekey, Proband, Mother, Father)

# 20,978 entries

Complete_numeric <- transform(Completetriodata, Mother_allele_1 = as.numeric(Mother_allele_1), Mother_allele_2 = as.numeric(Mother_allele_2), Proband_allele_1 = as.numeric(Proband_allele_1), Proband_allele_2 = as.numeric(Proband_allele_2), Father_allele_1 = as.numeric(Father_allele_1), Father_allele_2 = as.numeric(Father_allele_2))
Mother_1_above35 <- Complete_numeric %>% filter(Complete_numeric$Mother_allele_1 > 35)
Mother_2_above35 <- Complete_numeric %>% filter(Complete_numeric$Mother_allele_2 > 35)
Proband_1_above35 <- Complete_numeric %>% filter(Complete_numeric$Proband_allele_1 > 35)
Proband_2_above35 <- Complete_numeric %>% filter(Complete_numeric$Proband_allele_2 > 35)
Father_1_above35 <- Complete_numeric %>% filter(Complete_numeric$Father_allele_1 > 35)
Father_2_above35 <- Complete_numeric %>% filter(Complete_numeric$Father_allele_2 > 35)

Mother_allele_35 <- Mother_2_above35 %>% select(MothersPlatekey)
Proband_allele_35 <- Proband_2_above35 %>% select(ProbandPlatekey)
Father_allele_35 <- Father_2_above35 %>% select(FathersPlatekey)
write.csv(Mother_allele_35, "/mnt/ovd/slaveserver/o1627901491GvDAB/sharedFolder_5e06ef7d85ae976dbb6b59b83c35d909/shared_allGeCIPs/VG_DH_AT/Complete_trio_analysis/Mother_allele_35.csv")
write.csv(Proband_allele_35, "/mnt/ovd/slaveserver/o1627901491GvDAB/sharedFolder_5e06ef7d85ae976dbb6b59b83c35d909/shared_allGeCIPs/VG_DH_AT/Complete_trio_analysis/Proband_allele_35.csv")
write.csv(Father_allele_35, "/mnt/ovd/slaveserver/o1627901491GvDAB/sharedFolder_5e06ef7d85ae976dbb6b59b83c35d909/shared_allGeCIPs/VG_DH_AT/Complete_trio_analysis/Father_allele_35.csv")

# Exclusion of families with questionable QC from pile-up plots from genomes_over_35.R
# Stack proband, mother and father 2 alleles against corresponding platekeys

Alleleperrow <- Alleleperrow[c(3,1,2,4,5,6,7,8,9)]
Alleleperrow_analysis_proband <- data.frame(Alleleperrow[1], stack(Alleleperrow[4:5]))
Alleleperrow_analysis_mother <- data.frame(Alleleperrow[2], stack(Alleleperrow[6:7]))
Alleleperrow_analysis_father <- data.frame(Alleleperrow[3], stack(Alleleperrow[8:9]))
names(Alleleperrow_analysis_proband)[1] <- "Platekey"
names(Alleleperrow_analysis_mother)[1] <- "Platekey"
names(Alleleperrow_analysis_father)[1] <- "Platekey"
names(Alleleperrow_analysis_proband)[2] <- "EH_Allele"
names(Alleleperrow_analysis_mother)[2] <- "EH_Allele"
names(Alleleperrow_analysis_father)[2] <- "EH_Allele"
Alleleperrow_analysis_proband <- Alleleperrow_analysis_proband %>% select(Platekey, EH_Allele)
Alleleperrow_analysis_mother <- Alleleperrow_analysis_mother %>% select(Platekey, EH_Allele)
Alleleperrow_analysis_father <- Alleleperrow_analysis_father %>% select(Platekey, EH_Allele)

Bind_complete_genomes_per_row <- rbind(Alleleperrow_analysis_proband, Alleleperrow_analysis_mother, Alleleperrow_analysis_father)

Completetriodata_plusQC <- full_join(Bind_complete_genomes_per_row, Alleleperrow_analysis_complete_over35, by = "Platekey")
names(Completetriodata_plusQC)[3] <- "EH_Allele_samples_over35"
names(Completetriodata_plusQC)[4] <- "QC_via_pileup"

Completetriodata_plusQC <- full_join(Bind_complete_genomes_per_row, Alleleperrow_analysis_complete_over35, by = "Platekey")
names(Completetriodata_plusQC)[3] <- "EH_Allele_samples_over35"
names(Completetriodata_plusQC)[4] <- "QC_via_pileup"

# QC added, need to add to original Completetriodata where each family has its own row