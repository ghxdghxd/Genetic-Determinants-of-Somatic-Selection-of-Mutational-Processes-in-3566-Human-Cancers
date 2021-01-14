# filter the mRNA expressions

```R
clin <-read.csv("/share/data4/TCGA/TCGA_The_Immune_Landscape_of_Cancer/TCGA-Clinical-Data-Resource.csv",header=T,sep="\t",stringsAsFactors=F,check.names=F)
rownames(clin) = clin$bcr_patient_barcode
clin = clin[clin$type %in% c("BRCA","COAD","ESCA","GBM","KIRC","LIHC","LUAD","OV","PRAD","STAD","THCA","UCEC"), ]
clin$clinical_stage[-grep("Stage", clin$clinical_stage)] = NA
clin$clinical_stage = gsub("Stage ", "", clin$clinical_stage)

clin$ajcc_pathologic_tumor_stage[-grep("Stage", clin$ajcc_pathologic_tumor_stage)] = NA
clin$ajcc_pathologic_tumor_stage = gsub("Stage ", "", clin$ajcc_pathologic_tumor_stage)

# An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics
# based on Table 1TCGA Pan-Cancer Cohort Characteristics, PRAD GBM without stage
# BRCA stage X ===> NA
clin$ajcc_pathologic_tumor_stage[which(clin$ajcc_pathologic_tumor_stage=="X")] = NA

clin$stage = NA
# AJCC stage for BRCA COAD ESCA KIRC LIHC LUAD STAD THCA
c_list = c("BRCA", "COAD", "ESCA", "KIRC", "LIHC", "LUAD", "STAD", "THCA")
clin$stage[clin$type %in% c_list] = clin$ajcc_pathologic_tumor_stage[clin$type %in% c_list]
# clinical stages for OV, UCEC
clin$stage[clin$type %in% c("OV","UCEC")] = clin$clinical_stage[clin$type %in% c("OV","UCEC")]
clin$stage = gsub("[ABC0-9]","",clin$stage)

clin$age_at_initial_pathologic_diagnosis = as.numeric(clin$age_at_initial_pathologic_diagnosis)

escc <- read.table("~/Projects/TCGA/TCGA.ESCC.txt", header = T, sep = "\t", stringsAsFactors = F)
clin$type[which(clin$bcr_patient_barcode %in% escc$barcode)] = 'ESCC'

purity = read.table("/share/data4/TCGA/TCGA_The_Immune_Landscape_of_Cancer/TCGA_mastercalls.abs_tables_JSedit.fixed.txt", 
    header = T, sep = "\t", stringsAsFactors = F)

filterSample = function(){
    Sample = purity$array
    TriSample = substr(Sample, 1, 12)

    DupTriSample = unique(TriSample[duplicated(TriSample)], 1,12)

    ID = Sample[-which(TriSample %in% DupTriSample)]

    DupSample = sapply(DupTriSample, function(x){list(Sample[grep(x, Sample)])})

    DupSample = sapply(DupSample, function(x){
        if(length(x) > 1){
            y = x[which(substr(x, 14,15) == "01")]
            if(length(y) < 1){
                return(x)
            }else{
                return(y)
            }
        }else{
            return(x)
        }
    })

    ID1 = as.character(DupSample[which(sapply(DupSample, length)==1)])
    DupSample = DupSample[-which(sapply(DupSample, length)==1)]
    return(c(ID, ID1, DupSample))
}

purity_sample = filterSample()
purity = purity[which(purity$array %in% purity_sample), ]
purity$ID = substr(purity$array, 1, 12)
clin$purity = purity$purity[match(clin$bcr_patient_barcode, purity$ID)]

save(clin, file = 'clin.RData')

```
