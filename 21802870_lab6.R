
BiocManager::install("biomaRt")
library(biomaRt)
globin_family <- read.csv("C:/Users/User/Desktop/class notes-hws/mbg326/globin_family.csv",sep=",")
head(globin_family)
approved_symbol = as.character(globin_family[,2])
approved_symbol
listEnsembl()

genes_mart = useEnsembl(biomart = "genes")
genes_mart

searchDatasets(mart = genes_mart,pattern = "hsapiens")
hsapiens_dataset = useDataset(dataset = "hsapiens_gene_ensembl",mart = 
                                genes_mart)

single_command_mart  = useEnsembl(biomart = "genes",dataset = "hsapiens_gene_ensembl")
single_command_mart

 

query=getBM(attributes = c("hgnc_symbol",    
            "ensembl_peptide_id","peptide"),filters = "hgnc_symbol",
            values = c("HBA1","HBA2","HBD"), 
            mart = single_command_mart)
query
hba1_peptide = as.character(query[1,])
hba2_peptide = as.character(query[4,])
hbd_peptide = as.character(query[10,])
data_frame = data.frame(rbind(hba1_peptide,hba2_peptide,hbd_peptide))
data_frame

BiocManager::install("Biostrings")
BiocManager::install("msa")
library(Biostrings)
library(msa)


data(PAM250) #### Import the PAM250 matrix 
data(package="Biostrings") #### See for other substitution matrices
substitution_matrix = "PAM250"
gap_open = 4
gap_extend = 1
alignment_type = "global"
alignment_1 = pairwiseAlignment(hba1_peptide[2],hbd_peptide[2],
                                            substitutionMatrix = substitution_matrix,
                                            gapOpening = gap_open,
                                            gapExtension = gap_extend,
                                            type = "global",
                                            scoreOnly = FALSE
)
alignment_1

alignment_2 = pairwiseAlignment(hba1_peptide[2],hbd_peptide[2],
                                substitutionMatrix = substitution_matrix,
                                gapOpening = gap_open,
                                gapExtension = gap_extend,
                                type = "local",
                                scoreOnly = FALSE
)
alignment_2




aa_Stringset = AAStringSet(c(hba1_peptide[2],hba2_peptide[2],hbd_peptide[2]),use.names = T)
aa_msa = msa(aa_Stringset)

aa_msa

print(aa_msa, show="complete")
