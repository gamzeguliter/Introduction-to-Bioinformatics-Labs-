install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)
listEnsembl()

genes_mart = useEnsembl(biomart = "genes") 
genes_mart

datasets_genes_mart = listDatasets(mart = genes_mart)
dim(datasets_genes_mart)

head(datasets_genes_mart)

searchDatasets(mart = genes_mart,pattern = "hsapiens")
hsapiens_dataset = useDataset(dataset = "hsapiens_gene_ensembl",mart = genes_mart)

hsapiens_dataset

single_command_mart  = useEnsembl(biomart = "genes",dataset = "hsapiens_gene_ensembl")
single_command_mart

filters = listFilters(single_command_mart) # We could have used the "hsapiens_dataset" Mart object as well 
filters[1:5,]

attributes = listAttributes(single_command_mart)
attributes[1:5,]
affyids = c("202763_at","209310_s_at","207500_at")
query=getBM(attributes = c('affy_hg_u133_plus_2', 'entrezgene_id',"gene_exon_intron","hgnc_symbol"),
            filters = 'affy_hg_u133_plus_2',
            values = affyids, 
            mart = single_command_mart)
head(query)

BiocManager::install("msa")
library(Biostrings)
library(msa)

axolotl_pax7 = "MACLPGAVPRMMRPGPGQNYPRTGFPLEGFAVSTPLGQGRVNQLGGVFINGRPLPNHVRHKIVEMAHHGIRPCVISRQLRVSHGCVSKILCRYQETGSIRPGAIGGSKPRQVATPDVEKK*"

hs_gene_biomart = useEnsembl(biomart = "genes",dataset = "hsapiens_gene_ensembl")
hs_gene_biomart_query = getBM(attributes = c("hgnc_symbol", 
                                             "ensembl_peptide_id","peptide"),filters = "hgnc_symbol",values = "PAX7",mart = hs_gene_biomart)
hs_gene_biomart_query
human_pax7 = as.character(hs_gene_biomart_query[1,2])

data(PAM250)
data(package="Biostrings") 
substitution_matrix = "PAM250"
gap_open = 4
gap_extend = 1
alignment_type = "global"

human_axoltl_aln_global = pairwiseAlignment(axolotl_pax7,human_pax7,
                                            substitutionMatrix = substitution_matrix,
                                            gapOpening = gap_open,
                                            gapExtension = gap_extend,
                                            type = alignment_type,
                                            scoreOnly = FALSE
)
human_axoltl_aln_global


writePairwiseAlignments(x = human_axoltl_aln_global,file = "")

pid(human_axoltl_aln_global,type = "PID1")
pid(human_axoltl_aln_global,type = "PID2")

human_pax4 = getBM(attributes = c("hgnc_symbol",    
                                  "ensembl_peptide_id","peptide"),filters = "hgnc_symbol",values = "PAX4",mart = hs_gene_biomart)
human_pax4 = as.character(human_pax4 [1,2])
hs_px4_px7_axol_msa
print(hs_px4_px7_axol_msa, show="complete")