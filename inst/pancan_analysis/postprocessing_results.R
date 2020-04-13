require(ggplot2)

n_cancer_types <- 23

Num_Cell_Lines_Total <- 0
for(i in 1:n_cancer_types)
    Num_Cell_Lines_Total <- Num_Cell_Lines_Total + length(tc_results_separate_gene_sets$categorization_results_list_of_lists[[i]][[1]])

TC_Results_All_CLs_Cancer_Type_MSK10NN_Categorization <- matrix(nrow = Num_Cell_Lines_Total,ncol= 4)
colnames(TC_Results_All_CLs_Cancer_Type_MSK10NN_Categorization) <- c("Cell_Line_Name","Cancer_Type","MSK_10NN","Categorization")

j <- 0

for(i in 1:n_cancer_types)
{
  num_cell_lines_current <- length(tc_results_separate_gene_sets$categorization_results_list_of_lists[[i]][[1]])
  TC_Results_All_CLs_Cancer_Type_MSK10NN_Categorization[(j+1):(j+num_cell_lines_current),1] <- names(tc_results_separate_gene_sets$categorization_results_list_of_lists[[i]][[1]])
  TC_Results_All_CLs_Cancer_Type_MSK10NN_Categorization[(j+1):(j+num_cell_lines_current),2] <- rep(tc_results_separate_gene_sets$cancers_list[i],num_cell_lines_current)
  TC_Results_All_CLs_Cancer_Type_MSK10NN_Categorization[(j+1):(j+num_cell_lines_current),3] <- as.numeric(tc_results_separate_gene_sets$categorization_results_list_of_lists[[i]][[1]])
  TC_Results_All_CLs_Cancer_Type_MSK10NN_Categorization[(j+1):(j+num_cell_lines_current),4] <- tc_results_separate_gene_sets$categorization_results_list_of_lists[[i]][[3]][match(names(tc_results_separate_gene_sets$categorization_results_list_of_lists[[i]][[1]]),tc_results_separate_gene_sets$categorization_results_list_of_lists[[i]][[3]]$Sample_ID),2]
  j <- j + num_cell_lines_current
}

ggplot(as.data.frame(TC_Results_All_CLs_Cancer_Type_MSK10NN_Categorization), aes(x=Cancer_Type, y=as.numeric(MSK_10NN)/1000,color = Categorization)) + geom_point() + ylab("Mean Similarity to 10-nearest tumors") + xlab("Cancer Type of Cell Lines") + labs(title = "Similarity of Cell Lines to Tumors for 23 Cancer Types",subtitle = "679 COSMIC Cell Line Project Cell Lines compared to 7939 TCGA Tumors", caption = "Mean weighted similarity, using mutation, CNA and gene expression data")
