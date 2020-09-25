#!/usr/bin/env Rscript


##### to run diff ana like a black box
# AIM : incluse analyse diff in the pipeline
#input = file.inputR, nb_replicate

#setwd("~/Documents/these/complex_events/Analyse_diff/repo_git_ana_diff/")
source("analyse_diff_eve_cpx.R")


#argument parser :
args = commandArgs(trailingOnly=TRUE)

print(args)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==3) {
  input_file <- args[1]
  nb_replicate <- args[2]
  output_file <- args[3]
}


#run diff analysis : 
print("reading input file...")
quantif_df = read_quantif_file(input_file)

schema_expe = c()
#we do it twice for 2 edges
for (i in seq_len(nb_condition)) {
  schema_expe <- c(schema_expe, rep(paste("Cond", i, sep="", collapse=""), nb_replicate))
}


one_line_per_edge_matrix = construct_one_line_per_edge(quantif_df, schema_expe)
counts_matrix = construct_sum_matrix(one_line_per_edge_matrix, schema_expe)
data_counts_Norm = normalization(counts_matrix, one_line_per_edge_matrix, schema_expe)
phi_locals = phi_estimation(data_counts_Norm[[2]], schema_expe) # sur les comptage normalises sommes par condition et replicats
one_line_per_edge_matrix_Norm = data_counts_Norm[[1]]
print("one line per edge matrix done")

counts_df_list = one_edge_per_line_matrix_2_Test_df(one_line_per_edge_matrix_Norm,schema_expe)
results_tests = test_diff(phi_locals, counts_df_list,rownames(counts_matrix)) # sur les comptages normalises NON sommes par condition et replicats
results_tests_cor = correct_p_values(results_tests)
PSI_for_each_event = cumpute_PSI(counts_df_list)
delta_PSI_for_each_event = cumpute_delta_PSI(PSI_for_each_event)

print("write output")
write_output_file(results_tests_cor, delta_PSI_for_each_event, output_file)
