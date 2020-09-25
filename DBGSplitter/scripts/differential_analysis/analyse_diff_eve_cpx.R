# but = appliquer le modele de kissDE aux eve complexe = points de controle de l'epissage.
#On test les quantifs des arretes = k+1mers. on considere qu'a chaque noeuds 
#de degres sortant ou entrant >2 le spliceosome fait un choix.


library(Biobase)
library(aod)
library(DSS)
library(DESeq2)


################ varaiables globales #########################
nb_replicate = 2  #now it works with more than 2 replicates per condition. useless ?
nb_condition = 2
#rq : pour linstant on considere le meme nb de replicat par condition
nb_all_replicat = nb_condition * nb_replicate # le nombre total d'echantillion # = nbAll dans kissDE
##############################################################


##################### Lecture du fichier d'entree #################

read_quantif_file <- function(quantif_file_name) {
  
  quantifs = read.table(quantif_file_name, h=T)
  head(quantifs)
  
  # filter line when allmost one quantif < 1.
  #for the moment I do like in kissDE with the option filterLowCountsVariants.
  
  return(quantifs)
}

construct_sum_matrix_old<- function (quantif_df, schema_expe){
  matrix_sum = c()
  
  for (nb_line in 1:length(quantif_df$source)){
    row_to_add = c()
    
    for (cond in 1:4) { # pour chaque condition_replicat on ajoute une valeur a la ligne
      #on somme les quantifs pour chaque conditions_replicats
      sum = 0
      indice_values = c(cond+1, cond+5 ,cond+9, cond+13)
      real_indice_values = indice_values[1:quantif_df[nb_line,18]]
      for (indice in real_indice_values) {
        #print(c("indice = ",as.character(indice)))
        sum = sum + as.numeric(as.character(quantif_df[nb_line, indice])) 
        #print(c("to add = ",as.character(quantif_df[nb_line, indice])))
      }
      row_to_add = c(row_to_add, sum)
    }
    matrix_sum = rbind(matrix_sum, row_to_add)
  }
  rownames(matrix_sum)<- quantif_df$source
  return(matrix_sum)
}

construct_one_line_per_edge<- function(quantif_df, schema_expe){
  # we need for deSeq2 and to normalize the raw counts a matrix with one line per edge.
  # we need to keep the info of how many edges are involved in the event = nb_target
  
  one_edge_per_line_matrix = c()
  names_for_row=c()
  for (nb_line in 1:length(quantif_df$source)){
    matrix_add = c()
    
    #cumpute yhe indice of the column we need to sum. depends of the number of replicates(=len(schema_expe))
    #and of the number of target (=nb_target)
    indice_values = c()
    indice = 2 #first column we need is the col 2
    indice_values = c(indice_values,indice) 
    for (i in 1:(quantif_df[nb_line,"nb_targets"]-1)){ # -1 because we already put 2 for the first one
      indice = indice + length(schema_expe)
      indice_values = c(indice_values,indice)
    }
    #print(indice_values)
    
    for (i in indice_values){
      id_event = nb_line # source node
      names_for_row = c(names_for_row,as.character(quantif_df[nb_line,1]))
      vec = c(as.numeric(as.matrix(quantif_df[nb_line,i:(i+length(schema_expe)-1)])), id_event)
      names(vec) <- names(matrix_add)
      matrix_add = rbind(matrix_add, vec)
      
    } 
    one_edge_per_line_matrix = rbind(one_edge_per_line_matrix, matrix_add)
  }
  #print(names_for_row)
  colnames(one_edge_per_line_matrix) <- c(schema_expe, "event_id")
  rownames(one_edge_per_line_matrix) <- names_for_row 
  return(one_edge_per_line_matrix)
}

construct_sum_matrix <- function(one_line_per_edge_matrix, schema_expe){
  counts_sum_matrix = c()
  sv_row_same <- unique(rownames(one_line_per_edge_matrix)) 
  rownames(one_line_per_edge_matrix) = 1:length(one_line_per_edge_matrix[,1])
  one_line_per_edge_df = as.data.frame(one_line_per_edge_matrix)
  
  for (i in unique(one_line_per_edge_df$event_id)){
    event_df = one_line_per_edge_matrix[one_line_per_edge_df$event_id == i,1:length(schema_expe)]
    line_to_add = colSums(event_df)
    #print(line_to_add)
    counts_sum_matrix = rbind(counts_sum_matrix, line_to_add)
  }
  rownames(counts_sum_matrix) <- sv_row_same
  return(counts_sum_matrix)
}

one_edge_per_line_matrix_2_Test_df <- function(one_line_per_edge_matrix_Norm, schema_expe) {
  # We extract the needed data.frame to make the statistical test.
  
  one_line_per_edge_matrix_Norm = as.data.frame(one_line_per_edge_matrix_Norm)
  counts_df_list = list()
  i=1
  for (event in unique(one_line_per_edge_matrix_Norm[,length(schema_expe)+1])){
    counts_matrix = as.matrix(one_line_per_edge_matrix_Norm[one_line_per_edge_matrix_Norm$id_event == event,])
    
    counts = as.vector(t(counts_matrix[,1:length(schema_expe)])) #to read the matrix row by row
    condition = rep(1:nb_condition, each = nb_replicate, length = length(counts))
    replicat = rep(1:nb_replicate,length = length(counts))
    edge = rep(1:nrow(counts_matrix), each = length(schema_expe) ,length = length(counts))
    
    data_counts_df = as.data.frame(cbind(counts,condition, edge, replicat))
    names(data_counts_df) = c("counts","cond","edge","replicate")
    
    
    counts_df_list[[i]] <- data_counts_df
    i=i+1
  }
  
  # retourner une liste de data.frame pour ne lire qu'une seule fois le tableau.
  return(counts_df_list)
}


######################## functions to ana diff ###############################

counts_matrix_2_counts_df <- function(counts_matrix, nb_edge = 2, nb_condition = 2, nb_replicate = 2){
  #nb_edge = 3 # to define for each event : the number of edges we test. 
  counts = as.vector(t(counts_matrix)) #to read the matrix row by row
  condition = rep(1:nb_condition, each = nb_replicate, length = length(counts))
  replicat = rep(1:nb_replicate,length = length(counts))
  edge = rep(1:nb_edge, each = nb_all_replicat ,length = length(counts))
  
  data_counts_df = as.data.frame(cbind(counts,condition, edge, replicat))
  names(data_counts_df) = c("counts","cond","edge","replicate")
  
  return(data_counts_df)
}

normalization <- function(counts_matrix, one_line_per_edge_matrix, schema_expe) {

  
  # on fait la normalisation sur les comptage somme par edge directement.
  #group by condition and replicat in the matrix_counts
  #dataCountsEventSum = construct_sum_matrix(counts_matrix)
  #print(head(dataCountsEventSum))
  
  #Create DeSeq2 object for counts group by condition and replicate
  dds_sum_counts <- DESeqDataSetFromMatrix(
    countData= counts_matrix, 
    colData=data.frame(condition=schema_expe),
    design=~ condition)
  
  
  #cumpute sizeFactor
  dds_sum_counts <- estimateSizeFactors(dds_sum_counts)
  sizeFactorsSum <- sizeFactors(dds_sum_counts)
  
  #Create DeSeq2 object for raw counts
  # on the raw counts, without sum the counts for each cond, rep
  
  suppressMessages(dds_raw_counts <- DESeqDataSetFromMatrix(
  countData=one_line_per_edge_matrix[,1:length(schema_expe)], 
  colData=data.frame(condition=schema_expe),
  design=~ condition))
  
  #Apply sizes facors on raw counts
  sizeFactors(dds_raw_counts) <- sizeFactorsSum
  sizeFactors(dds_sum_counts) <- sizeFactorsSum
  # Create the matrix we will use to do the test.
  Norm_matrix_raw_counts = round(counts(dds_raw_counts, normalized= TRUE))
  #add event_id in the dataframe
  Norm_matrix_raw_counts = cbind(Norm_matrix_raw_counts, one_line_per_edge_matrix[,(length(schema_expe)+1)])
  colnames(Norm_matrix_raw_counts) <- c(schema_expe, "id_event")
  
  
  Norm_matrix_counts_sum = round(counts(dds_sum_counts, normalized= TRUE))
  # /!\ par defaut normalized = FALSE ! 
  
  return(list(Norm_matrix_raw_counts, Norm_matrix_counts_sum))
}

phi_estimation <- function(data_counts_Norm_sum, schema_expe){
  #calcul de phiDSS = phi local, parametre de sur-disperssion : 
  #dans les cas ou réplicats biologiques :
  #dataNormCountsEventSum : on somme les comptages pour chaque echantilons (C1R1, C2R2...)
  #une ligne par evenement et une colonne par echantillion et par chemin
  #dans KISSDE : c1r1UP c1r1LOW c1R2UP etc.

  designs <-schema_expe
  colnames(data_counts_Norm_sum) = rownames(designs)
  # Attention : il faut avoir le paquet biobase installe
  dispData <- newSeqCountSet(data_counts_Norm_sum, as.data.frame(designs), 
                             normalizationFactor = rep(1, length(schema_expe)))
  
  
  ## fix the seed to avoid the stochastic outputs of the 
  ## DSS:estDispersion function
  set.seed(40)
  dispData <- estDispersion(dispData)
  #calcul du phi pour chaque evenement
  phi_locals = dispersion(dispData)
  
  return(phi_locals)
}

test_diff <- function(phi_locals_foreach_event, data_counts_Norm, vector_of_sources){
  to_store_result = c()
  
  for (i in 1:length(phi_locals_foreach_event) ) {
    #print(i)
    phiDSS = phi_locals_foreach_event[i]
    event_id = i
    data_counts_Norm_df = data_counts_Norm[[i]]
    
    #modele additif 
    nbglmA <- negbin(counts~cond + edge, data=data_counts_Norm_df, random=~1, fixpar=list(4, phiDSS))
    #modele avec interaction
    nbglmI <- negbin(counts~cond * edge, data=data_counts_Norm_df, random=~1, fixpar=list(5, phiDSS))
    
    #le modele additif est imbriqué dans le modèle avec interaction => critere AIC
    #si modèle avec interaction est le meilleur modèle alors l'interaction est significative
    nbAnov <- anova(nbglmA, nbglmI)
    #pour récuperer la p-value : 
    p_value = nbAnov@anova.table$'P(> Chi2)'[2]
    #print(class(p_value))
    to_store_result = rbind(to_store_result,c(vector_of_sources[i],p_value))
  }
  #print(class(to_store_result[,2]))
  df_to_store_result = as.data.frame(to_store_result)

  names(df_to_store_result)=c("event_id","p_value")
  df_to_store_result$p_value = as.numeric(as.matrix(df_to_store_result$p_value))
  return(df_to_store_result)
}


cumpute_PSI <- function(counts_df_list, schema_expe){
  # counts df list is a list of data_frame.
  # each data frame contains the counts for one event.
  
  list_of_PSI = list()
  i=1 #indices in the list where we stock the PSI vectors
  for (event_counts in counts_df_list) {
    #print(i)
    PSI_cond1 = c()
    PSI_cond2 = c()
    
    #Cumputing sum for each replicate CiRi : 
    (sum_per_replicate_condition = aggregate(counts ~ replicate + cond, data=event_counts, FUN = sum))
  
    #Cumputing PSI for each replicate :
    PSI_cond1 = list()
    PSI_cond2 =list()
    for (rep in 1:nb_replicate){
      PSI_cond1_rep=c()
      PSI_cond2_rep=c()

      for (edge in unique(event_counts$edge)){
        
        #to filter low counts (ie sum of all edge <10 for one replicate): 
        sum_per_replicate_condition_c1 = sum_per_replicate_condition[sum_per_replicate_condition$replicate == rep & sum_per_replicate_condition$cond == 1,"counts"]
        if (sum_per_replicate_condition[sum_per_replicate_condition$replicate == rep & sum_per_replicate_condition$cond == 1,"counts"] <=  10) {
          PSI_cond1_rep=c(PSI_cond1_rep,NA)
        }else {
          PSI_cond1_rep = c(PSI_cond1_rep,(event_counts[event_counts$cond == 1 & event_counts$edge == edge & event_counts$replicate==rep,"counts"]/sum_per_replicate_condition[sum_per_replicate_condition$replicate == rep & sum_per_replicate_condition$cond == 1,"counts"]))
        }
        if (sum_per_replicate_condition[sum_per_replicate_condition$replicate == rep & sum_per_replicate_condition$cond == 2,"counts"] <= 10) {
          PSI_cond2_rep=c(PSI_cond2_rep,NA)
        }else {
          PSI_cond2_rep = c(PSI_cond2_rep,(event_counts[event_counts$cond == 2 & event_counts$edge == edge & event_counts$replicate==rep,"counts"]/sum_per_replicate_condition[sum_per_replicate_condition$replicate == rep & sum_per_replicate_condition$cond == 2,"counts"]))
        }
      }
    # vector of vector of psi for each condition
     PSI_cond1 = append(PSI_cond1,list(PSI_cond1_rep))
     PSI_cond2= append(PSI_cond2, list(PSI_cond2_rep))
    }    
    # to cumpute the mean PSI for each condition.
    #construct data.frame to make a colmean : 
    #on row per edge and one line per replicate, we want a mean PSI for each edge.
    df_psi_cond1 =c()
    for (rep in unique(event_counts$rep)){
      df_psi_cond1 = cbind(df_psi_cond1,PSI_cond1[[rep]])
    }
    
    
    df_psi_cond2 =c()
    for (rep in unique(event_counts$rep)){
      df_psi_cond2 = cbind(df_psi_cond2,PSI_cond2[[rep]])
    }
    
    #print(df_psi_cond2)
    
    #if more than a half of the psi are NA -> PSI = NA else we cumpute the psi without the NA
    #si on a plus de la moitie des psi=NA : 
    #cond1
    if (sum(is.na(df_psi_cond1[1,])) >= (length(df_psi_cond1[1,])/2) ) {
      PSI_cond1 = rowMeans(df_psi_cond1, na.rm=FALSE)
    }    else {#si on a moins de la moitie des psi=NA :
      PSI_cond1 = rowMeans(df_psi_cond1, na.rm=TRUE)
    }
    #cond2
    if (sum(is.na(df_psi_cond2[1,])) >= (length(df_psi_cond2[1,])/2) ) {
      PSI_cond2 = rowMeans(df_psi_cond2, na.rm=FALSE)
    } else { #si on a moins de la moitie des psi=NA
      PSI_cond2 = rowMeans(df_psi_cond2, na.rm=TRUE)
    }
    
    list_PSI_event = list(PSI_cond1, PSI_cond2)
    #print(list_PSI_event)
    
    list_of_PSI[[i]] = list_PSI_event
    i=i+1
  }
  return(list_of_PSI)
}




cumpute_delta_PSI<-function(PSI_for_each_event){
  # we generalise the delta psi for more than 2 values.
  # we simply do the substraction between the two vectors of psi for each event.
  # Nb : we have one vector of PSI for each condition. 
  #The length of the vectors corresponds to the number of edges we tested
  
  delta_PSI_list = list()
  i=1
  for (event_psi in PSI_for_each_event){
    delta_PSI = event_psi[[1]] - event_psi[[2]] 
    #print(delta_PSI)
    
    # we want all the length(delta_PSI) to be 4 
    if (length(delta_PSI) == 2) {
      delta_PSI = c(delta_PSI, NA, NA)
    } 
      
    if (length(delta_PSI) == 3) {
      delta_PSI = c(delta_PSI, NA)
    }
    delta_PSI_list[[i]] = delta_PSI  
    i=i+1
  }
  
  return(delta_PSI_list)
    
}

correct_p_values<-function(df_results_tests){
  
  df_results_tests$correct_p_value = p.adjust(df_results_tests$p_value, method = "fdr")
  
  return(df_results_tests)
}

write_output_file <- function(results_tests, delta_PSI_for_each_event,output_file_name){
  df_to_output = cbind(results_tests,as.matrix(delta_PSI_for_each_event))
  #names(df_to_output)<- c("event_id","p_value","delta_PSI_target1","delta_PSI_target2","delta_PSI_target3","delta_PSI_target4")
  #print(delta_PSI_for_each_event)
  write.table(x = as.matrix(df_to_output), file = output_file_name)
  
}