rm(list = ls())

#number_minority = 4
#number_iteration = 22

args = commandArgs(trailingOnly=TRUE)
number_minority = as.integer(args[1])
number_iteration   = as.integer(args[2])
data_file_name = paste("data/replication_",number_iteration,"_cluster_",number_minority ,".csv", sep = "")

print(paste("number of minority group is:",number_minority))
print(paste("number of replication is:",number_iteration))

############################
# Functions
############################

# Function to sample n data from each category

#install.packages("scales",repos = "http://cran.us.r-project.org")
#install.packages("tidyverse",repos = "http://cran.us.r-project.org")
#install.packages("fastDummies",repos = "http://cran.us.r-project.org")
#install.packages("pROC",repos = "http://cran.us.r-project.org")
#install.packages("dplyr",repos = "http://cran.us.r-project.org")
#install.packages("lattice",repos = "http://cran.us.r-project.org")
#install.packages("fungible",repos = "http://cran.us.r-project.org")
#install.packages("cluster",repos = "http://cran.us.r-project.org")
#install.packages("fastcluster",repos = "http://cran.us.r-project.org")

library(scales)
library(tidyverse)
library(fastDummies)
library(pROC)
library(dplyr)
library(lattice)
#library(fungible)
library(cluster)

#library(RColorBrewer)
#library(ggplot2)

############################################
# Sampling functions
############################################

sample_index_from_each_cluster <- function(df, n, cluster_index) {
  # Create a new column of row indices
  df$row_indices <- 1:nrow(df)
  
  # Perform the sampling and extract row indices
  result_indices <- df %>%
    group_by_at(cluster_index) %>%
    sample_n(n) %>%
    pull(row_indices)
  
  return(result_indices)
}

sample_index_randomly <- function(df, n) {
  
  result_df_indices <- sample(1:nrow(df), size = n , replace = FALSE)
  
  return(result_df_indices)
}

# Function to sample n data from each category
sample_from_each_cluster <- function(df, n, cluster_index) {
  result_df <- df %>%
    group_by_at(cluster_index) %>%
    sample_n(n)
  
  return(as.data.frame(result_df))
}
sample_randomly <- function(df, n) {
  result_df <- df %>%
    sample_n(n)
  return(as.data.frame(result_df))
}
n_smallest <- function(v, a, b) {
  if (a > b || a > length(v) || b > length(v)) {
    stop("Invalid indices!")
  }
  
  sorted_vector <- sort(v, decreasing = FALSE)
  return(sorted_vector[a:b])
}

inverse_logit <- function(x) {
  1 / (1 + exp(-x))
}


########################################################
# Simulation explanatory variable(clustered) data
########################################################

simulate_cluster_data <- function(in.nvar,in.clus.size,in.eta,in.compact,seed_number){
  # in.var: number of variables
  # in.clus: number of clusters
  # BR: vector of cluster proportion
  # in.eta: Indicator validities
  # BigN: Groups sizes for Population
  
  in.nclus <- length(in.clus.size)
  sample.data <- monte(
    seed=seed_number,
    nvar=in.nvar,
    nclus = in.nclus,
    clus.size = in.clus.size,
    eta2 = in.eta,
    secor = .7,
    compactness=in.compact,
    sortMeans = TRUE,
    skew.list = list( rep(1,10),rep(1,10),rep(1,10),rep(0,10))
    #kurt.list = kurtosis_list
  )
  
  return(sample.data)
}

# Simulate datasets with 2 to 5 cluster structure

in.nvar <- 10  ##Number of variables
unbalanced_BR <- list(c(2500,2500,2500,2500),c(1000,3000,3000,3000),
                      c(1000,1000,4000,4000),c(1000,1000,1000,7000)
)

########################################################
# Simulation outcome data
########################################################

gen_event <- function(beta0,beta, df,sample_size){
  # This function generates the outcome event
  
  number_variable <- dim(df)[2]-1
  
  # generate dummy variable for cluster
  df <- dummy_cols(df,select_columns = c("id"),remove_selected_columns = FALSE )
  
  
  X_dummy <- as.matrix(df[,(number_variable+2):(number_variable+5)])
  
  # Pick the explanaory variable matirx
  X <- as.matrix(df[,2:(number_variable+1)])
  
  
  # Combine the original features and interaction terms
  X_prime <- cbind(X, X_dummy)
  
  # compute the event probability based on inverse logit function
  eta <- beta0 + X_prime%*%beta
  pi <- exp(eta )/(1+exp(eta))
  E <- rbinom(n=sample_size,size=1,prob=pi)
  #pi<-runif(n=sample_size, min = 0, max = 1)
  #E <- rbinom(n=sample_size,size=1,prob=0.5)
  result <- data.frame(event = E)
  
  return(cbind(result,df))
}

compute_minimum_distance <- function(data_set,true_cluster,dm,sample_index,min_person_a,min_person_b){
  cluster_index_list = list()
  n_cluster=4
  for (i in 1:n_cluster){
    cluster_index_list[[i]] <- which(true_cluster == i)
  }
  
  cluster_sample_distance <- list()
  for (j in 1:n_cluster){
    current_cluster_sample_distance <- list()
    for ( i in 1:length(cluster_index_list[[j]]) ){
      
      
      current_cluster_sample_distance[[i]] <- n_smallest(dm[cluster_index_list[[j]][i],sample_index],
                                                         min_person_a,min_person_b)
    }
    cluster_sample_distance[[j]] = Reduce("+", current_cluster_sample_distance)
  }
  
  cluster_result <- do.call("rbind", cluster_sample_distance)
  colnames(cluster_result)<- min_person_a:min_person_b
  return(colMeans(cluster_result) )
}

compute_cluster_sample_minimum_distance <- function(cluster_result,k,dm,simulated_data,sample_size,
                                                    n_replications,min_person_a,min_person_b){
  
  simulated_data$id <- cutree(cluster_result,k)
  cluster_index <- 2
  n_cluster <- length(unique(simulated_data$id))
  cluster_size <- ceiling(sample_size/n_cluster)
  simulated_data$unique_id <- 1:dim(simulated_data)[1]
  cluster_index_list = list()
  for (i in 1:n_cluster){
    cluster_index_list[[i]] <- which(simulated_data$id == i)
  }
  cluster_minimum_distance <- list()
  for (k in 1:n_replications){
    
    cluster_sample <- sample_from_each_cluster(simulated_data,cluster_size,cluster_index)
    cluster_samples_index <- match(cluster_sample$unique_id, simulated_data$unique_id)
    cluster_sample_distance <- list()
    
    for (j in 1:n_cluster){
      current_cluster_sample_distance <- list()
      for ( i in 1:length(cluster_index_list[[j]]) ){
        
        current_cluster_sample_distance[[i]] <- n_smallest(dm[cluster_index_list[[j]][i],cluster_samples_index],
                                                           min_person_a,min_person_b)
      }
      cluster_sample_distance[[j]] = Reduce("+", current_cluster_sample_distance)
    }
    
    cluster_minimum_distance[[k]] <- Reduce("+", cluster_sample_distance)
  }
  cluster_result <- do.call("rbind", cluster_minimum_distance)
  colnames(cluster_result)<- min_person_a:min_person_b
  return(colMeans(cluster_result[,c(1,5,10)]) )
}


logistic_model_train <- function(train_set,test_set,number_var){
  #This function train a logistic model based on train dataset and return AUC
  test_outcome <- test_set$event
  train_glm <- glm(event ~ ., data = train_set[,c(1,3:(2+number_var))], family = "binomial")
  predicted_probs <- predict(train_glm, newdata = test_set[,3:(2+number_var)], type = "response")
  test_set$predicted <- predicted_probs
  test_set$class <- ifelse(test_set$predicted > 0.5, 1, 0)
  test_set_1 <- test_set[which(test_set$id %in% c(1)),]
  test_set_2 <- test_set[which(test_set$id %in% c(2)),]
  test_set_3 <- test_set[which(test_set$id %in% c(3)),]
  test_set_4 <- test_set[which(test_set$id %in% c(4)),]
  acc_overall <- mean(test_set$class == test_set$event)
  acc_1 <-  mean(test_set_1$class == test_set_1$event)
  acc_2 <-  mean(test_set_2$class == test_set_2$event)
  acc_3 <-  mean(test_set_3$class == test_set_3$event)
  acc_4 <-  mean(test_set_4$class == test_set_4$event)
  
  auc_overall <- auc(roc(test_set$event, test_set$predicted))
  
  auc_1 <- auc(roc(test_set_1$event, test_set_1$predicted))
  auc_2 <- auc(roc(test_set_2$event, test_set_2$predicted))
  auc_3 <- auc(roc(test_set_3$event, test_set_3$predicted))
  auc_4 <- auc(roc(test_set_4$event, test_set_4$predicted))
  
  return(list( c(acc_overall,acc_1,acc_2,acc_3,acc_4),  c(auc_overall,auc_1,auc_2,auc_3,auc_4) ))
}

generate_AUC <- function(simulated_df,cluster_result,dm, training_size,number_var ){
  # step 1: decide the k
  n_replication = 50
  true_group <- simulated_df$id # save the true cluster label
  
  #test_result = data.frame(distance_1 = numeric(0),distance_5 = numeric(0),distance_10 = numeric(0) )
  
  # step 2: Compute the cluster sample chosen from cluster structure
  optimal_k <- c(2:6,10)
  
  random_1_10_distance_list <- list()
  random_distance_list <- list()
  random_train_index <- list()
  
  
  random_distance_vec<- c()
  cluster_distance_vec1<-c()
  cluster_1_10_distance_list1 <- list()
  cluster_distance_vec2<-c()
  cluster_1_10_distance_list2 <- list()
  cluster_distance_vec3<-c()
  cluster_1_10_distance_list3 <- list()
  cluster_distance_vec4<-c()
  cluster_1_10_distance_list4 <- list()
  cluster_distance_vec5<-c()
  cluster_1_10_distance_list5 <- list()
  cluster_distance_vec6<-c()
  cluster_1_10_distance_list6 <- list()
  cluster_distance_vec7<-c()
  cluster_1_10_distance_list7 <- list()
  #cluster_distance_vec8<-c()
  #cluster_1_10_distance_list8 <- list()
  
  cluster_train_index1 <- list()
  cluster_train_index2 <- list()
  cluster_train_index3 <- list()
  cluster_train_index4 <- list()
  cluster_train_index5 <- list()
  cluster_train_index6 <- list()
  cluster_train_index7 <- list()
  cluster_train_index8 <- list()
  
  
  for (i in 1:n_replication){
    train_index_list <- list()
    for (j in 1:(length(optimal_k)+2) ){
      current_k = optimal_k[j]
      # obtain the training index for cluster sample and random sample
      
      if (j == length(optimal_k)+2 ){
        
        #generate random sample index
        train_index_list[[j]] = sample_index_randomly(simulated_df,training_size)
        random_train_index[[i]] = train_index_list[[j]]
        
        #train_index_distance = compute_minimum_distance(simulated_df,true_group,dm,train_index_list[[j]],1,20)
        #random_distance_vec[i] = mean(compute_minimum_distance(simulated_df,true_group,dm,train_index_list[[j]],1,10))
        random_1_10_distance_list[[i]] = compute_minimum_distance(simulated_df,true_group,dm,train_index_list[[j]],1,10)
        random_distance_vec[i] = mean(random_1_10_distance_list[[i]] )
      }
      else if (j == length(optimal_k)+1 ){
        
        simulated_df$id <- true_group
        train_index_list[[j]] = sample_index_from_each_cluster(simulated_df,training_size/4,2)
        train_index_distance = compute_minimum_distance(simulated_df,true_group,dm,train_index_list[[j]],1,10)
        cluster_train_index7[[i]] = train_index_list[[j]]
        cluster_distance_vec7[i]= mean(train_index_distance)
        cluster_1_10_distance_list7[[i]] <- train_index_distance
      }else{
        
        #generate cluster sample index
        simulated_df$id <- cutree(cluster_result,current_k)
        train_index_list[[j]] = sample_index_from_each_cluster(simulated_df,training_size/current_k,2)
        train_index_distance = compute_minimum_distance(simulated_df,true_group,dm,train_index_list[[j]],1,10)
        if (j==1){cluster_train_index1[[i]] = train_index_list[[j]]
        cluster_distance_vec1[i]= mean(train_index_distance)
        cluster_1_10_distance_list1[[i]] <- train_index_distance
        }
        if (j==2){cluster_train_index2[[i]] = train_index_list[[j]]
        cluster_distance_vec2[i]= mean(train_index_distance)
        cluster_1_10_distance_list2[[i]] <- train_index_distance
        }
        if (j==3){cluster_train_index3[[i]] = train_index_list[[j]]
        cluster_distance_vec3[i]= mean(train_index_distance)
        cluster_1_10_distance_list3[[i]] <- train_index_distance
        }
        if (j==4){cluster_train_index4[[i]] = train_index_list[[j]]
        cluster_distance_vec4[i]= mean(train_index_distance)
        cluster_1_10_distance_list4[[i]] <- train_index_distance
        }
        if (j==5){cluster_train_index5[[i]] = train_index_list[[j]]
        cluster_distance_vec5[i]= mean(train_index_distance)
        cluster_1_10_distance_list5[[i]] <- train_index_distance
        }
        if (j==6){cluster_train_index6[[i]] = train_index_list[[j]]
        cluster_distance_vec6[i]= mean(train_index_distance)
        cluster_1_10_distance_list6[[i]] <- train_index_distance
        }
      }
    }
    
  }
  
  cluster_distance_list = c(cluster_distance_vec1,cluster_distance_vec2,cluster_distance_vec3,
                            cluster_distance_vec4,cluster_distance_vec5,cluster_distance_vec6)
  
  cluster_1_10_distance_list = c(cluster_1_10_distance_list1,cluster_1_10_distance_list2,cluster_1_10_distance_list3,
                                 cluster_1_10_distance_list4,cluster_1_10_distance_list5,cluster_1_10_distance_list6)
  
  cluster_train_index = c(cluster_train_index1,cluster_train_index2,cluster_train_index3,
                          cluster_train_index4,cluster_train_index5,cluster_train_index6)
  
  # obtain the train index for random sample,cluster sample, true sample
  optimal_cluster_index=cluster_train_index[[which.min(cluster_distance_list)]]
  optimal_true_index = cluster_train_index7[[which.min(cluster_distance_vec7)]]
  #worst_random_index = random_train_index[[which.max(random_distance_vec)]]
  random_sample_pick = sample(1:n_replication,1)
  worst_random_index = random_train_index[[random_sample_pick]]
  final_train_index_list = list(worst_random_index,optimal_cluster_index,optimal_true_index)
  
  # compute the mean of 1 to 10 nearest neighbor distance
  random_1_10_distance = random_1_10_distance_list[[random_sample_pick]]
  optimal_1_10_distance = cluster_1_10_distance_list[[which.min(cluster_distance_list)]]
  optimal_true_1_10_distance = cluster_1_10_distance_list7[[which.min(cluster_distance_vec7)]]
  # print information
  print(paste("random sample distance is:",random_distance_vec[ random_sample_pick],sep = " " ))
  print( paste("cluster sample distance is: ",min(cluster_distance_list), sep = " "))
  print( paste("true sample distance is: ",min(cluster_distance_vec7), sep = " "))
  
  print(paste("cluster index is: ",which.min(cluster_distance_list),sep = " "))
  print(paste( "number of cluster is :",ceiling(which.min(cluster_distance_list)/n_replication), sep = " " ) )
  
  distance_df = data.frame(random_sample_distance=random_distance_vec[ random_sample_pick],
                           cluster_sample_distance=min(cluster_distance_list),
                           true_sample_distance=min(cluster_distance_vec7))
  
  all_train_index = unique(unlist(final_train_index_list))
  
  simulated_df$id <- true_group
  df.test = simulated_df[-all_train_index,]
  
  # Obtain the AUC for all cluster sample model and random sample model
  current_acc <- matrix(0,nrow=0, ncol=5)
  current_auc <- matrix(0,nrow=0, ncol=5)
  
  #overall
  for (m in 1:3){
    current_train_set <- simulated_df[final_train_index_list[[m]],]
    current_result <- logistic_model_train(current_train_set,df.test,number_var )
    current_model_acc <- current_result[[1]]
    current_model_auc <- current_result[[2]]
    
    current_acc <-  rbind(current_acc,current_model_acc)
    current_auc <-  rbind(current_auc,current_model_auc)
  }
  
  colnames(current_acc) = c("overall","cluster 1","cluster 2","cluster 3","cluster 4")
  colnames(current_auc) = c("overall","cluster 1","cluster 2","cluster 3","cluster 4")
  
  rownames(current_acc) = c("random sample","cluster sample","true cluster sample")
  rownames(current_auc) = c("random sample","cluster sample","true cluster sample")
  print(paste("potential k are:",optimal_k,sep = ""))
  
  
  
  return(list(current_acc,current_auc,distance_df,
              random_1_10_distance,optimal_1_10_distance,optimal_true_1_10_distance))
}

#######################################################
# Simulated
#######################################################

df <- read.csv(data_file_name)

true_group = df$id
min_value <- min(df)

df_positive <- df + abs(min_value)+0.5
df_positive$id <- true_group

#set.seed(1234)
# set the effect for logistic model
ahpha_0 <- 1 #intercept
#ahpha_i<-c(-1,-2,2,3)
ahpha_i<-c(-1,-0.5,0.5,1)

#beta_l <- rep(0.5,in.nvar)

beta_l <- c(0.5,0.5,-2,0.5,0.5,0.5,1,-1,-0.5,-0.5)

effect_vector <- as.matrix(c(beta_l,ahpha_i ))

#full_df <- gen_event_interaction(ahpha_0,effect_vector,df,10000 )

full_df <- gen_event(ahpha_0,effect_vector,df_positive,10000 )

table(full_df$event)
table(full_df[full_df$id==1,1])
table(full_df[full_df$id==2,1])
table(full_df[full_df$id==3,1])
table(full_df[full_df$id==4,1])


euclidean_dm <- as.matrix(daisy(full_df[,3:(2+5)] , metric = "euclidean")) # use 5 variable in clusteirng analysis
euclidean_hierarchical_clusters <- fastcluster::hclust(as.dist(euclidean_dm),method = "ward.D2")
set.seed(235)
result=generate_AUC(full_df[,1:(2+5)],euclidean_hierarchical_clusters,euclidean_dm,120,5)
print(result[[2]])
#file_name = paste("result/",number_minority,"_minority_",number_iteration,"_iteration.RData",sep = "")

file_name = paste("result/","replication_",number_iteration,
                  "_cluster_", number_minority,".RData",sep = "")


#save(result,file = file_name)