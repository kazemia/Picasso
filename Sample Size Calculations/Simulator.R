library(dplyr)
library(tidyr)
library(purrr)
library(doParallel)
library(ggplot2)
library(ggpubr)
library(latex2exp)

#Determine a directory to save data
setwd("C:/Users/kazem/Dropbox/PICASSO/Statistics/Sample Size Calculations")

#OBS: The following chunck of code generates the simulated data sets and 
#takes some time to run. Run Only if new simulations are needed. Otherwise,
#Scroll down to read the already simulated data.

data_list <- list()
#Number of simulations pr hypothesis configuration
n_sim = 100
#Clinically significant difference
r_diff = 12
#Standard deviation
VAS_sd = 25

p_value_calculator <- function(sub_data){
  #This function takes in a data set and returns p-values the hypothesis tests,
  #described in the protokoll.
  
  #Compare the outcome at 4 weeks between the stroid and the placebo arm with a
  #linear regression model and extract the p-value.
  p_w4 <- sub_data %>%
    filter(WT %in% c("VAS_w4P", "VAS_w4S")) %>%
    lm(VAS ~ WT, data = .) %>%
    summary() %>%
    pluck("coefficients",8)
  
  #Do a global test at 12 weeks among all arms with a
  #linear regression model and anova. Extract the p-value.
  p_w12G <- sub_data %>%
    filter(WT %in% c("VAS_w12P", "VAS_w12S", "VAS_w12T")) %>%
    lm(VAS ~ WT, data = .) %>%
    anova() %>%
    pluck("Pr(>F)", 1)
  
  #Compare the outcome at 12 weeks between the steroid and the therapy arm with a
  #linear regression model and extract the p-value.
  p_w12ST <- sub_data %>%
    filter(WT %in% c("VAS_w12S", "VAS_w12T")) %>%
    lm(VAS ~ WT, data = .) %>%
    summary() %>%
    pluck("coefficients", 8)
  
  #Compare the outcome at 12 weeks between the placebo and the therapy arm with a
  #linear regression model and extract the p-value.
  p_w12TP <- sub_data %>%
    filter(WT %in% c("VAS_w12P", "VAS_w12T")) %>%
    lm(VAS ~ WT, data = .) %>%
    summary() %>%
    pluck("coefficients", 8)
  
  #Compare the outcome at 12 weeks between the steroid and the placebo arm with a
  #linear regression model and extract the p-value.
  p_w12PS <- sub_data %>%
    filter(WT %in% c("VAS_w12S", "VAS_w12P")) %>%
    lm(VAS ~ WT, data = .) %>%
    summary() %>%
    pluck("coefficients", 8)
  
  data.frame(p_w4, p_w12G, p_w12ST, p_w12TP, p_w12PS)
}
for (counter in 1:7) {
  #Determine the number of patients
  n_p = (counter + 7)*30
  #Create a grid for all possible values of mean outcome
  my_data <- expand.grid(VAS_w4P = rep(NA, n_p/3),
                         mu_w4P = c(0,1,2)*r_diff,
                         mu_w4S = c(0,1,2)*r_diff,
                         mu_w12P = c(0,1,2)*r_diff,
                         mu_w12S = c(0,1,2)*r_diff,
                         mu_w12T = c(0,1,2)*r_diff,
                         sim = 1:n_sim)
  n <- nrow(my_data)
  #Determine which alternative hypotheses are true under each configuration
  #Generate data from a normal distribution
  my_data <- my_data %>% mutate(H_w4 = (mu_w4P != mu_w4S),
                                H_w12ST = (mu_w12S != mu_w12T),
                                H_w12TP = (mu_w12P != mu_w12T),
                                H_w12PS = (mu_w12P != mu_w12S),
                                VAS_w4P = rnorm(n, mu_w4P, VAS_sd), 
                                VAS_w12P = rnorm(n, mu_w12P, VAS_sd),
                                VAS_w4S = rnorm(n, mu_w4S, VAS_sd),
                                VAS_w12S = rnorm(n, mu_w12S, VAS_sd),
                                VAS_w12T = rnorm(n, mu_w12T, VAS_sd))
  #Generate a unique ID for each experiment
  my_data <- my_data %>% 
    mutate(experiment = as.numeric(as.factor(paste(
      mu_w4P, mu_w4S, mu_w12P, mu_w12S, mu_w12T, sim))))
  #Transform the data to a long format
  my_data_long <- my_data %>% 
    pivot_longer(cols = colnames(my_data)[grepl("VAS", colnames(my_data))],
                 names_to = "WT", values_to = "VAS")
  #Remove the actual data from the output
  my_data <- my_data %>% 
    select(-all_of(colnames(my_data)[grepl("VAS", colnames(my_data))])) %>% 
    unique()
  my_exps <- unique(my_data$experiment)
  
  #Calculate p-values for all experiments
  myCluster <- makeCluster(15)
  registerDoParallel(myCluster)
  my_pvalues <- foreach(experi = my_exps, .combine = rbind,
          .packages = c("dplyr", "purrr")) %dopar% {
            p_value_calculator(my_data_long %>% filter(experiment == experi))
            }
  stopCluster(myCluster)
  my_pvalues$experiment <- my_exps
  remove(my_data_long)
  #Add the p-values to the output
  my_data <- merge(my_data, my_pvalues)
  remove(my_pvalues)
  #Save the results to memory and disk
  data_list[[counter]] <- my_data
  save(my_data, file = paste0("sim_data", counter, ".Rdata"))
}


#OBS: The following chunck of code reads the already simulated data sets and 
#produces power and error calculations.
data_list <- list()
for (counter in 1:7) {
  load(paste0("sim_data", counter, ".Rdata"))
  data_list[[counter]] <- my_data
}

#Determine which null hypotheses are rejected and
#where a type I or a type II error was committed.
data_list <- lapply(data_list, function(x){
  x %>% mutate(GoEasy = (p_w4 < 0.025),
               H_w4_r = (p_w4 < 0.025),
               H_w12G_r = (GoEasy & (p_w12G < 0.05)) | ((p_w12G < 0.025)),
               H_w12ST_r = (GoEasy & H_w12G_r & (p_w12ST < 0.05)) | 
                 (H_w12G_r & (p_w12ST < 0.025)),
               H_w12TP_r = (GoEasy & H_w12G_r & (p_w12TP < 0.05)) | 
                 (H_w12G_r & (p_w12TP < 0.025)),
               H_w12PS_r = (GoEasy & H_w12G_r & (p_w12PS < 0.05)) | 
                 (H_w12G_r & (p_w12PS < 0.025)),
               H_w12G = H_w12PS | H_w12ST | H_w12TP,
               H_w4_t1 = H_w4_r & (!H_w4),
               H_w4_t2 = H_w4 & (!H_w4_r),
               H_w12G_t1 = H_w12G_r & (!H_w12G),
               H_w12G_t2 = H_w12G & (!H_w12G_r),
               H_w12ST_t1 = H_w12ST_r & (!H_w12ST),
               H_w12ST_t2 = H_w12ST & (!H_w12ST_r),
               H_w12TP_t1 = H_w12TP_r & (!H_w12TP),
               H_w12TP_t2 = H_w12TP & (!H_w12TP_r),
               H_w12PS_t1 = H_w12PS_r & (!H_w12PS),
               H_w12PS_t2 = H_w12PS & (!H_w12PS_r),
               fwt1 = (H_w4_t1 | H_w12G_t1 | H_w12ST_t1 | H_w12TP_t1 | H_w12PS_t1),
               FDR = (H_w4_t1 + H_w12G_t1 + H_w12ST_t1 + H_w12TP_t1 + H_w12PS_t1)/5)
})

#Calculate the FWER, the FDR and the marginal and conjunctive power
error_df <- data.frame(n_p = (1:7 + 7)*30,
                   fwer_weak = sapply(data_list,
                                      function(x){
                                        (x %>% filter((!H_w4) & (!H_w12G)) %>% 
                                           summarise_all(mean))$fwt1
                                        }),
                   fwer_strong = sapply(data_list,
                                      function(x){
                                        (x %>% summarise_all(mean))$fwt1
                                      }),
                   fdr_weak = sapply(data_list,
                                      function(x){
                                        (x %>% filter((!H_w4) & (!H_w12G)) %>% 
                                           summarise_all(mean))$FDR
                                      }),
                   fdr_strong = sapply(data_list,
                                        function(x){
                                          (x %>% summarise_all(mean))$FDR
                                        }),
                   H_w4_power = 1 - sapply(data_list,
                                           function(x){
                                             (x %>% filter(H_w4) %>% 
                                                summarise_all(mean))$H_w4_t2
                                             }),
                   H_w12G_power = 1 - sapply(data_list,
                                             function(x){
                                               (x %>% filter(H_w12G) %>% 
                                                  summarise_all(mean))$H_w12G_t2
                                               }),
                   H_w12ST_power = 1 - sapply(data_list,
                                              function(x){
                                                (x %>% filter(H_w12ST) %>% 
                                                   summarise_all(mean))$H_w12ST_t2
                                                }),
                   
                   H_w12TP_power = 1 - sapply(data_list, 
                                             function(x){
                                               (x %>% filter(H_w12TP) %>% 
                                                  summarise_all(mean))$H_w12TP_t2
                                               }),
                   
                   H_w12PS_power = 1 - sapply(data_list,
                                          function(x){
                                            (x %>% filter(H_w12PS) %>%
                                               summarise_all(mean))$H_w12PS_t2
                                            }),
                   con_power = sapply(data_list,
                                      function(x){
                                        (x %>% filter(H_w4 & H_w12ST & H_w12TP &
                                                        H_w12PS) %>% 
                                           mutate(con_t2 = H_w4_r & H_w12G_r & 
                                                    H_w12ST_r & H_w12TP_r & 
                                                    H_w12PS_r) %>% 
                                           summarise_all(mean))$con_t2
                                            })
                   )

power_df <- error_df %>% 
  pivot_longer(cols = colnames(error_df)[grepl("power", colnames(error_df))],
               names_to = "hypothesis", values_to = "power")

#Plot the results
ggarrange(plotlist = list(ggplot() + 
                            geom_line(data = power_df %>% 
                                        filter(hypothesis != "con_power"), 
                                      aes(x = n_p, y = power, 
                                          color = hypothesis)) + 
                            geom_vline(xintercept = 319) + 
                            scale_x_continuous(n.breaks = 6) + 
                            scale_y_continuous(n.breaks = 6) + 
                            ylab("Marginal power") + 
                            xlab("Number of patients") +
                            scale_color_discrete(name = "Hypothesis:", 
                                                 labels = c(TeX(r'($H^{W12G}$)'), 
                                                            TeX(r'($H^{W12PS}$)'), 
                                                            TeX(r'($H^{W12ST}$)'),
                                                            TeX(r'($H^{W12TP}$)'),
                                                            TeX(r'($H^{W4}$)'))),
                          ggplot() + 
                            geom_line(data = power_df %>% 
                                        filter(hypothesis == "con_power"), 
                                      aes(x = n_p, y = power)) + 
                            geom_vline(xintercept = 319) + 
                            scale_x_continuous(n.breaks = 6) + 
                            scale_y_continuous(n.breaks = 6) + 
                            ylab("Conjunctive power") +
                            xlab("Number of patients")), nrow = 1, ncol = 2)

ggplot() + geom_line(data = error_df, 
                     aes(x = n_p, y = fwer_weak)) + 
  geom_hline(yintercept = 0.05) + 
  ylab(TeX(r'(FWER $| \ H_0^{W12G} \cap H_0^{W4}$)')) +
  xlab("Number of patients")

ggarrange(plotlist = list(ggplot() + 
                            geom_line(data = error_df, 
                                      aes(x = n_p, y = fwer_strong)) + 
                            ylab("FWER") +
                            xlab("Number of patients"),
                          ggplot() + 
                            geom_line(data = error_df, 
                                      aes(x = n_p, y = fdr_strong)) + 
                            ylab("FDR") +
                            xlab("Number of patients")), nrow = 1, ncol = 2)



