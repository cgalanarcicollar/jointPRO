# Example 
source("sim_jmCOPD.R")
source("sim_jmPRO.R")
source("PROjm_plot.R")
source("JMBB_model.R")
source("TSBB_model.R")
source("dynpred_jmPRO.R")
source("plot_dynpred.R")
set.seed(8)
n <- 200 # number of subjects
t.max <- 6.5 # maximum follow-up time
nu <- 1.2 # shape for the Weibull baseline hazard
gamma <- c(-3.5) 

cens <- 0.1 # the random censoring

betas <- c("Intercept" = -0.19 , "time" = 0.03) 
ntrial <- 24
D <- matrix(c(1.6,0,0,0.05),ncol=2)

phi <- 0.05
alpha <- 2


# Simulate jm with PRO longitudinal response 

jm_data <- sim_jm.PRO(n = n,  t.max = t.max, minobs = 2, betas = betas, 
                      phi = phi, ntrial = ntrial, alpha = alpha, nu = nu,
                      gamma = gamma, D = D, type = "p")
  
long_data <- jm_data$longitudinal
surv_data <- jm_data$survival
head(jm_data$longitudinal)
head(jm_data$survival)

# Simulate jm data with PRO longitudinal response similar to COPD study

jm_data <- sim_jm.COPD(n = n, t.max = t.max, betas = betas, phi = phi,
                       ntrial = ntrial, alpha = alpha, nu = nu, gamma = gamma,
                       D = D, cens = cens, type = "p")

long_data <- jm_data$longitudinal
surv_data <- jm_data$survival
head(jm_data$longitudinal)
head(jm_data$survival)

# VISUALIZE THE DATA TOO 

plot_longitudinal_survival(data_long = jm_data$longitudinal, outcome_col = "y", time_col = "time", data_surv = jm_data$survival,Time_col = "Time", 
                           status_col = "status", id_col = "id", y_max = ntrial, ids = c(seq(5)), max_time = 6.5)
  
#Fit TSBB model 

fixed.formula <- list ( y ~ time )
random.formula <- list(  ~ -1 + id + time: id)


TSBB_ex <- TSBB(fixed.formula = fixed.formula, random.formula = random.formula,  
     m=ntrial, data_long = long_data, time_col = "time", outcome_col = "y",
     id_col = "id", Time = surv_data$Time, status = surv_data$status, 
     data_surv = surv_data, Time_col = "Time", event_col = "status",nDim = 1, type = "p") 



# Fit JMBB model 

long_data$id <- as.numeric(long_data$id)
JMBB_ex <-  JMBB(data_long = long_data, id_col= "id", outcome_col = "y", m = ntrial,
                 meas_col ="time", Time_col = "Time", status_col = "status", 
                 data_surv = surv_data, nDim = 1, chains = 3, iter = 2000,
                 warmup = 1000 , type ="p") 



# Perform dynamic predictions : 

# First divide training and test samples 

set.seed(35)
test_id <- sample(1:n, 2) #two test id 

test <- long_data[(long_data$id %in% test_id),]
test_1 <- surv_data[(surv_data$id %in% test_id),]
test_1

train <- long_data[!(long_data$id %in% test_id),]
train$idx<-match(train$id, unique(train$id))
train_1 <-  surv_data[!(surv_data$id %in% test_id),]
train_1$idx<-match(train_1$id, unique(train_1$id))


# Fit the model to train data set 

JMBB_train <-  JMBB(data_long = train, id_col= "idx", outcome_col = "y", m=ntrial, meas_col ="time",
                 Time_col = "Time", status_col = "status", data_surv = train_1, nDim = 1,
                 chains = 3, iter = 2000, warmup = 1000 , type ="p") 

# The id data to which we aim to perform the dynamic prediction

Long <- test[which(test$id==test_id[1]),]
surv <- test_1[which(test_1$id==test_id[1]),]

data_id <- list(longitudinal=Long,survival=surv)

dynpred <- dynpred_jmPRO(id_data = data_id, m=ntrial, jm_fit=JMBB_train, future_time= 3)
  

plots <- plot_dynpred(dynpred,id_data = data_id, m = ntrial, future_time = 3 )

plots[[1]]
plots[[2]]
plots[[3]]
