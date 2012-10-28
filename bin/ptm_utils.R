#
# Misc. utils for working with ptm data
#

# Libraries
library(lme4)
library(ordinal)

# Convenience libraries
library(ez)
library(languageR)

#
# complement of the %in% operator
#
`%ni%` <- Negate(`%in%`)


# x := name of the target translation frame
# returns: an R data frame
#
ptm.targetframe <- function(x){
  # Load the source frame
  d.src <- read.csv("source_frame.csv",header=T)

  # Load the target frame
  d.tgt <- read.csv(x,header=T)

  # User data
  d.user <- read.csv("user_frame.csv",header=T)

  # Merge these three frames
  d.tmp1 <- merge(d.src,d.tgt,by.x = "src_id", by.y = "src_id")
  my.data <- merge(d.tmp1,d.user,by.x = "user_id", by.y = "user_id")

  # Set factor levels
  my.data$src_id <- as.factor(my.data$src_id)
  levels(my.data$src_id) <- paste("src",levels(my.data$src_id),sep="")
  my.data$user_id <- as.factor(my.data$user_id)
  levels(my.data$user_id) <- paste("usr",levels(my.data$user_id),sep="")
  my.data$ui_id <- as.factor(my.data$ui_id)
  levels(my.data$ui_id) <- paste("ui",levels(my.data$ui_id),sep="")
  my.data$en_level <- as.factor(my.data$en_level)
  levels(my.data$en_level) <- paste("en",levels(my.data$en_level),sep="")
  
  # Change the raw token counts to proportions
  my.data$cnt_entity_tokens <-  my.data$cnt_entity_tokens / my.data$src_len
  my.data$cnt_nn <- my.data$cnt_nn / my.data$src_len
  my.data$cnt_jj <- my.data$cnt_jj / my.data$src_len
  my.data$cnt_adv <- my.data$cnt_adv / my.data$src_len
  my.data$cnt_vb <- my.data$cnt_vb / my.data$src_len
  
  return(my.data)
}

ptm.allframes <- function(){
  ar.d <- ptm.targetframe("trans_frame.ar.csv")
  de.d <- ptm.targetframe("trans_frame.de.csv")
  fr.d <- ptm.targetframe("trans_frame.fr.csv")
  my.data <- rbind(ar.d,de.d,fr.d)
  return(my.data)
}


#
# Exploratory data analysis on a target frame
#
ptm.edatarget <- function(data, lang){
  attach(data)
  # Diagnostics -- time
  hist(data$time, main=paste(lang,"time"))
  hist(log(data$time), main=paste(lang,"log time"))
  boxplot(time ~ ui_id, xlab="UI condition",ylab="time (ms)",main=paste(lang,"time by ui id"))
  boxplot(time ~ user_id, xlab="UI condition",ylab="time (ms)",main=paste(lang,"time by user id"))

# Diagnostics -- ranking
  boxplot(rank ~ ui_id, xlab="UI condition",ylab="Ranking",main=paste(lang,"ranking by ui id"))
  boxplot(rank ~ user_id, xlab="UI condition",ylab="Ranking",main=paste(lang,"ranking by user id"))

  # Correlation tests
  print(cor.test(data$rank, data$time))
  print(cor.test(data$rank, log(data$time)))

  detach(data)
}

#
# Gelman standardization 
#
# x := a vector to standardize
# returns: standardized vector
#
gscale <- function(x){
  y <- scale(x, center=T, scale=2*sd(x))
  return(y)
}


#
# Chi-square test
#
# x := chiquare value
# df := degrees of freedom
#
chisqtest <- function(x,df){
  return(1-pchisq(x,df))
}


#
# LINEAR MIXED EFFECTS MODELS
#
# Significance testing tips (from Barr et al., 2012):
#
#
# 1) Use likelihood ratio test vv. a model with the fixed effect
#    of interest removed.
# 2) Do *not* remove the random effects associated with the fixed effect
#    from the reference model.
#
#


#
# Response variable: log(time)
#
ptm.fittime <- function(fixed,random,data,random_ref="(1|user_id) + (1|src_id)"){
  tm.form.0 <- as.formula(paste("log(time) ~ ui_id + ",random,"+",fixed))
  lmm.tm.0 <- lmer(tm.form.0, data = data, REML=F)
  #print(summary(lmm.tm.0))

  tm.ref <- as.formula(paste("log(time) ~ ",random_ref,"+",fixed))
  lmm.tm.ref <- lmer(tm.ref, data = data, REML=F)

  print(anova(lmm.tm.ref, lmm.tm.0))
  plot(fitted(lmm.tm.0), resid(lmm.tm.0), main="(log) time",xlab="Fitted values",ylab="Residuals")

  ## Model fitting
  print("Coefficients of determination (ref then model)")
  print(cor(fitted(lmm.tm.ref), log(data$time))^2)
  print(cor(fitted(lmm.tm.0), log(data$time))^2)

  return(lmm.tm.0)
}

#
# Do per factor LR testing of each covariate in the fixed effects structures
#
ptm.lrtime <- function(fixed,random,data) {
  factors <- strsplit(fixed,"+",fixed=T)[[1]]
  factors <- c("ui_id",factors)
  tm.form.0 <- as.formula(paste("log(time) ~ ui_id + ",random,"+",fixed))
  lmm.tm.0 <- lmer(tm.form.0, data = data, REML=F)
  for (i in 1:length(factors)){
    ref.factors <- c(factors[1:(i-1)],factors[(i+1):length(factors)])
    if (i == 1) {
      ref.factors <- factors[2:length(factors)]
    } else if (i == length(factors)) {
      ref.factors <- factors[1:(i-1)]
    }

    print(paste("TEST FACTOR: ",factors[i]))
    
    ref.str <- paste0(ref.factors,collapse="+")
    tm.form.ref <- as.formula(paste("log(time) ~ ",random,"+",ref.str))
    lmm.tm.ref <- lmer(tm.form.ref, data = data, REML=F)
    print(anova(lmm.tm.ref, lmm.tm.0))
  }

}

#
# Response variable: 
#
ptm.fitpause <- function(fixed,random,data,random_ref="(1|user_id) + (1|src_id)"){
  print("Pause count (300ms)")
  form.0 <- as.formula(paste("log(pause_cnt) ~ ui_id + ",random,"+",fixed))
  lmm.0 <- lmer(form.0, data = data, REML=F)
  print(summary(lmm.0), corr=F)

  form.ref <- as.formula(paste("log(pause_cnt) ~ ",random_ref,"+",fixed))
  lmm.ref <- lmer(form.ref, data = data, REML=F)

  print(anova(lmm.ref, lmm.0))

  print("Pause count (1000ms)")
  form.0 <- as.formula(paste("log(pause_cnt1s) ~ ui_id + ",random,"+",fixed))
  lmm.0 <- lmer(form.0, data = data, REML=F)
  print(summary(lmm.0), corr=F)

  form.ref <- as.formula(paste("log(pause_cnt1s) ~ ",random_ref,"+",fixed))
  lmm.ref <- lmer(form.ref, data = data, REML=F)

  print(anova(lmm.ref, lmm.0))
  
  print("Pause mean duration (300ms cutoff)")
  form.0 <- as.formula(paste("log(pause_mean) ~ ui_id + ",random,"+",fixed))
  lmm.0 <- lmer(form.0, data = data, REML=F)
  print(summary(lmm.0), corr=F)

  form.ref <- as.formula(paste("log(pause_mean) ~ ",random_ref,"+",fixed))
  lmm.ref <- lmer(form.ref, data = data, REML=F)

  print(anova(lmm.ref, lmm.0))

  print("Pause ratio (300ms)")
  form.0 <- as.formula(paste("log(pause_ratio) ~ ui_id + ",random,"+",fixed))
  lmm.0 <- lmer(form.0, data = data, REML=F)
  print(summary(lmm.0), corr=F)

  form.ref <- as.formula(paste("log(pause_ratio) ~ ",random_ref,"+",fixed))
  lmm.ref <- lmer(form.ref, data = data, REML=F)

  print(anova(lmm.ref, lmm.0))

  print("Pause ratio (1000ms)")
  form.0 <- as.formula(paste("log(pause_ratio1s) ~ ui_id + ",random,"+",fixed))
  lmm.0 <- lmer(form.0, data = data, REML=F)
  print(summary(lmm.0), corr=F)

  form.ref <- as.formula(paste("log(pause_ratio1s) ~ ",random_ref,"+",fixed))
  lmm.ref <- lmer(form.ref, data = data, REML=F)

  print(anova(lmm.ref, lmm.0))
}


ptm.hcimodels <- function(data){
  # Standard HCI models
  # WARNING: These models are subject to language-as-fixed-effect fallacy,
  # which is evident by the inclusion of src_id as an independent variable.
  
  # Model 0: Simple one-way layout
  # Univariate, mixed-design models
  #
  # This is not an RM-ANOVA (within-subjects) because every subject did not
  # see every treatment level.
  print("RESPONSE: time (2x2 ANOVA with interactions)")
  lm.0 <- aov(log(time) ~ ui_id*src_id, data=data)
  print(summary(lm.0))

  print("RESPONSE: Rank (2x2 ordinal regression)")
  print(summary(clm(as.ordered(rank) ~ ui_id*src_id, data=data, Hess=T)))

  # Multivariate models
  # Note: treating rank as a numeric response here since we can't specify
  # a different family for each response.
  print("RM-MANOVA: time and rank")
  lm.3 <- manova(cbind(log(time),rank) ~ ui_id*src_id, data=data)
  print(summary(lm.3, test="Hotelling-Lawley"))
  print(summary(lm.3, test="Wilks"))
}
