#
# Target-side data analysis. Now includes both human and source covariates.
#
# NOTE: Edit distance and pause features are heavily correlated with
#       time. Don't add those since they screw up the model fitting.
#       Don't add them since we aren't making any assumptions about
#       the generation process.
#

# Libraries
library(ez)
library(languageR)

# PTM stuff
source("ptm_utils.R")


###########################################################
# Arabic data
###########################################################
ar.d <- ptm.targetframe("trans_frame.ar.csv")

# Step 0: Raw data inspection
xylowess.fnc(log(time) ~ src_id | user_id, d=ar.d, ylabel = "Log normalized time")
xylowess.fnc(rank ~ src_id | user_id, d=ar.d, ylabel = "Ranking")
ptm.edatarget(ar.d, "Arabic")
ptm.hcimodels(ar.d)

# Step 1: Tests for normality and correlation of response variables
qqmath(~log(time) | user_id, data=ar.d)
qqmath(~rank | user_id, data=ar.d)


# Step 2: Light data filtering
# TODO: if needed

# ME time models
random.ar.tm = "(1 + ui_id|user_id) + (1 + ui_id|src_id)"
fixed.ar.tm = "gscale(hours_per_week) + gscale(hourly_rate) + en_level + gscale(en_spell) + gscale(en_vocab) + gscale(en_skills) + gscale(en_usage) + gscale(en_ar_trans) + gscale(log(syn_complexity)) + gscale(cnt_entity_tokens) + gscale(log(mean_lexical_freq)) + gscale(cnt_vb) + gscale(cnt_nn) + gscale(cnt_jj) + gscale(cnt_adv) + gscale(log(src_len))"
lmm.ar.tm.0 <- ptm.fittime(fixed.ar.tm,random.ar.tm,ar.d, random.ar.tm)

# Model criticism
ar.d.2 <- ar.d[abs(scale(resid(lmm.ar.tm.0))) < 2.5, ]
1.0 - (nrow(ar.d.2) / nrow(ar.d))
dotchart(sort(xtabs(~ ar.d.2$user_id)), cex=0.7)
dotchart(sort(xtabs(~ ar.d.2$src_id)), cex=0.7)

# Fit final time models
lmm.ar.tm.0 <- ptm.fittime(fixed.ar.tm, random.ar.tm, ar.d.2, random.ar.tm)
print(summary(lmm.ar.tm.0),corr=F)
ptm.lrtime(fixed.ar.tm,random.ar.tm,ar.d.2)

# Rank models (ordinal)
# NOTE: Random effect for src_id not justified (use clmm for more that one random effect, but this method seems very unstable...and very slow).
# Also, this model has convergence issues so I pruned the number of
# covariates
cmm.ar.rn.0 <- clmm(as.ordered(rank) ~ ui_id + gscale(hourly_rate) + gscale(cnt_nn) + gscale(cnt_jj) + gscale(cnt_adv) + (1|user_id), data=ar.d, Hess=T)
summary(cmm.ar.rn.0)
cmm.ar.rn.ref <- clmm(as.ordered(rank) ~ gscale(hourly_rate) + gscale(cnt_nn) + gscale(cnt_jj) + gscale(cnt_adv) + (1|user_id), data=ar.d, Hess=T)
anova(cmm.ar.rn.ref, cmm.ar.rn.0)


##############################################################
# German data
##############################################################

de.d <- ptm.targetframe("trans_frame.de.csv")
xylowess.fnc(log(time) ~ src_id | user_id, d=de.d, ylabel = "Log normalized time")
xylowess.fnc(rank ~ src_id | user_id, d=de.d, ylabel = "Ranking")
ptm.edatarget(de.d, "German")
ptm.hcimodels(de.d)

# Step 1: Tests for normality and correlation of response variables
qqmath(~log(time) | user_id, data=de.d)
qqmath(~rank | user_id, data=de.d)


# ME Time models
random.de.tm = "(1 + ui_id|user_id) + (1 + ui_id|src_id)"
fixed.de.tm = " gscale(hours_per_week) + gscale(hourly_rate) + gscale(en_spell) + gscale(en_vocab) + gscale(en_skills) + gscale(en_usage) + gscale(de_spell) + gscale(de_vocab) + gscale(en_de_trans) + gscale(log(syn_complexity)) + gscale(cnt_entity_tokens) + gscale(log(mean_lexical_freq)) + gscale(cnt_vb) + gscale(cnt_nn) + gscale(cnt_jj) + gscale(cnt_adv) + gscale(log(src_len))"
lmm.de.tm.0 <- ptm.fittime(fixed.de.tm,random.de.tm,de.d, random.de.tm)

# Model criticism (for time models)
de.d.2 <- de.d[abs(scale(resid(lmm.de.tm.0))) < 2.5, ]
1.0 - (nrow(de.d.2) / nrow(de.d))
dotchart(sort(xtabs(~ de.d.2$user_id)), cex=0.7)
dotchart(sort(xtabs(~ de.d.2$src_id)), cex=0.7)

# Fit new time models
lmm.de.tm.0 <- ptm.fittime(fixed.de.tm, random.de.tm, de.d.2, random.de.tm)
print(summary(lmm.de.tm.0),corr=F)
ptm.lrtime(fixed.de.tm,random.de.tm,de.d.2)

# ME Rank models
cmm.de.rn.0 <- clmm(as.ordered(rank) ~ ui_id + gscale(en_vocab) + gscale(en_skills) + gscale(de_spell) + gscale(de_vocab) + (1|user_id), data=de.d, Hess=T)
summary(cmm.de.rn.0)
cmm.de.rn.ref <- clmm(as.ordered(rank) ~ 1 + gscale(en_vocab) + gscale(en_skills) + gscale(de_spell) + gscale(de_vocab) + (1|user_id), data=de.d, Hess=T)
anova(cmm.de.rn.ref, cmm.de.rn.0)


######################################################################
# French data
######################################################################
fr.d <- ptm.targetframe("trans_frame.fr.csv")
xylowess.fnc(log(time) ~ src_id | user_id, d=fr.d, ylabel = "Log normalized time")
xylowess.fnc(rank ~ src_id | user_id, d=fr.d, ylabel = "Ranking")
ptm.edatarget(fr.d, "French")
ptm.hcimodels(fr.d)

# Step 1: Tests for normality and correlation of response variables
qqmath(~log(time) | user_id, data=fr.d)
qqmath(~rank | user_id, data=fr.d)


# ME Time models
random.fr.tm = "(1 + ui_id|user_id) + (1 + ui_id|src_id)"
fixed.fr.tm = "gscale(hours_per_week) + gscale(hourly_rate) + en_level + gscale(en_spell) + gscale(en_vocab) + gscale(en_skills) + gscale(en_usage) + gscale(fr_spell) + gscale(fr_vocab) + gscale(fr_usage) + gscale(en_fr_trans) + gscale(log(syn_complexity)) + gscale(cnt_entity_tokens) + gscale(log(mean_lexical_freq)) + gscale(cnt_vb) + gscale(cnt_nn) + gscale(cnt_jj) + gscale(cnt_adv) + gscale(log(src_len))"
lmm.fr.tm.0 <- ptm.fittime(fixed.fr.tm,random.fr.tm,fr.d,random.fr.tm)

# Model criticism (time models)
fr.d.2 <- fr.d[abs(scale(resid(lmm.fr.tm.0))) < 2.5, ]
1.0 - (nrow(fr.d.2) / nrow(fr.d))
dotchart(sort(xtabs(~ fr.d.2$user_id)), cex=0.7)
dotchart(sort(xtabs(~ fr.d.2$src_id)), cex=0.7)

# Fit new time models
lmm.fr.tm.0 <- ptm.fittime(fixed.fr.tm, random.fr.tm, fr.d.2,random.fr.tm)
print(summary(lmm.fr.tm.0),corr=F)
ptm.lrtime(fixed.fr.tm,random.fr.tm,fr.d.2)

# ME ranking models
cmm.fr.rn.0 <- clmm(as.ordered(rank) ~ ui_id + gscale(en_skills) + gscale(en_usage) + gscale(en_spell) + gscale(fr_vocab) + (1|user_id), data=fr.d, Hess=T)
summary(cmm.fr.rn.0)
cmm.fr.rn.ref <- clmm(as.ordered(rank) ~ 1 + gscale(en_skills) + gscale(en_usage) + gscale(en_spell) + gscale(fr_vocab) + (1|user_id), data=fr.d, Hess=T)
anova(cmm.fr.rn.ref, cmm.fr.rn.0)


################################################################
# All languages -- language modeled as a random effect
################################################################
ptm.d <- ptm.allframes()

#
# The baseline HCI models are so clearly wrong that we don't even
# report them. Here they are for completeness.
#
# Baseline time model: 3x3 ANOVA
#summary(aov(log(time) ~ ui_id*src_id*tgt_lang, data=ptm.d))

# Baseline rank ordinal model
#summary(clm(as.ordered(rank) ~ ui_id*src_id*tgt_lang, data=ptm.d, Hess=T))

# ME time model
random.all <- "(1 + ui_id|user_id) + (1 + ui_id|src_id) + (1 + ui_id|tgt_lang)"
fixed.all <- "gscale(hours_per_week) + gscale(hourly_rate) + en_level + gscale(en_spell) + gscale(en_vocab) + gscale(en_skills) + gscale(en_usage) + gscale(log(syn_complexity)) + gscale(cnt_entity_tokens) + gscale(log(mean_lexical_freq)) + gscale(cnt_vb) + gscale(cnt_nn) + gscale(cnt_jj) + gscale(cnt_adv) + gscale(log(src_len))"
lmm.all.0 <- ptm.fittime(fixed.all,random.all,ptm.d,random.all)

# Model criticism
ptm.d.2 <- ptm.d[abs(scale(resid(lmm.all.0))) < 2.5, ]
lmm.all.0 <- ptm.fittime(fixed.all,random.all,ptm.d.2,random.all)
print(summary(lmm.all.0),corr=F)
ptm.lrtime(fixed.all,random.all,ptm.d.2)

# ME rank model
# src_id does not change logLik, so remove it.
# Also, can't do the LR test if we add additional numeric covariates
# It would be great to add some additional covariates, but they don't
# seem to work.
cmm.all.rn.0 <- clmm(as.ordered(rank) ~ ui_id + (1|user_id) + (1|tgt_lang), data=ptm.d, Hess=T)
summary(cmm.all.rn.0)
cmm.all.rn.ref <- clmm(as.ordered(rank) ~ 1  + (1|user_id) + (1|tgt_lang), data=ptm.d, Hess=T)
anova(cmm.all.rn.ref, cmm.all.rn.0)

# Pause models
# Hack: one data point has a pause_ratio1s value of 0, which prevents
# convergence
ptm.d.2 <- ptm.d[ptm.d$pause_ratio1s > 0,]
ptm.fitpause(fixed.all,random.all,ptm.d.2,random.all)
