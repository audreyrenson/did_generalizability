# libs --------------------------------------------------------------------

library(tidyverse)
library(dtplyr)
library(nlme)
library(geepack)
library(splines)
library(furrr)

plan(multisession)
set.seed(seed <- 7437437)

# funs --------------------------------------------------------------------

wmean = function(x, w,...) {
  #simple weighted mean
  w = w/mean(w)
  sum(w*x,...)/sum(w,...)
}

# data --------------------------------------------------------------------

#Air quality data from NYC SFH study
sfh = read_csv('data/sfh_air_quality.csv')

#smoking data from baseline NYC SFH study basline survey and from AHS (joined: s=1 NYC SFH, s=0 AHS)
ahs_sfh = read_csv('data/ahs_sfh_join.csv')

df_trans = sfh %>%
  filter(unit_type2 == 2, TimeEvent %in% c(0, 6)) %>%
  mutate(time = 1*(TimeEvent==6),
         treat = 1-(arm-1)) %>%
  group_by(treat,building_c, time) %>%
  summarise(log_nic = mean(log(nicotine), na.rm=TRUE), .groups='drop') %>%
  group_by(treat, building_c) %>%
  summarise(dY = diff(log_nic), .groups = 'drop') %>%
  right_join(ahs_sfh) %>%
  filter(!is.na(SMKR_bin)) %>%
  mutate(I00 = (1-treat)*(1-s),
         I10 = treat*(1-s),
         I01 = (1-treat)*s,
         I11 = treat*s,
         dY = replace_na(dY, 0))


# exploratory analysis ------------------------------------

df_s1 = df_trans %>% filter(s==1)

#almost replicating Anastasiou (building-level analysis)
lm(dY ~ treat, data=df_s1 %>% select(treat, building_c, dY) %>% unique())

#SATT (person-level analysis)
satt = coef(lm(dY ~ treat, data=df_s1))[2]

#W-SATT (person-level)
lm(dY ~ treat*SMKR_bin, data=df_s1)

#SATT adjusted for W
g_satt = 1 - predict(glm(treat ~ SMKR_bin, family=binomial(), data=df_s1), newdata=df_s1)
m0_satt = predict(lm(dY ~ treat*SMKR_bin, data=df_s1), newdata=df_s1 %>% filter(treat==1) %>% mutate(treat=0))
with(df_s1, mean(treat * dY / g_satt, na.rm=TRUE)) #IOW
with(df_s1[df_s1$treat==1, ], mean(dY - m0_satt, na.rm=TRUE)) #g-comp

#weighted W distribution
#weighted and unweighted W dist in treated
mw_samp = mean(df_adj$SMKR_bin[df_adj$treat==1])
mw_surv = mean(ahs_sfh$SMKR_bin[ahs_sfh$treat==1], na.rm=TRUE)
mw_pop  = wmean(ahs_sfh$SMKR_bin[ahs_sfh$treat==1], w=ahs_sfh$SP2WEIGHT[ahs_sfh$treat==1], na.rm=TRUE)


# transport ----------------------------------------------------------------

est_g = function(df) {
  mod_a = lm(treat ~ SMKR_bin, data=df, weights=SP2WEIGHT)
  mod_s = lm(s ~ SMKR_bin, data=df, weights=SP2WEIGHT)

  list(g10 = predict(mod_a, newdata=df %>% mutate(s==0)) * (1-predict(mod_s)),
       g01 = (1-predict(mod_a, newdata=df %>% mutate(s==1))) * predict(mod_s),
       g11 = predict(mod_a, newdata=df %>% mutate(s==1)) * predict(mod_s)) %>%
    bind_cols()
}
est_m = function(df) {
  #no weights in the outcome model because it only uses data from the sample, so weights=1
  list(m0 = 0,
       m1 = 1) %>%
    map(function(a) predict(lm(dY ~ SMKR_bin*treat, data=df %>% filter(s==1)), newdata=df %>% mutate(treat=a))) %>%
    bind_cols()
}

df_trans_mod = df_trans %>% bind_cols(est_g(.), est_m(.))


#normalizing constant used in all estimators
p10_w = function(df) wmean(df$I10, df$SP2WEIGHT)

#IOW
iow = with(df_trans_mod, wmean( dY * (I11*g10/g11 - I01*g10/g01), SP2WEIGHT)/p10_w(df_trans_mod))

#gcomp
gc = with(df_trans_mod, wmean(I10*(m1-m0), SP2WEIGHT)/p10_w(df_trans_mod))

#dr
dr = with(df_trans_mod, wmean( (dY-m1)*(I11*g10/g11) - (dY-m0)*(I01*g10/g01) + I10*(m1-m0), SP2WEIGHT)/p10_w(df_trans_mod))

# jackknife + bootstrap variance
df_repwgt_wt = df_trans %>% select(starts_with('SP2REP'))
df_boot_wt = tibble(w = as.list(df_repwgt_wt),
                    df_trans = map(w, function(W) mutate(df_trans, SP2WEIGHT=W))) %>%
  mutate(g = future_map(df_trans, est_g),
         m = future_map(df_trans, est_m),
         df_trans = pmap(list(df_trans, g, m), bind_cols),
         satt = map_dbl(df_trans, ~coef(lm(dY ~ treat, data=filter(., s==1), weights=SP2WEIGHT))[2]),
         gcomp = map_dbl(df_trans, ~with(., wmean(I10*(m1-m0), SP2WEIGHT)/p10_w(.))),
         iow   = map_dbl(df_trans, ~with(., wmean( dY * (I11*g10/g11 - I01*g10/g01), SP2WEIGHT)/p10_w(.))),
         dr    = map_dbl(df_trans, ~with(., wmean( (dY-m1)*(I11*g10/g11) - (dY-m0)*(I01*g10/g01) + I10*(m1-m0), SP2WEIGHT)/p10_w(.))))
df_se_wt = df_boot_wt %>%
  select(satt, gcomp, iow, dr) %>%
  summarise(across(everything(), sd)) %>%
  pivot_longer(everything(), names_to = 'par', values_to='se')

df_results_wt = tibble(par = c('satt','gcomp','iow','dr'),
                       est = c(satt, gc, iow, dr)) %>%
  left_join(df_se_wt) %>%
  mutate(lwr = est - 1.96*se, upr = est + 1.96*se)


