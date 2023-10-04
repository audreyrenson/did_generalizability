library(tidyverse)
library(simstudy)
library(furrr)


set.seed(seed <-  5256397)

plan(multisession)

# Data generating model --------------

dgm = defData(varname = 's', formula = 0, dist='binary', link='logit') %>%
  defData(varname = 'u', formula = '-1 + s', dist='binary', link='logit') %>%
  defData(varname = 'w', formula = '0.5 - 0.25*s', dist='binary', link='identity') %>%
  defData(varname = 'a', formula = '0.3 + 0.1*s + 0.1*w + 0.1*u', dist='binary', link='identity') %>%
  defData(varname = 'y0', formula = '1 + w + u', variance = .01) %>%
  defData(varname = 'y1', formula = '0.5*w + u + a + 0.5*w*a', variance=.1) %>%
  defData(varname = 'y1_0', formula = '0.5*w + u + 0 + 0.5*w*0', variance=.1) %>%
  defData(varname = 'y1_1', formula = '0.5*w + u + 1 + 0.5*w*1', variance=.1)


get_true_patt = function(def, n=1e6, seed=1111111, se=FALSE) {
  ncores =  parallel::detectCores()
  n_per_core = n / ncores

  df0 =future_map(2:ncores, ~genData(n=n_per_core, def), .options = furrr_options(seed=seed)) %>%
    bind_rows()

  I10 = df0$a==1 & df0$s==0
  true_patt = with(df0[I10,], mean(y1-y1_0))

  if(!se) return(true_patt) else tibble(true_patt = true_patt, true_patt_se = with(df0[I10,], sd(y1-y1_0)/sum(I10)))
}

# estimators --------------------------------------------------------------


dr = function(data, model) {

  stopifnot(model %in% c('true','gfal','qfal','bfal'))

  formulas = list(
    m0 = if(model %in% c('true','gfal')) y1-y0~w*a*s else y1-y0 ~ w+a,
    m1 = if(model %in% c('true','gfal')) y1-y0~w*a*s else y1-y0 ~ w+a,
    ga = if(model %in% c('true','qfal')) a~ w*s else a~w,
    gs = if(model %in% c('true','qfal')) s~w else s~1
  )

  models = map(formulas, lm, data=data)

  m0 = predict(models$m0, newdata=data %>% mutate(a=0,s=1))
  m1 = predict(models$m1, newdata=data %>% mutate(a=1,s=1))
  g11 = predict(models$ga, newdata=data %>% mutate(s=1)) * predict(models$gs, newdata=data)
  g01 = (1-predict(models$ga, newdata=data %>% mutate(s=1)))*predict(models$gs, newdata=data)
  g10 = predict(models$ga, newdata=data %>% mutate(s=0)) * (1-predict(models$gs, newdata=data))

  p10 = mean(     data$a*(1-data$s)  )

  I11 = 1*( data$a==1 & data$s==1  )
  I01 = 1*( data$a==0 & data$s==1  )
  I10 = 1*( data$a==1 & data$s==0  )

  dy = data$y1 - data$y0

  eif =  I11*g10*(dy - m1)/(p10*g11) -
    I01*g10*(dy - m0)/(p10*g01) +
    I10*(m1 - m0)/p10

  tibble(dr=mean(eif),
         gcomp = mean( I10*(m1 - m0)/p10 ),
         ipw = mean(  I11*g10*dy/(p10*g11) -
                        I01*g10*dy/(p10*g01))
  )
}


# sim ---------------------------------------------------------------------

nsims=200
nobs=1e4

true_patt = get_true_patt(dgm, seed=seed)

df_sims = tibble(sim = 1:nsims,
                 data = future_map(sim, ~genData(nobs, dgm), .options = furrr_options(seed=seed))) %>%
  expand_grid(model = c('true','gfal','qfal','bfal')) %>%
  mutate(est = future_map2(data, model, dr)) %>%
  unnest(est)


#estimate bias
df_sims %>%
  select(model, dr:ipw) %>%
  pivot_longer(dr:ipw, names_to = 'estimator') %>%
  group_by( model, estimator) %>%
  summarise(bias = mean(value - true_patt))

#boxplots
df_sims %>%
  select(model, dr:ipw) %>%
  pivot_longer(dr:ipw, names_to='estimator') %>%
  mutate(error = value - true_patt) %>%
  ggplot(aes(y=error, x=model, color=estimator, fill=estimator)) +
  geom_boxplot(position = position_dodge(width = .5), alpha=.6) +
  geom_jitter(position = position_jitterdodge(jitter.width = .2,dodge.width = .5),
              alpha=0.1)+
  theme_minimal() +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2")

