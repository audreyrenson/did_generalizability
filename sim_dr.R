if(!require(didsim)) {
  remotes::install_github('audreyrenson/didsim')
  library(didsim)
}
library(tidyverse)
library(simstudy)
library(furrr)


set.seed(seed <- 1451456)

plan(multisession)


#data generating model
def = defData(varname = 'u', formula = 0, variance=1) %>%
  defData(varname = 'w', formula = add_coefs_to_formula('u'), variance=.1) %>%
  defData(varname = 's', formula = '0.4 + 0.2*pnorm(w)', dist='binary',link='identity') %>%
  defData(varname = 'a', formula = '0.3 + 0.1*(1-s) - 0.1*pnorm(w) + 0.1*pnorm(w)*s - 0.1*pnorm(u)', dist='binary', link='identity') %>%
  defData(varname = 'y0', formula = add_coefs_to_formula('sin(w) + u', u_coef = .2), variance=.01) %>%
  defData(varname = 'y1', formula = add_coefs_to_formula('sin(w) + a + sin(w):a + u', u_coef=.2), variance=.01)

#estimate true PATT E[Y1(1)-Y1(0)|A=1, S=0]
def0 = def %>% mutate(formula = str_replace_all(formula, 'a', '0'))
def1 = def %>% mutate(formula = str_replace_all(formula, 'a', '1'))

df0 = genData(n=1e6, def0)
df1 = addColumns(def1 %>% filter(varname=='y1'), df0 %>% select(-y1))

true_patt = mean(  (df1$y1 - df0$y1)[df1$a==1& df1$s==0] )

var(df1$y1)
table(df1$a, df1$s) %>% prop.table()


# doubly robust estimator -------------------------------------------------


dr = function(data, model) {

  stopifnot(model %in% c('true','gfal','qfal','bfal'))

  formula_m = if(model %in% c('true','gfal')) y1-y0~sin(w)*a else y1-y0 ~a+w
  formula_g_as = if(model %in% c('true','qfal')) a~pnorm(w)*s else a~w+s
  formula_g_s = if(model %in% c('true','qfal')) s~pnorm(w) else s~w

  m_a = lm(formula_m, data=data)
  m1 = predict(m_a, newdata=data %>% mutate(a=1))
  m0 = predict(m_a, newdata=data %>% mutate(a=0))

  g_as = lm(formula_g_as, data)
  g_s = lm(formula_g_s, data)

  g01 = (1-predict(g_as, newdata = data %>% mutate(s=1), type='response')) *    predict(g_s, newdata=data, type='response')
  g10 =    predict(g_as, newdata = data %>% mutate(s=0), type='response')  * (1-predict(g_s, newdata=data, type='response'))
  g11 =    predict(g_as, newdata = data %>% mutate(s=1), type='response')  *    predict(g_s, newdata=data, type='response')

  p10 = mean(     data$a*(1-data$s)  )


  I11 = 1*( data$a==1 & data$s==1  )
  I01 = 1*( data$a==0 & data$s==1  )
  I10 = 1*( data$a==1 & data$s==0  )

  y1 = data$y1
  y0 = data$y0

  eif =  I11*g10*(y1 - y0 - m1)/(p10*g11) - I01*g10*(y1 - y0 - m0)/(p10*g01) + I10*(m1 - m0)/p10

  tibble(dr_est=mean(eif),
         dr_se = sqrt(var(eif)/nrow(data)),
         gcomp = mean( I10*(m1 - m0)/p10 ),
         ipw = mean(   ( I11*g10 / (p10*g11)  - I01*g10/(p10*g01) )*(y1-y0)   ))
}


# sim ---------------------------------------------------------------------

nsims=200

df_sims = tibble(sim = 1:nsims,
                 data = future_map(sim, ~genData(1e4, def), .options = furrr_options(seed=seed))) %>%
  expand_grid(model = c('true','gfal','qfal','bfal'))

df_sims_dr = df_sims %>%
  mutate(est = future_map2(data, model, dr)) %>%
  unnest(est) %>%
  mutate(error = dr_est - true_patt)





pdf(file='boxplots.pdf', width = 5, height=2.3)
df_sims_dr %>%
  select(model, double_robust=dr_est, gcomp, iow=ipw) %>%
  pivot_longer(-model, names_to='estimator') %>%
  mutate(error = value - true_patt) %>%
  ggplot(aes(y=error, x=model, color=estimator, fill=estimator)) +
  geom_boxplot(position = position_dodge(width = .5), alpha=.5) +
  theme_minimal()
dev.off()

df_sims_dr %>%
  group_by(model) %>%
  summarise(bias = mean(dr_est) - true_patt,
            median_bias = median(dr_est) - true_patt)

#positivity issues?
data = df_sims_dr %>%
  group_by(model) %>%
  filter(dr_est == min(dr_est)) %>%
  {.$data[[1]]}


df_sims_dr %>%
  mutate(error = gcomp - true_patt) %>%
  ggplot(aes(x=error)) +
  facet_wrap(~model) +
  geom_histogram() +
  ggtitle('gcomp')

df_sims_dr %>%
  filter(model == 'true') %>%
  mutate(error = ipw - true_patt) %>%
  ggplot(aes(x=error)) +
  geom_histogram() +
  ggtitle('ipw')
