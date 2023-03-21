if(!require(didsim)) {
  remotes::install_github('audreyrenson/didsim')
  library(didsim)
}
library(tidyverse)
library(simstudy)


set.seed(5)

#generate data
def = defData(varname = 's', formula = '0.1', dist = 'binary', variance = .1) %>%
  defData(varname = 'u', formula = 0, variance=1) %>%
  defData(varname = 'w', formula = add_coefs_to_formula('s + u'), variance=.1) %>%
  defData(varname = 'a', formula = add_coefs_to_formula('s + w + u'), dist='binary', link='logit') %>%
  defData(varname = 'y0', formula = add_coefs_to_formula('w + u', u_coef = 1), variance=.1) %>%
  defData(varname = 'y1', formula = add_coefs_to_formula('w + a + w:a + u', u_coef=1), variance=.1)


#estimate true PATT E[Y1(1)-Y1(0)|A=1, S=0]
def0 = def %>% mutate(formula = str_replace_all(formula, 'a', '0'))
def1 = def %>% mutate(formula = str_replace_all(formula, 'a', '1'))

df0 = genData(n=1e5, def0)
df1 = genData(n=1e5, def1)

true_patt = mean(df1$y1[df1$a==1 & df1$s==0]) - mean(df0$y1[df0$a==1 & df0$s==0])
true_patt

#g-computation
df = genData(n=1e5, def)


m1 = lm(y1-y0 ~ w, data=df %>% filter(a==1, s==1))
m0 = lm(y1-y0 ~ w, data=df %>% filter(a==0, s==1))

Ey1 = mean(predict(m1, newdata=df %>% filter(s==0, a==1)))
Ey0 = mean(predict(m0, newdata=df %>% filter(s==0, a==1)))
Ey1 - Ey0
