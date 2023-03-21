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
gcomp = Ey1 - Ey0

#doubly robust estimator
g_as = glm(a ~ s*w, binomial, df)
g_s = glm(s ~ w, binomial, df)

g00 = (1-predict(g_as, newdata = df %>% mutate(s=0), type='response')) * (1-predict(g_s, newdata=df, type='response'))
g01 = (1-predict(g_as, newdata = df %>% mutate(s=1), type='response')) *    predict(g_s, newdata=df, type='response')
g10 =    predict(g_as, newdata = df %>% mutate(s=0), type='response')  * (1-predict(g_s, newdata=df, type='response'))
g11 =    predict(g_as, newdata = df %>% mutate(s=1), type='response')  *    predict(g_s, newdata=df, type='response')

p00 = mean( (1-df$a)*(1-df$s)  )
p10 = mean(     df$a*(1-df$s)  )


I11 = 1*( df$a==1 & df$s==1  )
I01 = 1*( df$a==0 & df$s==1  )
I10 = 1*( df$a==1 & df$s==0  )

mm1 = predict(m1, newdata=df)
mm0 = predict(m0, newdata=df)

eif = with(df,  I11*g10*(y1 - y0 - mm1)/(p10*g11) - I01*g00*(y1 - y0 - mm0)/(p00*g01) + I10*(mm1 - mm0)/p10 )
dr = mean(eif)


