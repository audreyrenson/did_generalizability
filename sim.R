if(!require(didsim)) {
  remotes::install_github('audreyrenson/didsim')
  library(didsim)
}
library(tidyverse)
library(simstudy)

def = defData(varname = 's', formula = '0.1', dist = 'binary', variance = 1) %>%
  defData(varname = 'u', formula = 0, variance=1) %>%
  defData(varname = 'w', formula = add_coefs_to_formula('s + u'), variance=1) %>%
  defData(varname = 'a', formula = add_coefs_to_formula('s + w + u'), dist='binary', link='logit') %>%
  defData(varname = 'y0', formula = add_coefs_to_formula('w + u', u_coef = 1), variance=1) %>%
  defData(varname = 'y1', formula = add_coefs_to_formula('w + a + w:a + u', u_coef=1), variance=1)
