library(cowplot)
library(dplyr)

solarized <- c(base03 = '#002b36',
               base02 = '#073642',
               base01 = '#486e75',
               base00 = '#657b83',
               base0  = '#839496',
               base1  = '#93a1a1',
               base2  = '#eee8d5',
               base3  = '#fdf6e3',
               yellow = '#b58900',
               orange = '#cb4b16',
               red    = '#dc322f',
               magenta= '#d33682',
               violet = '#6c71c4',
               blue   = '#268bd2',
               cyan   = '#2aa198',
               green  = '#859900')

solpal5 <- solarized[c('red','blue','yellow','base02','base1')]
solpal5 <- unname(solpal5)
scale_colour_viridis = scale_colour_viridis_d(option="C", begin=.1, end=.9, direction= -1)

munsell <- function(x, ref) ifelse( (x/ref) <= (6. / 29) ^ 3,
                                    11.6 * 841. / 108 * x/ref + 4. / 29 - 1.6,
                                    11.6* (x/ref)^(1./3)                - 1.6)

function.zoo <- list('id'        = function(x) x,
                     'square'    = function(x) x^2,
                     'sqrt'      = function(x) x^0.5,
                     'saddle'    = function(x) 4*(x-0.5)^3 + 0.5,
                     'invsaddle' = function(x) sign(x/4 - 0.125) * (abs(x/4 - 0.125))^(1/3) + 0.5,
                     'munsell'   = function(x) (munsell(x, 1.0)-munsell(0.0, 1.0))/(munsell(1.0, 1.0)-munsell(0.0, 1.0)),
                     'log'       = function(x) (log(x+0.05) - log(0.05))/(log(1.05) - log(0.05)),
                     'exp'       = function(x) 0.05 * ( (1.05/0.05)^x-1),
                     'invlogit'  = function(x) atanh((tanh(5*0.5) - tanh(-0.5*5))*x + tanh(-0.5*5))/5 + 0.5,
                     'logit'     = function(x) (tanh(5*(x-0.5)) - tanh(-0.5*5))/(tanh(5*0.5) - tanh(-0.5*5))

)

lps.files <- list.files(path = "./data", pattern = 'lps[[:print:]]+.csv')

lps.priors.df <- data.frame()

for (f in lps.files) {
  
  dat <- read.csv(paste0('data/', f), sep = '\t')
  
  lps.priors.df <- rbind(lps.priors.df, dat)
}

sensitivity.files <- list.files(path = "./data", pattern = 'prec[[:print:]]+.csv')

sensitivity.priors.df <- data.frame()

for (f in sensitivity.files) {
  
  dat <- read.csv(paste0('data/', f), sep = '\t')
  
  sensitivity.priors.df <- rbind(sensitivity.priors.df, dat)
}