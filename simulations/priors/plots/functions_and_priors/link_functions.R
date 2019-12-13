library(cowplot)
library(dplyr)
library(reshape2)

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

a <- 0
b <- c(2.5, 5.0, 7.5, 10.0)
c <- seq(0, 20, len=101)


# a=0, b=a+z, c=a+x+z  =>  c - 2b = x - z
henaff <- function(x, z) pnorm((x-z)/sqrt(6))*pnorm(x+z/sqrt(2)) + pnorm(-(x-z)/sqrt(6))*pnorm(-x-z/sqrt(2))

gr <- expand.grid(b,c)
d <- data.frame(x = gr$Var2, g=gr$Var1, decision_noise=pnorm((gr$Var2 - gr$Var1)/2), sensory_noise=henaff(gr$Var2, gr$Var1), logit=1/(1+exp(-gr$Var2+gr$Var1)))

d <- melt(d, id.vars=c('x', 'g'))
ggplot(d %>% filter(g==10.0), aes(x=x, y=value, colour=variable)) +
  scale_colour_manual(values= solpal5) +
  geom_line(size=1)
