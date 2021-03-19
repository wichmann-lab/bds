library(bds)
library(reshape2)
library(dplyr)
library(parallel)
library(cowplot)

options(mc.cores=parallel::detectCores())

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

# pmean <- c(lapses=diff_scale$lapserate,sensitivity=diff_scale$sensitivity,diff(diff_scale$scale))

#log_posterior <- function(x) {
#  pars <- list(lapses=x[1], sensitivity=x[2], psi_diff=x[3:length(x)])
#  lp <- tryCatch(-log_prob(diff_scale$stanfit, unconstrain_pars(diff_scale$stanfit, pars)),error=function(cond) return(Inf))
  
#  return(lp)
#}

load.grid <- function(name, diff_scale) {
  log_posterior <- function(x) {
    pars <- list(lapses=x[1], sensitivity=x[2], psi_diff=x[3:length(x)])
    lp <- tryCatch(-log_prob(diff_scale$stanfit, unconstrain_pars(diff_scale$stanfit, pars)),error=function(cond) return(Inf))
    
    return(lp)
  }

  post.draws <- extract(diff_scale$stanfit,c("lapses", "sensitivity", "psi", "psi_diff"))
  sens <- seq(5, max(post.draws$sensitivity, 30), len=30)
  lps <- seq(0.001, max(post.draws$lapses, .1), len=30)
  lps.step <- lps[2] - lps[1]
  sens.step <- sens[2] - sens[1]

  if (!file.exists(paste0('eval/', name, '.tsv'))) {
    lp <- grid.eval(diff_scale, sensitivities = sens, lapses=lps)
    
    write.table(lp, file=paste0('eval/', name, '.tsv'), sep='\t')
  } else {
    lp <- read.table(paste0('eval/', name, '.tsv'), sep="\t")
  }
  
  marginal.lp <- lp %>% group_by(lapses, sensitivity) %>% summarise(density=matrixStats::logSumExp(density)+log(lps.step)+log(sens.step))

  marginal.lp
}

plot.diagnostics <- function(name, diff_scale) {
  log_posterior <- function(x) {
    pars <- list(lapses=x[1], sensitivity=x[2], psi_diff=x[3:length(x)])
    lp <- tryCatch(-log_prob(diff_scale$stanfit, unconstrain_pars(diff_scale$stanfit, pars)),error=function(cond) return(Inf))
    
    return(lp)
  }
  
  lp <- function(sl, psi_diff=pmean[3:length(pmean)]) {
    x <- c(sl[[1]], sl[[2]], psi_diff)
    log_posterior(x)
  }
  
  pmean <- c(lapses=diff_scale$lapserate,sensitivity=diff_scale$sensitivity,diff(diff_scale$scale))
  post.draws <- extract(diff_scale$stanfit,c("lapses", "sensitivity", "psi", "psi_diff"))
  
  map <- optim(pmean, log_posterior)$par
  map.local <- optim(map[1:2], lp)$par
  map2.fit <- optim(c(map.local, map[3:length(map)]), log_posterior)
  map <- map2.fit$par
  
  p.sens <- qplot(post.draws$sensitivity, binwidth=0.025) +
    geom_vline(data=data.frame(method=c("mean", "MAP"), val=c(pmean[2], map[2])), aes(xintercept=val, colour=method), size=1) +
    xlab("sensitivity") +
    scale_colour_manual(values = solpal5) +
    theme_classic()
  
  p.lps <- qplot(post.draws$lapses, binwidth=0.0025) +
    geom_vline(data=data.frame(method=c("mean", "MAP"), val=c(pmean[1], map[1])), aes(xintercept=val, colour=method), size=1) +
    xlab("lapserate") +
    xlim(0, NA) +
    scale_colour_manual(values = solpal5) +
    theme_classic()
  
  lp.sens <- function(s) {
    x <- map
    x[2] <- s
    -log_posterior(x)
  }
  
  sens <- seq(2.5, max(post.draws$sensitivity, 30), len=100)
  
  y <- exp(sapply(sens, lp.sens))
  
  p.sens.dens <- ggplot() +
    geom_line(data=data.frame(x=sens, y=y/sum(y*(sens[2]-sens[1]))), aes(x=x,y=y), color=solpal5[1]) +
    geom_density(data=data.frame(x=post.draws$sensitivity),aes(x=x), color=solpal5[2]) +
    geom_vline(data=data.frame(method=c("mean", "MAP"), val=c(pmean[2], map[2])), aes(xintercept=val, colour=method), size=1) +
    xlab('sensitivity') +
    ylab('posterior density') +
    scale_colour_manual(values = solpal5) +
    theme_classic()
  
  lp.lps <- function(l) {
    x <- map
    x[1] <- l
    -log_posterior(x)
  }
  
  lps <- seq(0.001, max(post.draws$lapses, .1), len=100)
  
  y <- exp(sapply(lps, lp.lps))
  
  p.lps.dens <- ggplot() +
    geom_line(data=data.frame(x=lps, y=y/sum(y*(lps[2]-lps[1]))), aes(x=x,y=y), color=solpal5[1]) +
    geom_density(data=data.frame(x=post.draws$lapses),aes(x=x), color=solpal5[2]) +
    geom_vline(data=data.frame(method=c("mean", "MAP"), val=c(pmean[1], map[1])), aes(xintercept=val, colour=method), size=1) +
    xlab('lapserate') +
    ylab('posterior density') +
    scale_colour_manual(values = solpal5) +
    theme_classic()
  
  v <- expand.grid(lps, sens)
  
  v$lp <- apply(v, 1, lp)
  
  fill_scale <- scale_fill_viridis_c(name="log likelihood", limits=c(-map2.fit$value - 20, -map2.fit$value), option="inferno", oob=scales::squish)
  
  p.density <- ggplot() + 
    geom_raster(data=v, aes(x=lapses, y=sensitivity, fill=-lp)) +
    geom_contour(data=v, aes(x=lapses, y=sensitivity, z=-lp), binwidth=5, colour=solpal5[5]) +
    geom_point(data=data.frame(lapses=post.draws$lapses, sensitivity=post.draws$sensitivity), aes(x=lapses, y=sensitivity), alpha=0.1) +
    geom_point(data=data.frame(lapses=map[1], sensitivity=map[2]), aes(x=lapses, y=sensitivity), color=solpal5[1]) +
    xlab('lapserate') +
    ylab('sensitivity') +
    theme_classic() +
    fill_scale
  
  ggsave(paste0(name, "_sensitivity.pdf"), p.sens)
  ggsave(paste0(name, "_lapses.pdf"), p.lps)
  ggsave(paste0(name, "_sens_dens.pdf"), p.sens.dens)
  ggsave(paste0(name, "_lps_dens.pdf"), p.lps.dens)
  ggsave(paste0(name, "_density.pdf"), p.density)
}

animate.posterior_landscape <- function(name, diff_scale, sc.lp) {
  log_posterior <- function(x) {
    pars <- list(lapses=x[1], sensitivity=x[2], psi_diff=x[3:length(x)])
    lp <- tryCatch(-log_prob(diff_scale$stanfit, unconstrain_pars(diff_scale$stanfit, pars)),error=function(cond) return(Inf))
  
    return(lp)
  }
  
  lp <- function(sl, psi_diff=pmean[3:length(pmean)]) {
    x <- c(sl[[1]], sl[[2]], psi_diff)
    log_posterior(x)
  }
  
  pmean <- c(lapses=diff_scale$lapserate,sensitivity=diff_scale$sensitivity,diff(diff_scale$scale))
  post.draws <- extract(diff_scale$stanfit,c("lapses", "sensitivity", "psi", "psi_diff"))

  map <- optim(pmean, log_posterior)$par
  map.local <- optim(map[1:2], lp)$par
  map2.fit <- optim(c(map.local, map[3:length(map)]), log_posterior)
  map <- map2.fit$par

  fill_scale <- scale_fill_viridis_c(name="log likelihood", limits=c(-map2.fit$value - 20, -map2.fit$value), option="inferno", oob=scales::squish)
  
  ct <- 0
  for (id in order(sc.lp)) {
    v <- read.table(paste0('eval/', name, '_', id, '.tsv'), sep="\t")
    
    map2 <- optim(map[1:2], function (x) lp(x, psi_diff=post.draws$psi_diff[id,]))$par
    
    smp <- c(post.draws$lapses[id], post.draws$sensitivity[id])
    p <- ggplot() + 
      geom_raster(data=v, aes(x=lapses, y=sensitivity, fill=lp)) +
      geom_contour(data=v, aes(x=lapses, y=sensitivity, z=lp), binwidth=5, colour=solpal5[5]) +
      geom_point(data=data.frame(lapses=c(map[1], pmean[1], map2[1], smp[1]), sensitivity=c(map[2], pmean[2], map2[2], smp[2]), method=c('MAP', 'mean', 'local max', 'HMC sample')), aes(x=lapses, y=sensitivity, colour=method, shape=method), size=3) +
      xlab('lapserate') +
      ylab('sensitivity') +
      ggtitle(sc.lp[id]) +
      fill_scale +
      scale_colour_manual(name="point estimates", values=solpal5[c(4,5,1,2)]) +
      scale_shape_manual(name="point estimates", values=c(4, 1, 17, 15)) +
      theme_classic()
  
    ggsave(paste0('anim/',name, "_", ct, '.jpg'), p)
    ct <- ct + 1
  }
  
#  system(paste0("ffmpeg -framerate 25 -i anim/", name,"%004d.jpg-c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p ", name, ".mp4 "))
}

summarize.posterior <- function(name, diff_scale, sc.lp, postfix='_full') {
  log_posterior <- function(x) {
    pars <- list(lapses=x[1], sensitivity=x[2], psi_diff=x[3:length(x)])
    lp <- tryCatch(-log_prob(diff_scale$stanfit, unconstrain_pars(diff_scale$stanfit, pars)),error=function(cond) return(Inf))
  
    return(lp)
  }
  
  lp <- function(sl, psi_diff=pmean[3:length(pmean)]) {
    x <- c(sl[[1]], sl[[2]], psi_diff)
    log_posterior(x)
  }
  
  pmean <- c(lapses=diff_scale$lapserate,sensitivity=diff_scale$sensitivity,diff(diff_scale$scale))
  post.draws <- extract(diff_scale$stanfit,c("lapses", "sensitivity", "psi", "psi_diff"))

  map <- optim(pmean, log_posterior)$par
  map.local <- optim(map[1:2], lp)$par
  map2.fit <- optim(c(map.local, map[3:length(map)]), log_posterior)
  map <- map2.fit$par

  fill_scale <- scale_fill_viridis_c(name="log likelihood", option="inferno")

  post.lp <- load.grid(name, diff_scale)
  
  full_dens <- sum(exp(post.lp$density))
  means <- c(lps.mean=sum(exp(post.lp$density)/full_dens*post.lp$lapses), sens.mean=sum(exp(post.lp$density)/full_dens*post.lp$sensitivity))
  
  p <- ggplot() + 
    geom_raster(data=post.lp, aes(x=lapses, y=sensitivity, fill=density)) +
    geom_point(data=data.frame(lapses=post.draws$lapses, sensitivity=post.draws$sensitivity), aes(x=lapses, y=sensitivity), alpha=0.1) +
    geom_contour(data=post.lp, aes(x=lapses, y=sensitivity, z=density), binwidth=5, colour=solpal5[5]) +
    geom_point(data=data.frame(lapses=c(map[1], pmean[1], means[1]), sensitivity=c(map[2], pmean[2], means[2]), method=c('MAP', 'sample mean', 'posterior mean')), aes(x=lapses, y=sensitivity, colour=method, shape=method), size=3) +
    xlab('lapserate') +
    ylab('sensitivity') +
    ggtitle(name) +
    fill_scale +
    scale_colour_manual(name="point estimates", values=solpal5[c(5,2,1)]) +
    scale_shape_manual(name="point estimates", values=c(1, 4, 17)) +
    theme_classic()
  
  ggsave(paste0('marginal/', name, postfix, ".pdf"), p)

  means
}