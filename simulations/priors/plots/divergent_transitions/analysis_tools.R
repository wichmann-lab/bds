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

grid.eval <- function(name, diff_scale) {
  log_posterior <- function(x) {
    pars <- list(lapses=x[1], sensitivity=x[2], psi_diff=x[3:length(x)])
    lp <- tryCatch(-log_prob(diff_scale$stanfit, unconstrain_pars(diff_scale$stanfit, pars)),error=function(cond) return(Inf))
    
    return(lp)
  }

  post.draws <- extract(diff_scale$stanfit,c("lapses", "sensitivity", "psi", "psi_diff"))
  sens <- seq(2.5, max(post.draws$sensitivity, 30), len=100)
  lps <- seq(0.001, max(post.draws$lapses, .1), len=100)
  lps.step <- lps[2] - lps[1]
  sens.step <- sens[2] - sens[1]

  v <- expand.grid(lps, sens)
  
  sc.lp <- data.frame()
  
  lp <- function(sl, psi_diff=pmean[3:length(pmean)]) {
    x <- c(sl[[1]], sl[[2]], psi_diff)
    -log_posterior(x)
  }
  
  grid.lps_sens <- function (id) {
    if (!file.exists(paste0('eval/', name, '_', id, '.tsv'))) {
    sc <- post.draws$psi_diff[id,]
    v$lp <- apply(v, 1, lp, psi_diff = sc)
    
    write.table(v, file=paste0('eval/', name, '_', id, '.tsv'), sep='\t')
    } else {
        v <- read.table(paste0('eval/', name, '_', id, '.tsv'), sep="\t")
    }
    
    matrixStats::logSumExp(v$lp+log(lps.step)+log(sens.step))
  }
  
  sc.lp <- mcmapply(grid.lps_sens, 1:length(post.draws$psi_diff[,1]), mc.cores=detectCores())
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
    geom_raster(data=v, aes(x=Var1, y=Var2, fill=-lp)) +
    geom_contour(data=v, aes(x=Var1, y=Var2, z=-lp), binwidth=5, colour=solpal5[5]) +
    geom_point(data=data.frame(Var1=post.draws$lapses, Var2=post.draws$sensitivity), aes(x=Var1, y=Var2), alpha=0.1) +
    geom_point(data=data.frame(Var1=map[1], Var2=map[2]), aes(x=Var1, y=Var2), color=solpal5[1]) +
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
      geom_raster(data=v, aes(x=Var1, y=Var2, fill=lp)) +
      geom_contour(data=v, aes(x=Var1, y=Var2, z=lp), binwidth=5, colour=solpal5[5]) +
      geom_point(data=data.frame(Var1=c(map[1], pmean[1], map2[1], smp[1]), Var2=c(map[2], pmean[2], map2[2], smp[2]), method=c('MAP', 'mean', 'local max', 'HMC sample')), aes(x=Var1, y=Var2, colour=method, shape=method), size=3) +
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

  post.lp <- data.frame()
  for (id in order(sc.lp)) {
    v <- read.table(paste0('eval/', name, '_', id, '.tsv'), sep="\t")
    
    if (nrow(post.lp) == 0) {
      post.lp <- v
    } else {
      post.lp$lp <- log(exp(post.lp$lp) + exp(v$lp))
    }
  }
  
  means <- c(lps.mean=sum(exp(post.lp$lp)/sum(exp(post.lp$lp))*post.lp$Var1), sens.mean=sum(exp(post.lp$lp)/sum(exp(post.lp$lp))*post.lp$Var2))
  
  p <- ggplot() + 
    geom_raster(data=post.lp, aes(x=Var1, y=Var2, fill=lp)) +
    geom_point(data=data.frame(Var1=post.draws$lapses, Var2=post.draws$sensitivity), aes(x=Var1, y=Var2), alpha=0.1) +
    geom_contour(data=post.lp, aes(x=Var1, y=Var2, z=lp), binwidth=5, colour=solpal5[5]) +
    geom_point(data=data.frame(Var1=c(map[1], pmean[1], means[1]), Var2=c(map[2], pmean[2], means[2]), method=c('MAP', 'sample mean', 'posterior mean')), aes(x=Var1, y=Var2, colour=method, shape=method), size=3) +
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