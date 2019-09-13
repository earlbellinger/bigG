#### Show beta marginal posterior distribution 
#### Author: Earl Bellinger ( bellinger@phys.au.dk )
#### Stellar Astrophysics Centre Aarhus

source(file.path('..', 'scripts', 'utils.R'))
source(file.path('..', 'scripts', 'seismology.R'))

library(showtext)
font_add_google("Crimson Text", "grace")
#font_add_google("Source Sans Pro", "droid")
showtext_auto()

#attach(read.table('beta_samples.dat', col.names=c('beta')) / 10**-12)
#DF <- read.table('beta_samples.dat', header=1)
DF <- read.table('chain.dat', col.names=c(
    'n', 'tau', 'M', 'Y', 'Z', 'alpha', 'beta', 't_0', 'r1', 'r2', 'r3', 'r4', 'r5'))
#beta <- DF$beta / 10**-12

G_0 = 6.67408e-8
#t_0 = 13.799e9
G <- with(DF, G_0 * ( tau / t_0 )**beta)
deltaG <- (G - G_0) / G_0
print(mad(deltaG))

breaks <- 79
x <- with(DF, abs(beta)/(t_0*1e9)/1e-12)
hst <- hist(x, plot=F, breaks=breaks)
cdf <- with(hst, cumsum(counts) / max(cumsum(counts)))
upper.limit.bin <- which(cdf >= 0.95)[1]
upper.limit <- hst$breaks[upper.limit.bin]
print(upper.limit)

#breaks <- 10
x <- with(DF, abs(beta)/(t_0*1e9)/1e-12)
hst <- hist(x, plot=F, breaks=breaks)

plot_beta <- function(..., text.cex = 1, mgp = utils.mgp, mar = utils.mar,
        font = utils.font, tcl = utils.tcl) {
    par(mar = mar + c(0.3, 0, 1, 1), las = 1, cex.axis = text.cex,
        lwd=1.33,
        fig = c(0,1,0,1))
    
    #ymax <- signif(max(hst$counts) / sum(hst$counts) * 100, 2) / 100
    ymax <- 0.03
    #xlim <- c(0, 6)
    xlim <- c(0, 10)
    ylim <- c(0, sum(hst$counts) * ymax)
    
    plot(NA, axes = F, #log='x',
        xaxs = 'i',
        yaxs = 'i',
        xlim = xlim,
        ylim = ylim,
        xlab = "",
        ylab = ""
    )
    
    #upper.limit <- 
    cdf <- with(hst, cumsum(counts) / max(cumsum(counts)))
    upper.limit.bin <- which(cdf >= 0.95)[1]
    upper.limit <- hst$breaks[upper.limit.bin]
    
    col.pal <- adjustcolor(c(blue, orange), alpha.f=0.75)
    
    par(lwd=0.1)
    hist(x, add=T, border=F, breaks=breaks, 
        col=col.pal[ifelse(hst$breaks < upper.limit, 1, 2)])
    par(lwd=1.33)
    
    lines(c(upper.limit, upper.limit), c(0, hst$counts[upper.limit.bin-1]*0.97),
        lty=2, lwd=1.5)
    
    #lines(c(upper.limit, upper.limit), c(0, hst$counts[upper.limit.bin-1]*1.7),
    #    lty=2, lwd=1.5)
    
    #par(family='Helvetica LT Std Light')
    #text(upper.limit*1, hst$counts[upper.limit.bin-1]*1.55,#*0.92, 
    #     pos=3, labels="95% CI", cex=0.8*text.cex, col=1)
    #par(family=font)
    
    rect(xlim[2], ylim[1]-0.05, xlim[2]+0.05, ylim[2]+0.05,
        col='white', bg='white', border=NA)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1*text.cex,
        family=font, majorn=3, labels=T, lwd.ticks=par()$lwd, lwd=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.45, 0), cex.axis=1*text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=par()$lwd, 
        lwd=par()$lwd)
    
    yaxt <- paste0(pretty(c(0, ymax), n=4)*100, '%')
    par(mgp=mgp+c(0, 0.35, 0))
    axis(2, pretty(ylim, n=4), yaxt, tcl=0, cex.axis=1*text.cex)
    
    magaxis(3, tcl=tcl, mgp=mgp+c(0, 0.45, 0), cex.axis=1*text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=0, lwd=par()$lwd)
    magaxis(4, tcl=tcl, mgp=mgp+c(0, 0.45, 0), cex.axis=1*text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=0, lwd=par()$lwd)
    
    box(lwd=par()$lwd)
    
    par(xpd=NA)
    points(3.19, -365, pch=20, cex=0.1, col=1, lwd=0.2)
    par(xpd=F)
    
    mtext(expression("Marginalized posterior probability"), 
        2, 2.5, outer=F, las=0, cex=1.15*text.cex)
    mtext(expression(abs(G/G)~""~"["*10^-12~yr^-1*"]"), 
        1, 2.2, outer=F, cex=1.15*text.cex)
    #mtext(expression("Absolute rate of change of G"~"["*10^-12~yr^-1*"]"), 
    #    1, 2.2, outer=F, cex=1.15*text.cex)
    
}

make_plots(
    plot_beta,
    'beta_posterior',
    paper_pdf_width = 4.17309 * (1.2 * 0.95),
    paper_pdf_height = 4.17309 * 0.95,
    cex.paper = 1.12,
    use.cairo = T,
    thin=F, make_png=F, short=F, slides=F, 
    font = 'grace' 
)

if (F) {
breaks <- 20
x <- abs(deltaG)
hst <- hist(x, plot=F, breaks=breaks)

plot_deltaG <- function(..., text.cex = 1, mgp = utils.mgp, mar = utils.mar,
        font = utils.font, tcl = utils.tcl) {
    par(mar = mar + c(0.3, 0, 1, 1), las = 1, cex.axis = text.cex,
        lwd=1.33,
        fig = c(0,1,0,1))
    
    ymax <- round(max(hst$counts) / sum(hst$counts) * 100, -1) / 100
    xlim <- c(0, 3)
    ylim <- c(0, sum(hst$counts) * ymax)
    
    plot(NA, axes = F, 
        xaxs = 'i',
        yaxs = 'i',
        xlim = xlim,
        ylim = ylim,
        xlab = "",
        ylab = ""
    )
    
    cdf <- with(hst, cumsum(counts) / max(cumsum(counts)))
    upper.limit.bin <- which(cdf >= 0.95)[1]
    upper.limit <- hst$breaks[upper.limit.bin]
    
    col.pal <- adjustcolor(c(blue, orange), alpha.f=0.75)
    
    hist(x, add=T, border=F, breaks=breaks,
        col=col.pal[ifelse(hst$breaks < upper.limit, 1, 2)])
    
    lines(c(upper.limit, upper.limit), c(0, hst$counts[upper.limit.bin-1]*0.97),
        lty=2, lwd=1.5)
    
    par(family='Helvetica LT Std Light')
    text(upper.limit*1.04, hst$counts[upper.limit.bin-1]*0.92, 
         pos=3, labels="95% CI", cex=0.8*text.cex, col=1)
    par(family=font)
    
    rect(xlim[2], ylim[1]-0.05, xlim[2]+0.05, ylim[2]+0.05,
        col='white', bg='white', border=NA)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=0.8*text.cex,
        family=font, majorn=3, labels=T, lwd.ticks=par()$lwd, lwd=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=0.8*text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=par()$lwd, 
        lwd=par()$lwd)
    
    yaxt <- paste0(pretty(c(0, ymax), n=3)*100, '%')
    par(mgp=mgp+c(0, 0.35, 0))
    axis(2, pretty(ylim, n=3), yaxt, tcl=0, cex.axis=0.8*text.cex)
    
    magaxis(3, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=0.8*text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=0, lwd=par()$lwd)
    magaxis(4, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=0.8*text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=0, lwd=par()$lwd)
    
    box(lwd=par()$lwd)
    
    mtext(expression("Marginalized posterior probability"), 
        2, 2.5, outer=F, las=0, cex=1.15*text.cex)
    mtext(expression(abs(delta*G/G)), 
        1, 2.2, outer=F, cex=1.15*text.cex)
    
}

make_plots(
    plot_deltaG,
    'deltaG',
    paper_pdf_width = 4.17309 * (1.2 * 0.95),
    paper_pdf_height = 4.17309 * 0.95,
    cex.paper = 1.12,
    use.cairo = T,
    thin=F, make_png=F, short=F, slides=F, 
    font = 'grace' 
)
}

