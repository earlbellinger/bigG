#### CD diagram for tracks of different values of \beta 
#### Author: Earl Patrick Bellinger ( bellinger@phys.au.dk )
#### Stellar Astrophysics Centre Aarhus

source(file.path('..', 'scripts', 'utils.R'))
source(file.path('..', 'scripts', 'seismology.R'))

library(showtext)
font_add_google("Crimson Text", "grace")
font_add_google("Source Sans Pro", "droid")
showtext_auto()

DFs <- Map(function(filename) 
        read.table(file.path('small_grid', filename), header=1),
    filename=list.files('small_grid'))

plot_CD <- function(..., text.cex = 1, mgp = utils.mgp, mar = utils.mar,
        font = utils.font, tcl = utils.tcl) {
    par(mar = mar + c(0.3, -0.5, 2, 1), las = 1, cex.axis = text.cex,
        lwd=1.33,
        fig = c(0,1,0,1))
    
    xlim <- c(100, 260)
    ylim <- c(5, 22)
    
    plot(NA, axes = F,
        xaxs = 'i',
        yaxs = 'i',
        xlim = xlim,
        ylim = ylim,
        xlab = "",
        ylab = ""
    )
    
    for (ii in 1:length(DFs)) {
        DF <- DFs[[ii]]
        
        ## recompute the dnu02s 
        ell0 <- DF[,grepl('nu_0', names(DF))]
        #print(ell0)
        ell2 <- DF[,grepl('nu_2', names(DF))]
        ell0n <- sapply(1:ncol(ell0), function(jj) strsplit(names(ell0)[jj], 'nu_0_')[[1]][2])
        ell2n <- sapply(1:ncol(ell2), function(jj) strsplit(names(ell2)[jj], 'nu_2_')[[1]][2])
        #print(ell0n)
        dnus <- c()
        for (jj in 1:nrow(DF)) {
            freqs <- rbind(data.frame(l=rep(0, ncol(ell0)), n=ell0n, nu=as.numeric(ell0[jj,])),
                           data.frame(l=rep(2, ncol(ell2)), n=ell2n, nu=as.numeric(ell2[jj,])))
            freqs <- freqs[as.numeric(freqs$nu) > DF$nu_max[jj] - 8.5*DF$Dnu0[jj] &
                           as.numeric(freqs$nu) < DF$nu_max[jj] + 8.5*DF$Dnu0[jj],]
            seis.DF <- seismology(freqs, nu_max=DF$nu_max[jj])
            dnus <- c(dnus, seis.DF$dnu02)
        }
        
        DFs[[ii]]$dnu02 <- dnus
    }
    
    
    par(family='Helvetica LT Std Light')#'Droid')#
    for (ii in 1:length(DFs)) {
        DF <- DFs[[ii]]
        name <- as.numeric(strsplit(names(DFs)[ii], '.dat')[[1]][1])
        with(DF, lines(Dnu0, dnu02, lwd=1.33, 
            col=ifelse(name>0, orange, ifelse(name<0, blue, 1)),
            lty=1))
    }
    for (ii in 1:length(DFs)) {
        DF <- DFs[[ii]]
        name <- as.numeric(strsplit(names(DFs)[ii], '.dat')[[1]][1])
        with(DF, points(Dnu0[1], dnu02[1], pch=20, cex=0.1)) 
        #with(DF, points(Teff[length(DF)], L[length(DF)], pch=1, cex=0.5)) 
        if (name == 0.05)
            with(DF, text(Dnu0[1]-0, dnu02[1], 
                pos=2, labels=bquote(.(name)), cex=0.8*text.cex, col=orange))
        if (name == 0)
            with(DF, text(Dnu0[1]+10, dnu02[1]-0.5, 
                pos=3, labels=expression(beta==0), cex=0.8*text.cex))
        if (name == -0.05)# || name == 0 || name == 0.05) {
            with(DF, text(Dnu0[1] - 2, dnu02[1] - 0.1, 
                pos=4, labels=bquote(.(name)), cex=0.8*text.cex, col=blue))
    }
    par(family=font)
    
    with(DFs[['0.00.dat']], lines(Dnu0, dnu02, lwd=1.66))
    
    xlimi <- c(173.541-0.064, 173.541+0.064)
    ylimi <- c(7.901 - 0.167, 7.901 + 0.167)
    segments(xlimi[1], ylimi[1], 121.1, 16.2, lty=2, col='gray', lwd=par()$lwd)
    segments(xlimi[2], ylimi[2], 156,   20.7, lty=2, col='gray', lwd=par()$lwd)
    rect(xlimi[1], ylimi[1], xlimi[2], ylimi[2])
    
    magaxis(1:2, tcl=tcl, labels=F, lwd.ticks=par()$lwd)
    #magaxis(3:4, tcl=tcl, labels=F, lwd.ticks=0)
    
    axis(1, tick=F, tcl=0, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.35, 0), lwd=par()$lwd)
    
    axis(2, tick=F, tcl=0, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.45, 0), lwd=par()$lwd)
    
    box(lwd=par()$lwd)
    
    mtext(expression("Small separation"~delta*nu[0][2]/mu*Hz), 
        2, 2, outer=F, las=0, cex=1.15*text.cex)
    mtext(expression("Large frequency separation"~Delta*nu/mu*Hz),
        1, 2.2, outer=F,        cex=1.15*text.cex)
    
    
    par(fig = c(0.1, 0.5, 0.45, 0.95), new = T)
    
    plot(NA, axes = F,
        xaxs = 'i', yaxs = 'i',
        xlim = xlimi, ylim = ylimi,
        xlab = "", ylab = "")
    for (ii in 1:length(DFs)) {
        DF <- DFs[[ii]]
        name <- as.numeric(strsplit(names(DFs)[ii], '.dat')[[1]][1])
        with(DF, lines(Dnu0, dnu02, lwd=1.33, 
            col=ifelse(name>0, orange, ifelse(name<0, blue, 1)),
            lty=1))
    }
    for (ii in 1:length(DFs)) {
        DF <- DFs[[ii]]
        name <- as.numeric(strsplit(names(DFs)[ii], '.dat')[[1]][1])
        with(DF, points(Dnu0[1], dnu02[1], pch=20, cex=0.1)) 
    }
    with(DFs[['0.00.dat']], lines(Dnu0, dnu02, lwd=1.66))
    magaxis(1, tcl=tcl, labels=F, lwd.ticks=par()$lwd, majorn=2, 
        mgp=mgp+c(0, 0.3, 0), family=font, cex.axis=0.8*text.cex)
    magaxis(2, tcl=tcl, labels=T, lwd.ticks=par()$lwd, majorn=2, 
        mgp=mgp+c(0, 0.4, 0), family=font, cex.axis=0.8*text.cex)
    par(mgp=mgp+c(0, 0.3, 0))
    axis(1, pretty(xlimi, n=2), tick=F, cex.axis=0.8*text.cex)
    #axis(1, tick=F, tcl=0, 
    #    cex.axis=text.cex, tcl=0, las=1, 
    #    mgp=mgp+c(0, 0.35, 0), lwd=par()$lwd)
    #axis(2, tick=F, tcl=0, 
    #    cex.axis=text.cex, tcl=0, las=1, 
    #    mgp=mgp+c(0, 0.45, 0), lwd=par()$lwd)
    
    box(lwd=par()$lwd)
    
}

make_plots(
    plot_CD,
    'CD',
    #paper_pdf_width = 4.17309 * (1.2 * 0.8),
    #paper_pdf_height = 4.17309 * 0.8,
    paper_pdf_width = 4.17309 * (1.2 * 0.95),
    paper_pdf_height = 4.17309 * 0.95,
    cex.paper = 1.12,
    use.cairo = T,
    thin=F, make_png=F, short=F, slides=F, 
    font = 'grace' #'Palatino Linotype'
)

