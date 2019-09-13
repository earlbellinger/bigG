#### HR diagram for tracks of different values of \beta 
#### Author: Earl Patrick Bellinger ( bellinger@phys.au.dk )
#### Stellar Astrophysics Centre Aarhus

source(file.path('..', 'scripts', 'utils.R'))

library(showtext)
font_add_google("Crimson Text", "grace")
font_add_google("Source Sans Pro", "droid")
showtext_auto()

DFs <- Map(function(filename) 
        read.table(file.path('small_grid', filename), header=1),
    filename=list.files('small_grid'))

plot_HR <- function(...,
        text.cex = 1, mgp = utils.mgp, mar = utils.mar, font = utils.font,
        tcl = utils.tcl) {
    par(mar = mar + c(0.3, -0.5, 2, 1),
        las = 1,
        cex.axis = text.cex,
        lwd=1.33,
        fig = c(0,1,0,1)
    )
    
    xlim <- c(6500, 5000)
    ylim <- c(0.2, 1.6)
    
    plot(NA, axes = F,
        #log = 'xy',
        xaxs = 'i',
        yaxs = 'i',
        xlim = xlim,
        ylim = ylim,
        xlab = "",
        ylab = ""
    )
    
    spectral.divs <- c(30000, 10000, 7500, 6000, 5200, 3700, 2400)
    rs <- c(175 / 255, 199 / 255, 1, 1, 1, 1, 1, 1)
    gs <- c(201, 216, 244, 229, 217, 199, 166) / 255
    bs <-
        c(1, 1, 243 / 255, 207 / 255, 178 / 255, 142 / 255, 81 / 255)
    cols <- c(
        rgb(175 / 255, 201 / 255, 1),
        # O
        rgb(199 / 255, 216 / 255, 1),
        # B
        rgb(1,       244 / 255, 243 / 255),
        # A
        rgb(1,       229 / 255, 207 / 255),
        # F
        rgb(1,       217 / 255, 178 / 255),
        # G
        rgb(1,       199 / 255, 142 / 255),
        # K
        rgb(1,       166 / 255, 81 / 255)
    )  # M
    #if (F) {
    
    par(family='Helvetica LT Std Light')#'Droid')#
    for (ii in 1:length(DFs)) {
        DF <- DFs[[ii]]
        name <- as.numeric(strsplit(names(DFs)[ii], '.dat')[[1]][1])
        with(DF, lines(Teff, L, lwd=1.33, 
            col=ifelse(name>0, orange, ifelse(name<0, blue, 1)),
            lty=1))
    }
    for (ii in 1:length(DFs)) {
        DF <- DFs[[ii]]
        name <- as.numeric(strsplit(names(DFs)[ii], '.dat')[[1]][1])
        with(DF, points(Teff[1], L[1], pch=20, cex=0.1)) 
        with(DF, points(Teff[length(DF)], L[length(DF)], pch=1, cex=0.5)) 
        if (name == 0.05)
            with(DF, text(Teff[1]-5, L[1]-0.004, 
                pos=2, labels=bquote(.(name)), cex=0.8*text.cex, col=orange))
        if (name == 0)
            with(DF, text(Teff[1], L[1], 
                pos=1, labels=expression(beta==0), cex=0.8*text.cex))
        if (name == -0.05)# || name == 0 || name == 0.05) {
            with(DF, text(Teff[1]+5, L[1]-0.004, 
                pos=4, labels=bquote(.(name)), cex=0.8*text.cex, col=blue))
    }
    par(family=font)
    
    with(DFs[['0.00.dat']], lines(Teff, L, lwd=1.66))
    
    par(xpd = NA)
    rect(xlim[2],
        ylim[1] * 0.1,
        xlim[2] - 100,
        ylim[2] * 10,
        col = 'white',
        border = NA)
    rect(xlim[1] * 1.1,
        ylim[1],
        xlim[2] * 0.9,
        ylim[1] * 0.9,
        col = 'white',
        border = NA)
    par(xpd = F)
    
    
    par(mgp = mgp + c(0, 0.43, 0))
    magaxis(
        2,
        tcl = tcl,
        mgp = mgp + c(0, 0.45, 0),
        cex.axis = text.cex,
        family = font,
        las = 1,
        majorn = 3,
        labels = T,
        lwd.ticks = par()$lwd
    )
    magaxis(
        1,
        tcl = tcl,
        mgp = mgp + c(0, 0.35, 0),
        cex.axis = text.cex,
        family = font,
        las = 1,
        majorn = 3,
        labels = T,
        lwd.ticks = par()$lwd
    )
    
    box(lwd = par()$lwd)
    
    mtext(
        expression(Luminosity~L / L["solar"]),
        2,
        2,
        outer = F,
        las = 0,
        cex = 1.15*text.cex
    )
    mtext(expression(Effective~temperature~T["eff"] / K),
        1,
        2.2,
        outer = F,
        cex = 1.15*text.cex)
    
    mtext(
        expression("Spectral type"),
        side = 3,
        line = 1.3,
        cex = 1.15*text.cex
    )
    
    spectral.labs <- c("O", "B", "A", "F", "G", "K", "M")
    selector <- 1:(length(spectral.divs))
    spectral.Teffs <- sapply(selector,
        function(ii) {
            div <- spectral.divs[ii]
            if (div >= xlim[1])
                return(Inf)
            if (div < xlim[2])
                div <- xlim[2]
            if (ii == 1)
                return((xlim[1] + div) / 2)
            prev <- spectral.divs[ii - 1]
            if (prev > xlim[1])
                prev <- xlim[1]
            (div + prev) / 2
        })
    axis(
        3,
        at = spectral.divs,
        tcl = tcl,
        labels = F,
        cex.axis = text.cex,
        lwd.ticks = par()$lwd
    )
    par(mgp = mgp + c(0, -0.05, 0))
    axis(
        3,
        at = spectral.Teffs,
        labels = spectral.labs[selector],
        cex.axis = text.cex,
        tcl = 0
    )
}

make_plots(
    plot_HR,
    'HR',
    #paper_pdf_width = 4.17309 * (1.2 * 0.8),
    #paper_pdf_height = 4.17309 * 0.8,
    paper_pdf_width = 4.17309 * (1.2 * 0.95),
    paper_pdf_height = 4.17309 * 0.95,
    cex.paper = 1.12,
    use.cairo = T,
    thin=F, make_png=F, short=F, slides=F, 
    font = 'grace'#'Palatino Linotype'
)

