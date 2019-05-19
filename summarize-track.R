#### Obtain model properties from evolutionary tracks 
#### Author: Earl Bellinger ( bellinger@phys.au.dk ) 
#### Stellar Astrophysics Centre Aarhus 

source('/home/earl/asteroseismology/scripts/seismology.R')
source('/home/earl/asteroseismology/scripts/utils.R')

Z_div_X_solar = 0.02293 # GS98
Rsol = 6.9599 * 10**10

directory <- commandArgs(1)[1]
print(directory)

trackfile <- file.path(directory, 'track')
params.DF <- read.table(trackfile, header=1)
params.DF <- params.DF[,-which(colnames(params.DF)=='age')]

hstry <- read.table(file.path(directory, 'LOGS', 'history.data'),
    col.names=c('Model', 'M', 'age', 'radius', 'Teff', 'L', 'X_c', 'qc'),
    header=FALSE, skip=4)

DF <- do.call(plyr:::rbind.fill, Map(function(model) {
    fgong <- read.table(file.path(directory, 'LOGS', 
            paste0('profile', model, '.FGONG.dat')),
        header=TRUE)[1,]
    
    hst <- hstry[model,]
    
    obs.DF <- NULL
    obs.DF['age'] <- hst$age
    obs.DF['Teff'] <- hst$Teff
    obs.DF['L'] <- hst$L
    obs.DF['radius'] <- hst$radius / Rsol
    obs.DF['Fe_H'] <- with(fgong, log10(Z/X/Z_div_X_solar))
    obs.DF['nu_max'] <- nu_max_scaling(hst$M, obs.DF['radius'], hst$Teff,
        Teff_sun=5777.54)
    
    freqs <- parse_freqs(file.path(directory, 'LOGS', 
            paste0('profile', model, '-freqs.dat')), 
        new.adipls=T)
    seis.DF <- seismology(freqs, nu_max=obs.DF['nu_max'], 
        all_freqs=T, all_ratios=T)
    
    cbind(params.DF, t(as.data.frame(obs.DF)), seis.DF)
}, model=hstry$Model))

write.table(DF, paste0(directory, '.dat'), quote=FALSE, sep='\t', 
    row.names=FALSE)
