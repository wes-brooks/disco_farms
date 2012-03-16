# 1) Put each farm's data into two bins: one with soil moisture above the breakpoint, the other below.
# 2) Find the precip breakpoint within each of the soil moisture bins.
# 3) Make two plots of RC versus precipitation: one for the low SM bin, one for the high SM bin.
#This all depends on the R object 'breakpoints' being a list of the breakpoints by site.
#That list of breakpoints is created by the code you've already copied.


#Find the precip breakpoint within each of the bins:
binned_precip_breakpoints2 = list()
limits2 = list()

for(site in farms)
{
    limits2[[site]] = c(-Inf, breakpoints[[site]], Inf)
    site_data = data[data$site==site,]
    bp = vector()
    cat(paste("Precipitation breakpoints at ", site, " when binned by soil moisture:\n", sep=""))
        
    for(i in 1:2)
    {
        bin = site_data[site_data$sm>=limits2[[site]][i] & site_data$sm<limits2[[site]][i+1] & !is.na(site_data$prec),]
        lower_lim = sort( bin$prec )[4]
        upper_lim = sort(bin$prec, decreasing=T )[4]
        th <- optimize(piecewise, c(upper_lim, lower_lim), x=bin$prec, y=bin$rc)$minimum
        cat( paste(limits2[[site]][i], " <= SM < ", limits2[[site]][i+1], ": ", round(th,2), "\n", sep=""))
        bp = c(bp, th)
    }
    binned_precip_breakpoints2[[site]] = bp
    cat("\n")
}


#Plot each bin separately, with its precip breakpoint.
layout(t(matrix(1:6,2,3)))

for(site in farms)
{    
    titles = c(paste("SM < ", limits2[[site]][2], "%", sep=""), paste(limits2[[site]][2], "% <= SM", sep=""))
    site_data = data[data$site==site,]
    
    for(i in 1:2)
    {
        bin = site_data[site_data$sm>=limits2[[site]][i] & site_data$sm<limits2[[site]][i+1] & !is.na(site_data$prec),]
        piecewise_plot(x=bin$prec, y=bin$rc, th=binned_precip_breakpoints2[[site]][i], ylim=range(data$rc, na.rm=T),
                         xlim=range(data$prec, na.rm=T), xlab="precip", ylab="runoff coefficient", title=titles[i])
    }
}