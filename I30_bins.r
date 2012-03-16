# 1) Put each farm's data into two bins: one with soil moisture above the breakpoint, the other below.
# 2) Find the I30 breakpoint within each of the soil moisture bins.
# 3) Make two plots of RC versus I30: one for the low SM bin, one for the high SM bin.
#This all depends on the R object 'breakpoints' being a list of the breakpoints by site.
#That list of breakpoints is created by the code you've already copied.


#Find the I30 breakpoint within each of the bins:
binned_I30_breakpoints2 = list()
limits2 = list()

for(site in farms)
{
    limits2[[site]] = c(-Inf, breakpoints[[site]], Inf)
    site_data = data[data$site==site,]
    bp = vector()
    cat(paste("Intensity breakpoints at ", site, " when binned by soil moisture:\n", sep=""))
    
    for(i in 1:2)
    {
        bin = site_data[site_data$sm>=limits2[[site]][i] & site_data$sm<limits2[[site]][i+1] & !is.na(site_data$I30),]
        lower_lim = sort( bin$I30 )[4]
        upper_lim = sort(bin$I30, decreasing=T )[4]
        th <- optimize(piecewise, c(upper_lim, lower_lim), x=bin$I30, y=bin$rc)$minimum
        cat( paste(limits2[[site]][i], " <= SM < ", limits2[[site]][i+1], ": ", round(th,1), "\n", sep=""))
        bp = c(bp, th)
    }
    binned_I30_breakpoints2[[site]] = bp
    cat("\n")
}

#Plot each bin separately, with its I30 breakpoint.
layout(t(matrix(1:6,2,3)))

for(site in farms)
{
    titles = c(paste("SM < ", limits2[[site]][2], "%", sep=""), paste(limits2[[site]][2], "% <= SM", sep=""))
    site_data = data[data$site==site,]
    
    for(i in 1:2)
    {
        bin = site_data[site_data$sm>=limits2[[site]][i] & site_data$sm<limits2[[site]][i+1] & !is.na(site_data$I30),]
        piecewise_plot(x=bin$I30, y=bin$rc, th=binned_I30_breakpoints2[[site]][i], ylim=range(data$rc, na.rm=T),
                         xlim=range(data$I30, na.rm=T), xlab="I30", ylab="runoff coefficient", title=titles[i])
    }
}