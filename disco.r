library(rpart)
library(ggplot2)
library(MASS)
library(lattice)

setwd('/afs/cs.wisc.edu/u/b/r/brook/private/STAT 998/disco_farms/code')
disco = read.csv('../data/disco.csv')

disco$soil_moisture = as.numeric(disco$soil_moisture)
disco$runoff_occured = ifelse(disco$runoff>0,T,F)
disco$r = ifelse(disco$runoff_occured,1,-1)
disco$absorbed = disco$precip - disco$runoff
disco$a60 = disco$i60 - disco$runoff
disco$a30 = disco$i30 - disco$runoff
disco$a15 = disco$i15 - disco$runoff
disco$discrete_precip = disco$precip - disco$precip%%0.5
disco$rounded_precip = round(disco$precip * 2) / 2

koepke = disco[disco$farm=='Koepke',]
pagel = disco[disco$farm=='Pagel',]
pioneer = disco[disco$farm=='Pioneer',]
riecher = disco[disco$farm=='Riecher',]
saxon = disco[disco$farm=='Saxon',]


stump_control <- rpart.control(maxdepth=1, maxsurrogate=0, usesurrogate=0, maxcompete=1, cp=0, xval=0) 
with( disco, rpart(runoff_occured~soil_moisture, control=stump_control, weights=ifelse(r==1,w,1)) ) 
stump1 = rpart(runoff_occured~soil_moisture, data=disco, control=stump_control)


xx=c(15,50)
yy=c(0,3.5)


#####
#Finding the breakpoint (Minimizing sum-of-squared errors)
#####
piecewise <- function(th, x, y) { # conditional minimum SSQ given theta
        X <- cbind(x, pmax(0, x - th))
        sum( abs(lsfit(X, y)$resid**2) )
}

#####
#Finding the breakpoint (Minimizing absolute errors)
#####
piecewise_abs <- function(th, x, y) { # conditional minimum SSQ given theta
        X <- cbind(x, pmax(0, x - th))
        sum( abs(lsfit(X, y)$resid) )
}


#####
#Finding the breakpoint (minimize average standard error)
#####
piecewise_std <- function(th, x, y) { # conditional minimum SSQ given theta
        X <- cbind(x, pmax(0, x - th))
        sd( lsfit(X, y)$resid )
}

xs = c(31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 42.5, 43.5, 44.5, 45.5)
th <- optimize(piecewise_abs, range(pagel$soil_moisture), x=pagel$soil_moisture, y=pagel$runoff_coef)$minimum
fm <- lm(y ~ x + pmax(0, x - th))

#####
#Finding the breakpoint 2-dimensional, use absolute errors
#####
piecewise_2dim_abs <- function(x, y, z, th, C) { # conditional minimum SSQ given theta
        X <- cbind(x, y, pmax(0, x + C*y - th))
        sum( abs(lsfit(X, z)$resid) )
}

piecewise_SM_abs <- function(th, C, x, y, z) { # conditional minimum SSQ given theta
        X <- cbind(x, pmax(0, x + C*y - th))
        sum( abs(lsfit(X, z)$resid) )
}

piecewise_precip_abs <- function(C, x, y, z) { # conditional minimum SSQ given theta
        optimize(piecewise_SM_abs, range(x), C=C, x=x, y=y, z=z)$minimum
}


th = with(disco, optimize( piecewise_precip_abs, c(0,10), x=soil_moisture, y=precip, z=runoff_coef))$objective

with(disco, optim(1, piecewise_SM_abs, C=1, x=soil_moisture, y=precip, z=runoff_coef))


pdf('../figures/moisture_boxplot.pdf')
boxplot(disco$soil_moisture~factor(disco$runoff_occured, labels=c('no runoff','runoff') ), frame.plot=F, ylab='soil moisture copntent (%)', boxwex=0.5, ylim=xx)
dev.off()

pdf('../figures/precip_boxplot.pdf')
boxplot(disco$precip~factor(disco$runoff_occured, labels=c('no runoff','runoff') ), frame.plot=F, ylab='rainfall (inches)', boxwex=0.5)
dev.off()

pdf('../figures/rainfall_v_moisture.pdf')
plot(x=jitter(disco$soil_moisture), y=disco$precip, bty='n', cex=0.4, xlab='soil moisture content (%)', ylab='rainfall (inches)', xlim=xx)
dev.off()


lik <- function(x0, n0, x1, n1) {
    -2*(x0*log(x0/n0) + (n0-x0)*log((n0-x0)/n0) + x1*log(x1/n1) + (n1-x1)*log((n1-x1)/n1))
}


fm1=rpart(runoff_occured ~ precip + soil_moisture, data=disco, method='class')
svm1 = ksvm(runoff_occured ~ precip + soil_moisture, data=disco, kernel='vanilladot', type='C-svc')
lda1 = lda(runoff_occured ~ precip + soil_moisture, data=disco)

#LDA-based classifier for each of the farms
w = min(dim(disco)[1]/sum(disco$runoff_occured) - 1, 3)
w_k = min(dim(koepke)[1]/sum(koepke$runoff_occured) - 1, 3)
w_pa = min(dim(pagel)[1]/sum(pagel$runoff_occured) - 1, 3)
w_pi = min(dim(pioneer)[1]/sum(pioneer$runoff_occured) - 1, 3)
w_r = min(dim(riecher)[1]/sum(riecher$runoff_occured) - 1, 3)
w_s = min(dim(saxon)[1]/sum(saxon$runoff_occured) - 1, 3)

fm = lm(r ~ precip + soil_moisture, data=disco, weights=ifelse(runoff_occured,w,1) )
fm_k= lm(r ~ precip + soil_moisture, data=koepke, weights=ifelse(runoff_occured,w_k,1) )
fm_pa= lm(r ~ precip + soil_moisture, data=pagel, weights=ifelse(runoff_occured,w_pa,1))
fm_pi= lm(r ~ precip + soil_moisture, data=pioneer, weights=ifelse(runoff_occured,w_pi,1))
fm_r= lm(r ~ precip + soil_moisture, data=riecher, weights=ifelse(runoff_occured,w_r,1))
fm_s= lm(r ~ precip + soil_moisture, data=saxon, weights=ifelse(runoff_occured,w_s,1))

a = -fm$coefficients[1]/fm$coefficients[2]
a_k = -fm_k$coefficients[1]/fm_k$coefficients[2]
a_pa = -fm_pa$coefficients[1]/fm_pa$coefficients[2]
a_pi = -fm_pi$coefficients[1]/fm_pi$coefficients[2]
a_r = -fm_r$coefficients[1]/fm_r$coefficients[2]
a_s = -fm_s$coefficients[1]/fm_s$coefficients[2]

b = -fm$coefficients[3]/fm$coefficients[2]
b_k = -fm_k$coefficients[3]/fm_k$coefficients[2]
b_pa = -fm_pa$coefficients[3]/fm_pa$coefficients[2]
b_pi = -fm_pi$coefficients[3]/fm_pi$coefficients[2]
b_r = -fm_r$coefficients[3]/fm_r$coefficients[2]
b_s = -fm_s$coefficients[3]/fm_s$coefficients[2]

farms = c("Koepke", "Pagel", "Pioneer", "Riecher", "Saxon")
aa = c(a_k, a_pa, a_pi, a_r, a_s)
bb = c(b_k, b_pa, b_pi, b_r, b_s)
coefs = data.frame(farms, aa, bb)

#Rpart for "breakpoints"
stump_control <- rpart.control(maxdepth=1, maxsurrogate=0, usesurrogate=0, maxcompete=1, cp=0, xval=0) 
stump1 = rpart(runoff_occured~soil_moisture, data=disco, control=stump_control, weights=ifelse(runoff_occured,w,1))
stump_k = rpart(runoff_occured~soil_moisture, data=koepke, control=stump_control)
stump_pa = rpart(runoff_occured~soil_moisture, data=pagel, control=stump_control)
stump_pi = rpart(runoff_occured~soil_moisture, data=pioneer, control=stump_control)
stump_r = rpart(runoff_occured~soil_moisture, data=riecher, control=stump_control)
stump_s = rpart(runoff_occured~soil_moisture, data=saxon, control=stump_control)

stump_k_weighted = rpart(runoff_occured~soil_moisture, data=koepke, control=stump_control, weights=ifelse(runoff_occured,w,1))
stump_pa_weighted = rpart(runoff_occured~soil_moisture, data=pagel, control=stump_control, weights=ifelse(runoff_occured,w,1))
stump_pi_weighted = rpart(runoff_occured~soil_moisture, data=pioneer, control=stump_control, weights=ifelse(runoff_occured,w,1))
stump_r_weighted = rpart(runoff_occured~soil_moisture, data=riecher, control=stump_control, weights=ifelse(runoff_occured,w,1))
stump_s_weighted = rpart(runoff_occured~soil_moisture, data=saxon, control=stump_control, weights=ifelse(runoff_occured,w,1))


xp <- seq(xx[1], xx[2], length = 50)
yp <- seq(yy[1], yy[2], length = 50)
np <- length(xp)
pt <- expand.grid(soil_moisture = xp, precip = yp)

predictions <- predict(fm1, pt)
zp.rp <- predictions[,2]

pdf('../figures/aggregate_tree.pdf')
plot(disco$soil_moisture, disco$precip, xlab = "soil moisture (%)", ylab = "rain (inches)", pch=as.integer(disco$farm), col=as.integer(disco$runoff_occured)+1, bty='n', xlim=xx, ylim=yy, cex=0.5)
contour(xp, yp, matrix(zp.rp, np), add = T, levels=c(0.1, 0.8), lwd=2, col=c('purple', 'orange'))
abline(a=a, b=b, lty=2, col='red', lwd=2)
legend(c('Runoff w/ confidence >0.1', 'Runoff w/ confidence >0.8', 'Linear decision rule'), lty=c(1,1,2), lwd=c(2,2,2), x='topleft', bty='n', col=c('purple', 'orange', 'red') )
dev.off()


pdf('../figures/farm_by_moisture_boxplot.pdf')
qplot(factor(disco$runoff_occured, labels=c('none','runoff') ), soil_moisture, data=disco, geom="boxplot") + facet_grid(farm~.) + coord_flip()+theme_bw() + ylab('\nsoil moisture content (%)') + xlab('')
dev.off()

pdf('../figures/farm_by_precip_boxplot.pdf')
qplot(factor(disco$runoff_occured, labels=c('none','runoff') ), precip, data=disco, geom="boxplot") + facet_grid(farm~.) + coord_flip()+theme_bw() + ylab('\nrainfall (inches)') + xlab('')
dev.off()

pdf('../figures/scatter_by_farm.pdf')
(p <- qplot(x=soil_moisture, y=precip, data=disco, geom="point", shape=ifelse(runoff_occured,17,1), color=ifelse(runoff_occured,1,0)) + facet_wrap(~farm, 3,2) + ylab('rainfall (inches)\n') + xlab('\nsoil moisture content (%)') + opts(legend.position='none') )
p + geom_abline(data=coefs, aes(intercept=a, slope=b))

dev.off()


pdf('../figures/aggregate_linear.pdf')
plot(disco$soil_moisture, disco$precip, xlab = "soil moisture (%)", ylab = "rain (inches)", pch=as.integer(disco$farm), col=ifelse(disco$runoff_occured,'red', 'blue'), bty='n', xlim=xx, ylim=yy, cex=0.5)
abline(a=a, b=b, lty=2)
dev.off()



layout(matrix(1:6,2,3))
xx=c(15,50)
yy=c(0,3.5)

#Koepke
k1 = koepke[koepke$runoff_occured==1,]
k2 = koepke[koepke$runoff_occured==0,]
plot(precip~soil_moisture, data=k1, pch=17, col='red', xlim=xx, ylim=yy, bty='n', cex.lab=1.5)
par(new=T)
plot(precip~soil_moisture, data=k2, pch=1, col='blue', xlim=xx, ylim=yy, bty='n', ann=F, xaxt='n', yaxt='n')
abline(a=a_k, b=b_k, lty=2)
title("Koepke", cex=1.5)

#Pagel
pa1 = pagel[pagel$runoff_occured==1,]
pa2 = pagel[pagel$runoff_occured==0,]
plot(precip~soil_moisture, data=pa1, pch=17, col='red', xlim=xx, ylim=yy, bty='n', cex.lab=1.5)
par(new=T)
plot(precip~soil_moisture, data=pa2, pch=1, col='blue', xlim=xx, ylim=yy, bty='n', ann=F, xaxt='n', yaxt='n')
abline(a=a_pa, b=b_pa, lty=2)
title("Pagel", cex=1.5)

#Pioneer
pi1 = pioneer[pioneer$runoff_occured==1,]
pi2 = pioneer[pioneer$runoff_occured==0,]
plot(precip~soil_moisture, data=pi1, pch=17, col='red', xlim=xx, ylim=yy, bty='n', cex.lab=1.5)
par(new=T)
plot(precip~soil_moisture, data=pi2, pch=1, col='blue', xlim=xx, ylim=yy, bty='n', ann=F, xaxt='n', yaxt='n')
abline(a=a_pi, b=b_pi, lty=2)
title("Pioneer", cex=1.5)

#Riecher
r1 = riecher[riecher$runoff_occured==1,]
r2 = riecher[riecher$runoff_occured==0,]
plot(precip~soil_moisture, data=r1, pch=17, col='red', xlim=xx, ylim=yy, bty='n', cex.lab=1.5)
par(new=T)
plot(precip~soil_moisture, data=r2, pch=1, col='blue', xlim=xx, ylim=yy, bty='n', ann=F, xaxt='n', yaxt='n')
abline(a=a_r, b=b_r, lty=2)
title("Riecher", cex=1.5)

#Saxon
s1 = saxon[saxon$runoff_occured==1,]
s2 = saxon[saxon$runoff_occured==0,]
plot(precip~soil_moisture, data=s1, pch=17, col='red', xlim=xx, ylim=yy, bty='n', cex.lab=1.5)
par(new=T)
plot(precip~soil_moisture, data=s2, pch=1, col='blue', xlim=xx, ylim=yy, bty='n', ann=F, xaxt='n', yaxt='n')
abline(a=a_s, b=b_s, lty=2)
title("Saxon", cex=1.5)

disco = within(disco, fitted_class <- ifelse(a+b*soil_moisture>precip,-1,1))
koepke = within(koepke, fitted_class <- ifelse(a_k+b_k*soil_moisture>precip,-1,1))
pagel = within(pagel, fitted_class <- ifelse(a_pa+b_pa*soil_moisture>precip,-1,1))
pioneer = within(pagel, fitted_class <- ifelse(a_pi+b_pi*soil_moisture>precip,-1,1))
riecher = within(riecher, fitted_class <- ifelse(a_r+b_r*soil_moisture>precip,-1,1))
saxon = within(saxon, fitted_class <- ifelse(a_s+b_s*soil_moisture>precip,-1,1))

fme = lmer(r ~ soil_moisture + precip+(precip+soil_moisture|farm) , data=disco, weights=ifelse(runoff_occured,w,1) )
fme2 = lm(r ~ soil_moisture+precip+soil_moisture:farm + precip:farm, data=disco, weights=ifelse(runoff_occured,w,1) )
fme3 = lmer(r ~ (soil_moisture + precip|farm), data=disco, weights=ifelse(runoff_occured,w,1) )




#Cross-validation (without farm factor, soil moisture only)
#Take CV training set
#randomize
f = 10
B = 200
f_pos = vector()
f_neg = vector()
t_pos = vector()
t_neg = vector()

for( j in 1:B ) {
    indexes = 1:length(disco[,1])
    r = sample(indexes, replace=F)
    pp = vector()
    tt = vector()
    
    folds = list()
    for(i in 1:f) {
        q_l = 1-i/f
        q_h = 1-(i-1)/f
        fold = indexes[r<quantile(indexes, q_h) & r>=quantile(indexes, q_l)]
    
        train = disco[-fold,]
        test = disco[fold,]

        threshold = with( train, optimize(piecewise_abs, range(soil_moisture), x=soil_moisture, y=runoff_coef)$minimum )

        pp = c(pp, ifelse(test$soil_moisture>=threshold,1,0))
        tt = c(tt, test$runoff_occured)
    }
    
    f_pos = c(f_pos, sum(pp==1&tt==0))
    f_neg = c(f_neg, sum(pp==0&tt==1))
    t_pos = c(t_pos, sum(pp==1&tt==1))
    t_neg = c(t_neg, sum(pp==0&tt==0))

}
no_farm_no_precip = data.frame(f_pos, f_neg, t_pos, t_neg)


#Cross-validation (with farm factor, soil moisture only)
#Take CV training set
#randomize
f = 10
B = 200
f_pos = vector()
f_neg = vector()
t_pos = vector()
t_neg = vector()

for( j in 1:B ) {
    indexes = 1:length(disco[,1])
    r = sample(indexes, replace=F)
    pp = vector()
    tt = vector()
    
    folds = list()
    for(i in 1:f) {
        q_l = 1-i/f
        q_h = 1-(i-1)/f
        fold = indexes[r<quantile(indexes, q_h) & r>=quantile(indexes, q_l)]
    
        train = disco[-fold,]
        test = disco[fold,]
    
        w = min( c(dim(train)[1]/sum(train$runoff_occured)-1, 3) )
    
        for(farm in farms) {
            train_farm = train[train$farm==farm,]
            threshold_farm = with( train_farm, optimize(piecewise_abs, range(soil_moisture), x=soil_moisture, y=runoff_coef)$minimum )
            indx = which(test$farm==farm)
            
            pp = c(pp, ifelse(test[indx,'soil_moisture']>=threshold,1,0))
            tt = c(tt, test[indx,'runoff_occured'])

        }
    }
    
    f_pos = c(f_pos, sum(pp==1&tt==0))
    f_neg = c(f_neg, sum(pp==0&tt==1))
    t_pos = c(t_pos, sum(pp==1&tt==1))
    t_neg = c(t_neg, sum(pp==0&tt==0))

}
farm_no_precip = data.frame(f_pos, f_neg, t_pos, t_neg)











#Cross-validation (without farm factor - soil moisture and precip)
#Take CV training set
#randomize
f = 10
B = 200
lr_f_pos = vector()
lr_f_neg = vector()
lr_t_pos = vector()
lr_t_neg = vector()

lda_f_pos = vector()
lda_f_neg = vector()
lda_t_pos = vector()
lda_t_neg = vector()

for( j in 1:B ) {
    indexes = 1:length(disco[,1])
    r = sample(indexes, replace=F)
    lr_pp = vector()
    lda_pp = vector()
    tt = vector()
    
    folds = list()
    for(i in 1:f) {
        q_l = 1-i/f
        q_h = 1-(i-1)/f
        fold = indexes[r<quantile(indexes, q_h) & r>=quantile(indexes, q_l)]
    
        train = disco[-fold,]
        test = disco[fold,]
    
        w = min( c(dim(train)[1]/sum(train$runoff_occured)-1, 3) )
    
        #lr_model = with(train, glm(runoff_occured~soil_moisture + precip, family='binomial', weights=ifelse(runoff_occured,2,1)))
        lda_model = with(train, lm(r~ precip + soil_moisture, weights=ifelse(runoff_occured,w,1)))

        intercept = -lda_model$coefficients[1]/lda_model$coefficients[2]
        slope = -lda_model$coefficients[3]/lda_model$coefficients[2]

        #lr_prediction = predict(lr_model, test)
        #lda_prediction = predict(lda_model, test)

        #lr_pp = c(lr_pp, ifelse(lr_prediction>=0,1,0))
        #lda_pp = c(lda_pp, ifelse(test$discrete_precip - (intercept+slope*test$soil_moisture)>=0,1,0) )
        lda_pp = c(lda_pp, ifelse(test$rounded_precip - (intercept+slope*test$soil_moisture)>=0,1,0) )
        tt = c(tt, test$runoff_occured)
    }
    
    lr_f_pos = c(lr_f_pos, sum(lr_pp==1&tt==0))
    lr_f_neg = c(lr_f_neg, sum(lr_pp==0&tt==1))
    lr_t_pos = c(lr_t_pos, sum(lr_pp==1&tt==1))
    lr_t_neg = c(lr_t_neg, sum(lr_pp==0&tt==0))

    lda_f_pos = c(lda_f_pos, sum(lda_pp==1&tt==0))
    lda_f_neg = c(lda_f_neg, sum(lda_pp==0&tt==1))
    lda_t_pos = c(lda_t_pos, sum(lda_pp==1&tt==1))
    lda_t_neg = c(lda_t_neg, sum(lda_pp==0&tt==0))
}
no_farm = data.frame(lr_f_pos, lr_f_neg, lr_t_pos, lr_t_neg, lda_f_pos, lda_f_neg, lda_t_pos, lda_t_neg)


#Cross-validation (with farm factor - soil moisture and precip)
#Take CV training set
#randomize
farms = levels(disco$farm)
f = 10
B = 200
lr_f_pos = vector()
lr_f_neg = vector()
lr_t_pos = vector()
lr_t_neg = vector()

lda_f_pos = vector()
lda_f_neg = vector()
lda_t_pos = vector()
lda_t_neg = vector()

for( j in 1:B ) {
    indexes = 1:length(disco[,1])
    r = sample(indexes, replace=F)
    lr_pp = vector()
    lda_pp = vector()
    tt = vector()
    
    folds = list()
    for(i in 1:f) {
        q_l = 1-i/f
        q_h = 1-(i-1)/f
        fold = indexes[r<quantile(indexes, q_h) & r>=quantile(indexes, q_l)]
    
        train = disco[-fold,]
        test = disco[fold,]
    
        w = min( c(dim(train)[1]/sum(train$runoff_occured)-1, 3) )
    
        #lr_model = with(train, glm(runoff_occured~soil_moisture:farm + precip:farm, family='binomial', weights=ifelse(runoff_occured,2,1)))


        for(farm in farms) {
            train_farm = train[train$farm==farm,]
            indx = which(test$farm==farm)
                
            lda_model = with(train_farm, lm(r~ precip + soil_moisture, weights=ifelse(runoff_occured,w,1)))
    
            intercept = -lda_model$coefficients[1]/lda_model$coefficients[2]
            slope = -lda_model$coefficients[3]/lda_model$coefficients[2]
    
            #lda_pp = c(lda_pp, ifelse(test[indx,'discrete_precip'] - (intercept+slope*test[indx,'soil_moisture'])>=0,1,0) )
            lda_pp = c(lda_pp, ifelse(test[indx,'rounded_precip'] - (intercept+slope*test[indx,'soil_moisture'])>=0,1,0) )
            tt = c(tt, test[indx,'runoff_occured'])

        }


        #lr_prediction = predict(lr_model, test)

        #lr_pp = c(lr_pp, ifelse(lr_prediction>=0,1,0))
    }
    
    #lr_f_pos = c(lr_f_pos, sum(lr_pp==1&tt==0))
    #lr_f_neg = c(lr_f_neg, sum(lr_pp==0&tt==1))
    #lr_t_pos = c(lr_t_pos, sum(lr_pp==1&tt==1))
    #lr_t_neg = c(lr_t_neg, sum(lr_pp==0&tt==0))

    lda_f_pos = c(lda_f_pos, sum(lda_pp==1&tt==0))
    lda_f_neg = c(lda_f_neg, sum(lda_pp==0&tt==1))
    lda_t_pos = c(lda_t_pos, sum(lda_pp==1&tt==1))
    lda_t_neg = c(lda_t_neg, sum(lda_pp==0&tt==0))
}
#farm = data.frame(lr_f_pos, lr_f_neg, lr_t_pos, lr_t_neg, lda_f_pos, lda_f_neg, lda_t_pos, lda_t_neg)
farm = data.frame(f_pos=lda_f_pos, f_neg=lda_f_neg, t_pos=lda_t_pos, t_neg=lda_t_neg)




#Compute the deviance for the chi-squared test of a significant contribution by farm:

predicted = ifelse(disco$precip - (a+b*disco$soil_moisture)>0,1,0)
null_dev = lik(sum(disco$r==1 & predicted==0), sum(predicted==0), sum(disco$r==1 & predicted==1), sum(predicted==1))

dev=0
for(f in farms) {
    indx = which(disco$farm==f)
    dev = dev + lik(sum(disco[indx,'r']==1 & predicted[indx]==0), sum(predicted[indx]==0), sum(disco[indx,'r']==1 & predicted[indx]==1), sum(predicted[indx]==1))
}

xx=c(15,50)
yy=c(0,3.5)

plot(disco$soil_moisture, disco$precip, xlab = "soil moisture (%)", ylab = "rain (inches)", pch=ifelse(disco$runoff_occured,17,1), col=ifelse(disco$runoff_occured,'black', 'grey50'), bty='n', xlim=xx, ylim=yy, cex=0.6)
legend(legend=c('Runoff event', 'Non-runoff event'), pch=c(17,1), col=c('black', 'gray50'), x='topright', bty='n')