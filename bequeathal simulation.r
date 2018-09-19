############################################
#
# Individual-based R simulation code
# to accompany Clarke, McElreath, Mabry, McEachern.
#
# To use: Enter the functions progtext() and beqsim() into R's memory,
#         then 

progtext <- function( now , min=0 , max=100 , start.time , prefix="Task" , interval=100 ) {
    if ( floor(now/interval) != now/interval ) return()
    dst <- (proc.time() - start.time)[3] # elapsed time in seconds
    pdone <- now/(max-min)
    dh <- floor(dst/(60*60))
    dm <- floor((dst - dh*60*60)/60)
    ds <- floor(dst - dh*360 - dm*60)
    rate <- (now-min)/dst # rate as tasks per second
    etr <- (max-now) / rate # estimated time remaining as # seconds
    etrh <- floor(etr/(60*60))
    etrm <- floor((etr - etrh*60*60)/60)
    etrs <- floor(etr - etrh*60*60 - etrm*60)
    cat( paste( "\r" , prefix , " " , now , "/" , max , " done (" , format(pdone*100,digits=1) , "%) | Elapsed: " , dh , "h " , dm , "m " , ds , "s" , " | Rate: " , format(rate,digits=2) , " o/s" , " | ETA: " , etrh , "h " , etrm , "m " , etrs , "s    " , collapse="" , sep="" ) )
    cat( "\r" )
}

# bequeathal sim
# pa: prob adult survives dispersal
# pj: prob juvenile survives dispersal
# va > 1: competitive ability of adult (Ca in paper)
# ma: adult mortality rate (1-sa)
# mj: juvenile mortality rate (1-sj)
# n: number of sites
# s: heritable prob that adult stays at site, forcing offspring to disperse
beqsim <- function( tmax=100 , n=10 , pa=0.75 , pj=0.5 , va=2 , ma=0 , mj=0 , sstart=0.9 , mu=0.01 , sa , sj , da , dj , local.dispersal=FALSE , progress=TRUE , sex=TRUE , discrete=FALSE ) {
    st <- proc.time()
    slotspersite <- 100
    shist <- rep(0,tmax)
    shistfull <- matrix(0,nrow=tmax,ncol=n)
    agehist <- matrix(0,nrow=tmax,ncol=n)
    occupancyhist <- rep(0,tmax) # contains history of proportion of sites occupied
    
    # check params
    if ( !missing(sa) ) ma <- 1 - sa
    if ( !missing(sj) ) mj <- 1 - sj
    if ( !missing(da) ) pa <- da
    if ( !missing(dj) ) pj <- dj
    
    # utility function for finding first empty slot at a site
    findslot <- function( slots ) {
        aslot <- 0
        for ( i in 1:length(slots) ) {
            if ( is.na(slots[i]) ) {
                aslot <- i
                break
            }
        }
        aslot
    }
    countfilled <- function( slots ) {
        num <- 0
        for ( i in 1:length(slots) ) {
            if ( !is.na(slots[i]) ) {
                num <- num + 1 
            } else { 
                break 
            }
        }
        num
    }
    fclip <- function(x) ifelse( x > 1 , 1 , ifelse( x < 0 , 0 , x ) )
    # init population
    if ( length(sstart) > 1 ) {
        # passed vector of initial site states
        s <- sstart
    } else {
        # passed initial frequency of Stay
        s <- fclip( rnorm( n , mean=sstart , sd=mu*10 ) )
        if ( discrete==TRUE ) {
            # only allow s=0 and s=1 as genotypes
            s <- sample( c(0,1) , size=n , replace=TRUE , prob=c( 1-sstart , sstart ) )
        }
    }
    
    ages <- rep( 1 , n ) # all start at age 1 (adults)
    # loop over generations
    for ( i in 1:tmax ) {
        # check for collapsed population
        n_adults <- sum( !is.na(s) )
        # record keeping
        shist[i] <- mean( s , na.rm=TRUE )
        shistfull[i,] <- s
        agehist[i,] <- ages
        occupancyhist[i] <- sum(!is.na(s))/n
        juvs <- matrix(NA,nrow=n,ncol=slotspersite) # track up to 100 juveniles at each site
        adults <- matrix(NA,nrow=n,ncol=slotspersite) # track up to 100 adults at each site
        # loop over sites and disperse or not
        dispersers <- matrix(NA,nrow=n,ncol=3) # col 1 is genotype, col 2 is age, col 3 is home site
        for ( j in 1:n ) {
            adults[j,1] <- s[j]
            if ( !is.na(adults[j,1]) ) {
                # there is an adult, so ask if it stays
                stay <- ifelse( runif(1) < adults[j,1] , 1 , 0 )
                # make offspring
                goff <- adults[j,1]
                if ( sex==TRUE & n_adults > 1 ) {
                    # sexual reproduction, so make juvenile intermediate between parents, plus mutation
                    # if only 1 adult in population, will throw an error
                    parent2 <- sample( s[-j][!is.na(s[-j])] , size=1 )
                    if ( discrete==FALSE ) {
                        goff <- ( adults[j,1] + parent2 )/2
                    } else {
                        goff <- sample( c(adults[j,1],parent2) , 1 )
                    }
                }
                # handle mutation
                if ( discrete==FALSE ) {
                    juvs[j,1] <- fclip( rnorm( 1, mean=goff , sd=mu ) )
                } else {
                    juvs[j,1] <- sample( c(goff,1-goff) , size=1 , prob=c(1-mu,mu) )
                }
                # dispersal
                if ( stay==1 ) {
                    # offspring disperses
                    offspring <- juvs[j,1] # copy genotype
                    juvs[j,1] <- NA # remove juvenile from birth site
                    # check if survives dispersal. if so, place in dispersal queue
                    if ( runif(1) < pj ) {
                        k <- findslot( dispersers[,1] )
                        dispersers[k,] <- c( offspring , 0 , j )
                    }
                } else {
                    # adult disperses
                    adult <- adults[j,1] # copy genotype
                    adults[j,1] <- NA # remove adult from own site
                    # check if survives dispersal. if so, place in dispersal queue
                    if ( runif(1) < pa ) {
                        k <- findslot( dispersers[,1] )
                        dispersers[k,] <- c( adult , 1 , j )
                    }
                }
            }
        }
        # loop over surviving dispersers and place at sites
        for ( j in 1:n ) {
            if ( !is.na(dispersers[j,1]) ) {
                genotype <- dispersers[j,1]
                age <- dispersers[j,2]
                home <- dispersers[j,3]
                if ( local.dispersal==FALSE ) {
                    destination <- sample( (1:n)[-home] , size=1 )
                } else {
                    site1 <- home-1
                    site2 <- home+1
                    site1 <- ifelse( site1 < 1 , n , site1 )
                    site2 <- ifelse( site2 > n , 1 , site2 )
                    destination <- sample( c(site1,site2) , size=1 )
                }
                if ( age==0 ) {
                    k <- findslot( juvs[destination,] )
                    juvs[destination,k] <- genotype
                } else {
                    k <- findslot( adults[destination,] )
                    adults[destination,k] <- genotype
                }
            }
        }
        # loop over sites and compete
        for ( j in 1:n ) {
            numadults <- countfilled( adults[j,] )
            numjuvs <- countfilled( juvs[j,] )
            genotypes <- c( adults[ j , !is.na(adults[j,]) ] , juvs[ j , !is.na(juvs[j,]) ] )
            rhps <- c( rep( va , length(adults[ j , !is.na(adults[j,]) ]) ) , rep( 1 , length(juvs[ j , !is.na(juvs[j,]) ]) ) )
            kages <- c( rep( 1 , length(adults[ j , !is.na(adults[j,]) ]) ) , rep( 0 , length(juvs[ j , !is.na(juvs[j,]) ]) ) )
            if ( length(genotypes) > 0 ) {
                # winner <- sample( genotypes , size=1 )
                winner <- sample( 1:length(genotypes) , size=1 , prob=rhps )
                winner.g <- genotypes[winner]
                winner.age <- kages[winner]
                # erase all adults and juvs at site, and place winner in adult slot
                adults[j,] <- rep(NA,slotspersite)
                juvs[j,] <- rep(NA,slotspersite)
                # store winner, if survives
                pr.die <- ma
                if ( winner.age == 0 ) pr.die <- mj
                if ( runif(1) < pr.die ) winner.g <- NA
                adults[j,1] <- winner.g
            } # else no one at site, so no winner
        }
        # loop over sites and process mortality
        #for ( j in 1:n ) {
        #    if ( runif(1) < ma ) {
        #        adults[j,1] <- NA
        #    }
        #}
        # update s record vector
        for ( j in 1:n ) {
            s[j] <- fclip( adults[j,1] )
        }
        if ( progress==TRUE ) progtext( i , 0 , tmax , st )
    } # tmax
    # analytical prediction
    shat <- ( pa * va * ( 2 + (-1+pa)*va ) - pj*(1+pa*va) )/( (-2 + pj - (-2 + pa)*va)*(pj-pa*va) )
    shat <- fclip( shat )
    # results
    if ( progress==TRUE ) system("say sim done")
    list( pa=pa , pj=pj , va=va , n=n , mu=mu , ma=ma , mj=mj , shist=shist , shistfull=shistfull , occupancyhist=occupancyhist , shat=shat )
}

library(compiler)
lbeqsim <- cmpfun(beqsim)

##########################################################
# example 

x <- lbeqsim( tmax=2000 , pa=0.77 , pj=0.55 , va=2^1 , ma=1-0.65 , mj=1-0.8 , n=1000 , sstart=0.7 , mu=0.001 , local.dispersal=FALSE , progress=TRUE , sex=TRUE , discrete=TRUE )

# plot full grid of occupancy, generations on x-axis, site on y-axis
plot_grid <- function( x ) {
    n <- ncol(x$shistfull)
    tmax <- length(x$shist)
    
    quartzFonts(serif = quartzFont(rep("Helvetica", 4)))
    par(family="serif")
    
    plot( NULL , col="white" , ylim=c(1,n) , xlim=c(1,tmax) , xlab="generation" , ylab="site" , yaxt='n' , xaxp=c(1,tmax,1) )
    for ( i in seq(from=1,to=tmax,by=2) ) { # loop over generations
        pt_col <- ifelse( is.na(x$shistfull[i,]) , "white" , ifelse( x$shistfull[i,]==0 , "orange" , rgb(0,0.55,1) ) )
        points( rep(i,n) , 1:n , cex=0.3 , pch=16 , col=pt_col )
    }
}

plot_grid( x )

# more typical time series plots

col.alpha <- function( acol , alpha=0.2 ) {
    acol <- col2rgb(acol)
    acol <- rgb(acol[1]/255,acol[2]/255,acol[3]/255,alpha)
    acol
}

n <- ncol(x$shistfull)
tmax <- length(x$shist)
quartzFonts(serif = quartzFont(rep("Helvetica", 4)))
par(family="serif")
plot( (x$shist) ,type="n", ylim=c(0,1) , xlim=c(0,tmax) , xlab="" , ylab="" , yaxt="n" , xaxp=c(0,tmax,1) )
axis( 2 , at=c(0,1) , labels=c(0,1) )

lines( 1:2000 , psim , col=col.alpha("red",1) , lty=1 , lwd=0.5 )
lines( 1:2000 , Rsim , col=col.alpha("black",1) , lty=1 , lwd=0.5 )

lines( 1:length(x$shist) , 1-x$shist , col=col.alpha("orange",0.6) , lty=1 , lwd=1 )
lines( 1:length(x$shist) , x$occupancyhist , col=col.alpha("black",0.6) , lty=1 , lwd=1 )
