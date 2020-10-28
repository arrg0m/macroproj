## LAG NUMBER p = 4 IS FIXED

library(readxl)
library(AER)
library(zoo)
library(mFilter)
## DATA PREPROCESSING: IMPORT
inv_gnp_work <- read_excel("C:/Users/hmbyun/Dropbox/Work/inv-gnp-work.xlsx", 
                           sheet = "SVAR Export")
dp <- inv_gnp_work$dp
da <- inv_gnp_work$da
h  <- inv_gnp_work$h
r  <- inv_gnp_work$r
pi <- inv_gnp_work$pi
# date <- inv_gnp_work$DATE

dp <- zoo(dp)
da <- zoo(da)
h <- zoo(h)
r <- zoo(r)
pi <- zoo(pi)

#start.date <- c(1955, 1); end.date <- c(1979, 2)
start.date <- c(1982, 3); end.date <- c(2000, 4)
strr <- paste("US, ", start.date[1], "Q", start.date[2]," - ", end.date[1], "Q", end.date[2], sep = "")

## DETRENDING AND NORMALIZING
myPrep <- function(x, dt = FALSE, nr = TRUE){
  if(dt == TRUE){
    x <- hpfilter(x, freq = 1600, type = "lambda")$cycle
  }
  if(nr == TRUE){
    x <- (x - mean(x))/sqrt(var(x))
  }
  return(x)
}


## DATA PREPROCESSING: CUT AND NORMALIZE

dp.ts <- ts(data = dp, start = c(1947, 2), deltat = 1/4); dp.ts <- window(dp.ts, start = start.date, end = end.date)
da.ts <- ts(data = da, start = c(1947, 2), deltat = 1/4); da.ts <- window(da.ts, start = start.date, end = end.date)
h.ts  <- ts(data = h , start = c(1947, 2), deltat = 1/4); h.ts <- window(h.ts, start = start.date, end = end.date)
r.ts  <- ts(data = r , start = c(1947, 2), deltat = 1/4); r.ts <- window(r.ts, start = start.date, end = end.date)
pi.ts <- ts(data = pi, start = c(1947, 2), deltat = 1/4); pi.ts <- window(pi.ts, start = start.date, end = end.date)
dp <- myPrep(dp.ts, dt = FALSE, nr = TRUE); dp <- zoo(dp)
da <- myPrep(da.ts, dt = FALSE, nr = TRUE); da <- zoo(da)
h  <- myPrep(h.ts, dt = FALSE , nr = TRUE); h  <- zoo( h)
r  <- myPrep(r.ts, dt = FALSE, nr = TRUE); r  <- zoo( r)
pi <- myPrep(pi.ts, dt = FALSE , nr = TRUE); pi <- zoo(pi)
# dp <- zoo(dp.ts); da <- zoo(da.ts); h <- zoo(h.ts); r <- zoo(r.ts); pi <- zoo(pi.ts)

## DATA PREPROCESSING: CONSTRUCT DERIVED VARIABLES

dp1 <- lag(dp, -1, na.pad = TRUE)
dp2 <- lag(dp, -2, na.pad = TRUE)
dp3 <- lag(dp, -3, na.pad = TRUE)
dp4 <- lag(dp, -4, na.pad = TRUE)

da1 <- lag(da, -1, na.pad = TRUE)
da2 <- lag(da, -2, na.pad = TRUE)
da3 <- lag(da, -3, na.pad = TRUE)
da4 <- lag(da, -4, na.pad = TRUE)
dda <- diff(zoo(da), na.pad = TRUE)
dda1 <- lag(dda, -1, na.pad = TRUE)
dda2 <- lag(dda, -2, na.pad = TRUE)
dda3 <- lag(dda, -3, na.pad = TRUE)

h1 <- lag(h, -1, na.pad = TRUE)
h2 <- lag(h, -2, na.pad = TRUE)
h3 <- lag(h, -3, na.pad = TRUE)
h4 <- lag(h, -4, na.pad = TRUE)
dh <- diff(zoo(h), na.pad = TRUE)
dh1 <- lag(dh, -1, na.pad = TRUE)
dh2 <- lag(dh, -2, na.pad = TRUE)
dh3 <- lag(dh, -3, na.pad = TRUE)

r1 <- lag(r, -1, na.pad = TRUE)
r2 <- lag(r, -2, na.pad = TRUE)
r3 <- lag(r, -3, na.pad = TRUE)
r4 <- lag(r, -4, na.pad = TRUE)
dr <- diff(zoo(r), na.pad = TRUE)
dr1 <- lag(dr, -1, na.pad = TRUE)
dr2 <- lag(dr, -2, na.pad = TRUE)
dr3 <- lag(dr, -3, na.pad = TRUE)

pi1 <- lag(pi, -1, na.pad = TRUE)
pi2 <- lag(pi, -2, na.pad = TRUE)
pi3 <- lag(pi, -3, na.pad = TRUE)
pi4 <- lag(pi, -4, na.pad = TRUE)
dpi <- diff(zoo(pi), na.pad = TRUE)
dpi1 <- lag(dpi, -1, na.pad = TRUE)
dpi2 <- lag(dpi, -2, na.pad = TRUE)
dpi3 <- lag(dpi, -3, na.pad = TRUE)

## RUN FIVE IV REGRESSIONS AS IN FISHER(2006)

eqn1iv <- ivreg(dp ~ -1 + dp1 + dp2 + dp3 + dp4 + dda + dda1 + dda2 + dda3 + dh + dh1 + dh2 + dh3 + dr + dr1 + dr2 + dr3 + dpi + dpi1 + dpi2 + dpi3 
                        | dp1 + dp2 + dp3 + dp4 + da1 + da2 + da3 + da4 + h1 + h2 + h3 + h4 + r1 + r2 + r3 + r4 + pi1 + pi2 + pi3 + pi4)
summary(eqn1iv)
eqn1iv.coef <- eqn1iv$coefficients
eqn1iv.zres <- zoo(eqn1iv$residuals)

eqn2iv <- ivreg(da ~ 0 + dp + dp1 + dp2 + dp3 + dp4 + da1 + da2 + da3 + da4 + dh + dh1 + dh2 + dh3 + dr + dr1 + dr2 + dr3 + dpi + dpi1 + dpi2 + dpi3 
                | dp1 + dp2 + dp3 + dp4 + da1 + da2 + da3 + da4 + h1 + h2 + h3 + h4 + r1 + r2 + r3 + r4 + pi1 + pi2 + pi3 + pi4 + eqn1iv.zres)
summary(eqn2iv)
eqn2iv.coef <- eqn2iv$coefficients
eqn2iv.zres <- zoo(eqn2iv$residuals)

eqn3iv <- ivreg(h ~ 0 + dp + dp1 + dp2 + dp3 + dp4
                + da + da1 + da2 + da3 + da4
                + h1 + h2 + h3 + h4
                + r + r1 + r2 + r3 + r4
                + pi + pi1 + pi2 + pi3 + pi4 | dp1 + dp2 + dp3 + dp4 + da1 + da2 + da3 + da4 + h1 + h2 + h3 + h4 + r1 + r2 + r3 + r4 + pi1 + pi2 + pi3 + pi4 + eqn1iv.zres + eqn2iv.zres)
summary(eqn3iv)
# GMM instead
eqn3iv.coef <- eqn3iv$coefficients
eqn3iv.zres <- zoo(eqn3iv$residuals)

eqn4iv <- ivreg(r ~ 0 + dp + dp1 + dp2 + dp3 + dp4
                + da + da1 + da2 + da3 + da4
                + h + h1 + h2 + h3 + h4
                + r1 + r2 + r3 + r4
                + pi + pi1 + pi2 + pi3 + pi4 | dp1 + dp2 + dp3 + dp4 + da1 + da2 + da3 + da4 + h1 + h2 + h3 + h4 + r1 + r2 + r3 + r4 + pi1 + pi2 + pi3 + pi4 + eqn1iv.zres + eqn2iv.zres + eqn3iv.zres)
summary(eqn4iv)
eqn4iv.coef <- eqn4iv$coefficients
eqn4iv.zres <- zoo(eqn4iv$residuals)

eqn5iv <- ivreg(pi ~ 0 + dp + dp1 + dp2 + dp3 + dp4
                + da + da1 + da2 + da3 + da4
                + h + h1 + h2 + h3 + h4
                + r + r1 + r2 + r3 + r4
                + pi1 + pi2 + pi3 + pi4 | dp1 + dp2 + dp3 + dp4 + da1 + da2 + da3 + da4 + h1 + h2 + h3 + h4 + r1 + r2 + r3 + r4 + pi1 + pi2 + pi3 + pi4 + eqn1iv.zres + eqn2iv.zres + eqn3iv.zres + eqn4iv.zres)
summary(eqn5iv)
eqn5iv.coef <- eqn5iv$coefficients
eqn5iv.zres <- zoo(eqn5iv$residuals)

## "REWIND" TO FIND VAR COEFFICIENTS

# Ay = Gamma(L)y + e
# Dimension Indexed by [to, from, p]
# e.g. Gamma[1,2,1]: coef of dp from da, with lag p = 1
A <- array(0, c(5,5))
Gamma <- array(0, c(5, 5, 4)) 

# INITIAL ASSUMPTIONS
A[1,1] = 1; A[2,2] = 1; A[3,3] = 1; A[4,4] = 1; A[5,5] = 1;

# FROM eqn1iv.coef
Gamma[1,1,] = eqn1iv.coef[c('dp1','dp2','dp3','dp4')]
A[1,c(2,3,4,5)] = -eqn1iv.coef[c('dda','dh','dr','dpi')]
Gamma[1,2,1] = -eqn1iv.coef['dda'] + eqn1iv.coef['dda1']
Gamma[1,2,2] = -eqn1iv.coef['dda1'] + eqn1iv.coef['dda2']
Gamma[1,2,3] = -eqn1iv.coef['dda2'] + eqn1iv.coef['dda3']
Gamma[1,2,4] = -eqn1iv.coef['dda3']
Gamma[1,3,1] = -eqn1iv.coef['dh'] + eqn1iv.coef['dh1']
Gamma[1,3,2] = -eqn1iv.coef['dh1'] + eqn1iv.coef['dh2']
Gamma[1,3,3] = -eqn1iv.coef['dh2'] + eqn1iv.coef['dh3']
Gamma[1,3,4] = -eqn1iv.coef['dh3']
Gamma[1,4,1] = -eqn1iv.coef['dr'] + eqn1iv.coef['dr1']
Gamma[1,4,2] = -eqn1iv.coef['dr1'] + eqn1iv.coef['dr2']
Gamma[1,4,3] = -eqn1iv.coef['dr2'] + eqn1iv.coef['dr3']
Gamma[1,4,4] = -eqn1iv.coef['dr3']
Gamma[1,5,1] = -eqn1iv.coef['dpi'] + eqn1iv.coef['dpi1']
Gamma[1,5,2] = -eqn1iv.coef['dpi1'] + eqn1iv.coef['dpi2']
Gamma[1,5,3] = -eqn1iv.coef['dpi2'] + eqn1iv.coef['dpi3']
Gamma[1,5,4] = -eqn1iv.coef['dpi3']

# FROM eqn2iv.coef
A[2,c(1,3,4,5)] = -eqn2iv.coef[c('dp','dh','dr','dpi')]
Gamma[2,1,] = eqn2iv.coef[c('dp1', 'dp2', 'dp3', 'dp4')]
Gamma[2,2,] = eqn2iv.coef[c('da1', 'da2', 'da3', 'da4')]
Gamma[2,3,1] = -eqn2iv.coef['dh'] + eqn2iv.coef['dh1']
Gamma[2,3,2] = -eqn2iv.coef['dh1'] + eqn2iv.coef['dh2']
Gamma[2,3,3] = -eqn2iv.coef['dh2'] + eqn2iv.coef['dh3']
Gamma[2,3,4] = -eqn2iv.coef['dh3']
Gamma[2,4,1] = -eqn2iv.coef['dr'] + eqn2iv.coef['dr1']
Gamma[2,4,2] = -eqn2iv.coef['dr1'] + eqn2iv.coef['dr2']
Gamma[2,4,3] = -eqn2iv.coef['dr2'] + eqn2iv.coef['dr3']
Gamma[2,4,4] = -eqn2iv.coef['dr3']
Gamma[2,5,1] = -eqn2iv.coef['dpi'] + eqn2iv.coef['dpi1']
Gamma[2,5,2] = -eqn2iv.coef['dpi1'] + eqn2iv.coef['dpi2']
Gamma[2,5,3] = -eqn2iv.coef['dpi2'] + eqn2iv.coef['dpi3']
Gamma[2,5,4] = -eqn2iv.coef['dpi3']

# FROM eqn3iv.coef
A[3,c(1,2,4,5)] = -eqn3iv.coef[c('dp','da','r','pi')]
Gamma[3,1,] = eqn3iv.coef[c('dp1','dp2','dp3','dp4')]
Gamma[3,2,] = eqn3iv.coef[c('da1','da2','da3','da4')]
Gamma[3,3,] = eqn3iv.coef[c('h1','h2','h3','h4')]
Gamma[3,4,] = eqn3iv.coef[c('r1','r2','r3','r4')]
Gamma[3,5,] = eqn3iv.coef[c('pi1','pi2','pi3','pi4')]

# FROM eqn4iv.coef
A[4,c(1,2,3,5)] = -eqn4iv.coef[c('dp','da','h','pi')]
Gamma[4,1,] = eqn4iv.coef[c('dp1','dp2','dp3','dp4')]
Gamma[4,2,] = eqn4iv.coef[c('da1','da2','da3','da4')]
Gamma[4,3,] = eqn4iv.coef[c('h1','h2','h3','h4')]
Gamma[4,4,] = eqn4iv.coef[c('r1','r2','r3','r4')]
Gamma[4,5,] = eqn4iv.coef[c('pi1','pi2','pi3','pi4')]

# FROM eqn5iv.coef
A[5,c(1,2,3,4)] = -eqn5iv.coef[c('dp','da','h','r')]
Gamma[5,1,] = eqn5iv.coef[c('dp1','dp2','dp3','dp4')]
Gamma[5,2,] = eqn5iv.coef[c('da1','da2','da3','da4')]
Gamma[5,3,] = eqn5iv.coef[c('h1','h2','h3','h4')]
Gamma[5,4,] = eqn5iv.coef[c('r1','r2','r3','r4')]
Gamma[5,5,] = eqn5iv.coef[c('pi1','pi2','pi3','pi4')]

Gamma[3,5,4] <- 0 # this value was of NA; Temporally set to zero till find problem (WO intercept)
# Gamma[3,5,3] <- 0; Gamma[3,5,4] <- 0; Gamma[4,5,4] <- 0; # W/ Intercept
solve(A - Gamma[,,1] - Gamma[,,2] - Gamma[,,3] - Gamma[,,4]) # Double check LR Identification schemes

B <- diag(sqrt(c(var(eqn1iv.zres), var(eqn2iv.zres), var(eqn3iv.zres), var(eqn4iv.zres), var(eqn5iv.zres))))

# EVALUATE IRFs
myIRF <- function(G0, G1, G2, G3, G4, nLags = 48){
  C <- array(0, c(5, 5, ifelse(nLags >= 4, nLags + 1, 5)))
  C[,,1] = solve(G0)
  C[,,2] = -solve(G0, G1%*%C[,,1])
  C[,,3] = -solve(G0, G2%*%C[,,2] + G1%*%C[,,1])
  C[,,4] = -solve(G0, G3%*%C[,,3] + G2%*%C[,,2] + G1%*%C[,,1])
  for(k in 4:nLags){
    C[,,k+1] = -solve(G0, G4%*%C[,,k-3] + G3%*%C[,,k-2] + G2%*%C[,,k-1] + G1%*%C[,,k])
  }
  if(nLags >= 4)
    return(C[,,])
  else
    return(C[,,1:nLags])
}
#C <- myIRF(solve(B,A), -solve(B,Gamma[,,1]), -solve(B,Gamma[,,2]), -solve(B,Gamma[,,3]), -solve(B,Gamma[,,4]))
C <- myIRF(A, -Gamma[,,1], -Gamma[,,2], -Gamma[,,3], -Gamma[,,4])
plot(C[2,2,])
