# install.packages(c("readxl", "vars", "ggplot2", "grid", "gridExtra"))

library(readxl)
library(vars)
library(ggplot2)
inv_gnp_work <- read_excel("C:/Users/hmbyun/Dropbox/Work/R/inv-gnp-work.xlsx", 
                           sheet = "SVAR Export")

dp <- inv_gnp_work$dp
da <- inv_gnp_work$da
h  <- inv_gnp_work$logh
r  <- inv_gnp_work$r
pi <- inv_gnp_work$pi

# date <- inv_gnp_work$DATE
start.date <- c(1955, 1); end.date <- c(1979, 2)
#start.date <- c(1982, 3); end.date <- c(2000, 4)
strr <- paste("US, ", start.date[1], "Q", start.date[2]," - ", end.date[1], "Q", end.date[2], sep = "")

dp.ts <- ts(data = dp, start = c(1947, 2), deltat = 1/4); dp.ts <- window(dp.ts, start = start.date, end = end.date)
da.ts <- ts(data = da, start = c(1947, 2), deltat = 1/4); da.ts <- window(da.ts, start = start.date, end = end.date)
h.ts  <- ts(data = h , start = c(1947, 2), deltat = 1/4); h.ts <- window(h.ts, start = start.date, end = end.date)
h.mean = mean(h.ts); h.ts <- h.ts - h.mean
# r.ts  <- ts(data = r , start = c(1947, 2), deltat = 1/4); r.ts <- window(r.ts, start = start.date, end = end.date)
# pi.ts <- ts(data = pi, start = c(1947, 2), deltat = 1/4); pi.ts <- window(pi.ts, start = start.date, end = end.date)

basic.ts <- ts(
  data = matrix(data = c(dp.ts, da.ts, h.ts), ncol = 3, byrow = FALSE),
  start = start.date, deltat = 1/4, names = c('dp', 'da', 'h')
)
basic.var <- VAR(y = basic.ts, p=4, type = 'none')
# basic.svar <- BQ(basic.var)
basic.svar <- SVAR(
  x <- basic.var,
  estmethod <- "direct",
  Amat <- matrix(c(NA, 0, 0, NA, NA, 0, NA, NA, NA), nrow=3, byrow=TRUE),
  Bmat <- matrix(c(NA, 0, 0, 0, NA, 0, 0, 0, NA), nrow=3, byrow=TRUE),
  lrtest = FALSE
)
## 2-variable case
basic.ts.2 <- ts(
  data = matrix(data = c(da.ts, h.ts), ncol = 2, byrow = FALSE),
  start = start.date, deltat = 1/4, names = c('da', 'h')
)
basic.var.2 <- VAR(y = basic.ts.2, p=4, type = 'none')
basic.svar.2 <- SVAR(
  x <- basic.var.2,
  estmethod <- "direct",
  Amat <- matrix(c(NA, 0, NA, NA), nrow=2, byrow=TRUE),
  Bmat <- matrix(c(NA, 0, 0, NA), nrow=2, byrow=TRUE),
  lrtest = FALSE
)
# BQ(basic.var.2)


#IRFS
Tmax = 48
#basic.irf <- irf(basic.var, n.ahead = Tmax, ortho = TRUE)
#basic.irf.2 <- irf(basic.var.2, n.ahead = Tmax, ortho = TRUE)
basic.irf <- irf(basic.svar, n.ahead = Tmax)
basic.irf.2 <- irf(basic.svar.2, n.ahead = Tmax)

dp.dp <- basic.irf$irf$dp[,1]; p.dp <- cumsum(dp.dp)
dp.da <- basic.irf$irf$dp[,2]; p.da <- cumsum(dp.da)
da.dp <- basic.irf$irf$da[,1]; a.dp <- cumsum(da.dp)
da.da <- basic.irf$irf$da[,2]; a.da <- cumsum(da.da)
h.dp <- basic.irf$irf$h[,1]
h.da <- basic.irf$irf$h[,2]

da.da.2 <- basic.irf.2$irf$da[,1]; a.da.2 <- cumsum(da.da.2)
h.da.2 <- basic.irf.2$irf$h[,1];


library(grid)
library(gridExtra)
# in place of Fig 4

timeIdx = 1:(Tmax+1)
x.grid.1 <- rep(c(0, timeIdx), 1)
x.grid.2 <- rep(c(0, timeIdx), 2)
x.grid.3 <- rep(c(0, timeIdx), 3)

df <- data.frame(x = x.grid.3, val = c(c(0, -h.dp), c(0, h.da), c(0, h.da.2)), variable = rep(c("Hour / Vt","Hour / At","Hour / At only"),each = Tmax + 2))
p1 <- ggplot(data = df, aes(x = x, y = val)) + geom_point(aes(colour=variable)) + theme(legend.position="bottom", legend.direction="vertical")

df <- data.frame(x = x.grid.3, val = c(c(0, - h.dp - a.dp), c(0, h.da + a.da), c(0,  h.da.2 + a.da.2)), variable = rep(c("Output / Vt","Output / At", "Output / At only"),each = Tmax + 2))
p2 <- ggplot(data = df, aes(x = x, y = val)) + geom_point(aes(colour=variable)) + theme(legend.position="bottom", legend.direction="vertical")

df <- data.frame(x = x.grid.3, val = c(c(0, p.dp), c(0,-p.da), 0*x.grid.1), variable = rep(c("Invt Price / Vt","Invt Price / At", "Invt Price / At only"),each = Tmax + 2))
p3 <- ggplot(data = df, aes(x = x, y = val)) + geom_point(aes(colour=variable)) + theme(legend.position="bottom", legend.direction="vertical")

df <- data.frame(x = x.grid.3, val = c(c(0, -a.dp), c(0, a.da), c(0, a.da.2)), variable = rep(c("Labor Prod / Vt","Labor Prod / At", "Labor Prod / At only"),each = Tmax + 2))
p4 <- ggplot(data = df, aes(x = x, y = val)) + geom_point(aes(colour=variable)) + theme(legend.position="bottom", legend.direction="vertical")

grid.arrange(p1, p2, p3, p4, ncol =4, top =strr)

# horizontal line at zero (center)
fevd(basic.svar)
