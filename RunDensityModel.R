# This piece of code runs the inla model described in Spatially structured eco-evolutionary dynamics in a host-pathogen interaction render isolated populations vulnerable to disease 
rm(list=ls())
library(sp)
library(INLA)
library(ggplot2)

INLA:::inla.dynload.workaround() 
rm(list = ls())


# load mesh for INLA
load(file = 'new_mesh.Rdata')
  

plot(mesh_c_islands)
mesh <- mesh_c_islands


# Load data for statistical modeling:
load(file = 'DataForDensityAnalyses.Rfile')   # contains files: d, sp1_gps_matrix

# plot connectivity values
png('Connectivity.png', height = 4, width = 5, units = 'in', res= 300)
hist(d$Sh, col= 'green', border = 'green', xlab = 'Connectivity', ylab = 'Frequency in monitoring data', main='')
lines(c(q1, q1), c(0, 3500))
lines(c(q2, q2), c(0, 3500))
dev.off()

# do data modifications:

sp1_gps_matrix <- sp1_gps_matrix[inds,]
d$rep <- as.numeric(d$PathPresence)+1 

d$Intercept <- 1 
d$AbsenceLowConnectivity <- 1*((d$PAPrev ==0))*(d$ShCat==1)
d$AbsenceHighConnectivity <- 1*((d$PAPrev ==0))*(d$ShCat==3)
d$AbsenceIntermediateConnectivity <- 1*((d$PAPrev ==0))*(d$ShCat==2)
d$PresenceHighConnectivity <- 1*((d$PAPrev ==1))*(d$ShCat==3)
d$PresenceIntermediateConnectivity <-1*((d$PAPrev ==1))*(d$ShCat==2)
d$PresenceLowConnectivity <-1*((d$PAPrev ==1))*(d$ShCat==1)
d$PathPresence11 <- 1*((d$PAPrev ==1))
d$PAPrev1 <- as.numeric(d$PAPrev)-1
d$rep <- as.numeric(d$PAPrev)+1


#################################################################
## make spde id's and stacks (needed for INLA):

sigma0 = 1
size = min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2]))))
range0 = size / 5
kappa0 = sqrt(8) / range0
tau0 = 1 / (sqrt(4 * pi) * kappa0 * sigma0)
spde = inla.spde2.matern(mesh, B.tau = cbind(log(tau0), -1, +1),   B.kappa = cbind(log(kappa0), 0, -1), theta.prior.mean = c(0, 0), theta.prior.prec = c(0.1, 1) )
n.years = length(unique(d$Year))

N_vars <- 7

d$Intercept = 1
group.years <-   as.integer(as.numeric(as.character(d$Year))) -1    #.....#    #HERE OUT A VECTOR DESCRIBING YEAR FOR EVERY OBSERVATION  

spde.idx <- inla.spde.make.index("spatial", n.spde=mesh$n,  n.group= n.years)
Atemp2 <- inla.spde.make.A(mesh=mesh, loc=sp1_gps_matrix, group = group.years)
effects2 <- list(spde.idx,  list(d))


h.spec <- list(theta=list(initial = 0.1, param =c(0,5)))
sigma <- 5
hyper=spde$f$hyper.default
hyper$theta1$initial = hyper$theta1$initial - log(sigma)
hyper$theta1$param[1] = hyper$theta1$param[1] - log(sigma)


####################################################################################################################################################


data <- inla.stack(data=list(Y=(d$RelChange)), A = list(Atemp2, 1), effects2, tag="obs")


save(file = 'Inputs_INLA_DensityModel.Rdata', mesh, spde.idx, sp1_gps_matrix, d, size, data, Atemp2)
inla.setOption(scale.model.default = TRUE)
hyper.prec = list(prec = list(param = c(1, 0.05)))


# run statistical model:

res <- inla(Y ~ -1+ f(spatial, hyper = hyper, model = spde, group = spatial.group, control.group =list(model='ar1', hyper=h.spec)) +
                 Rainfall_July+
                 Rainfall_August+
                 PresenceLowConnectivity+
                 PresenceHighConnectivity+
                 PresenceIntermediateConnectivity +
                 AbsenceLowConnectivity+
                 AbsenceHighConnectivity+
                 AbsenceIntermediateConnectivity +
                 LogPlant.drynessPrev+
                 LogPlant.dryness, family = 'gaussian', data = inla.stack.data(data), control.predictor =list(compute=TRUE,A=inla.stack.A(data)), control.fixed=list(expand.factor.strategy='inla'),control.inla = list(strategy = "gaussian", int.strategy = "eb"),  verbose = TRUE, control.compute=list(config = TRUE))



save(file = 'DensityModel.Rdata', res)


# explore results

summary(res)
plot(res)








