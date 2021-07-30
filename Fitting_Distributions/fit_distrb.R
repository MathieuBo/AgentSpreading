setwd("C:\\Users\\thoma\\Documents\\M1_Neurasmus\\NeuroBIM_M1\\Internship\\GitRepo\\AgentModel") 
basedir <- getwd()
#install.packages("fitdistrplus")#Added 
#install.packages("logspline")#Added 
library(fitdistrplus)
library(logspline)
library(R.matlab)


# SYNUCLEIN
synexp <- read.csv(".\\SncaExpression.csv", header = FALSE)
synexp_val <- synexp$V2
synexp_val

  # Saving as PNG
png(".\\Fitting_Distribution\\Cullen_Frey_Synuclein.png")
  # Code
descdist(synexp_val, discrete = FALSE)
  # Close device
dev.off()

# Fitting
## Normal?
fit.norm <- fitdist(synexp_val, "norm")
  # Saving as PNG
png("norm_Synuclein.png")
  # Code
plot(fit.norm)
# Close device
dev.off()

## Uniform?
fit.unif <- fitdist(synexp_val, "unif")
# Saving as PNG
png("unif_Synuclein.png")
# Code
plot(fit.unif)
# Close device
dev.off()

## Log-normal
## Normal?
fit.lognorm <- fitdist(synexp_val, "lnorm")
# Saving as PNG
png("log_norm_Synuclein.png")
# Code
plot(fit.lognorm)
# Close device
dev.off()



# VOLUME
volume <- read.csv(".\\Volume_per_regions.csv", header = FALSE)
volume_val <- volume$V2
volume_val
log_volume <- log10(volume$V2)
maxi_log <- max(log_volume)
mini_log <- min(log_volume)
norm_log <- (log_volume-mini_log)/(maxi_log-mini_log)
norm_log

maxi_vol <- max(volume_val)
mini_vol <- min(volume_val)
norm_vol <- (volume_val-mini_vol)/(maxi_vol-mini_vol)

# Saving as PNG
png("Cullen_Frey_Volume.png")
descdist(volume_val, discrete = FALSE)
dev.off()

png("Cullen_Frey_log_Volume.png")
descdist(log_volume, discrete = FALSE)
dev.off()

# Fitting
#Log
## beta
fit.beta.normlog <- fitdist(norm_log, "beta", method="mme") 
png("beta_normlogvolume.png")
plot(fit.beta.normlog)
dev.off()

fit.beta.normvol <- fitdist(norm_vol, "beta", method="mme") 
png("beta_normvolume.png")
plot(fit.beta.normvol)
dev.off()

## Normal
fit.norm.logvol <- fitdist(log_volume, "norm") 
png("norm_logvolume.png")
plot(fit.norm.logvol)
dev.off()

# Non Log
##gamma
fit.gamma.vol <- fitdist(volume_val, "gamma", method="mme") 
png("gamma_volume.png")
plot(fit.gamma.vol)
dev.off()




# CONNECTOME
Connectome <- read.csv(".\\W.csv", header = FALSE)
#Connectome <- Connectome +0.000000001
Connectome
Connectome_val <- unlist(Connectome, recursive = TRUE, use.names = FALSE)
Connectome_val <- Connectome_val[Connectome_val > 0]
Connectome_val
png("Cullen_Frey_Connectome_wX0.png")
descdist(Connectome_val, discrete = FALSE)
dev.off()

minus_log_Connectome_val <- -log10(Connectome_val)
minus_log_Connectome_val
minus_log_Connectome_val <- minus_log_Connectome_val[minus_log_Connectome_val != "Inf"] 
minus_log_Connectome_val
png("Cullen_Frey_minuslogConnectome.png")
descdist(minus_log_Connectome_val, discrete = FALSE)
dev.off()

##log connectome
#-> norm
fit.norm.log.connectome <- fitdist(minus_log_Connectome_val, "norm")
png("norm_logConnectome_val.png")
plot(fit.norm.log.connectome)
dev.off()

#-> beta
max_abs = max(minus_log_Connectome_val)
min_abs = min(minus_log_Connectome_val)
norm_minus_log_Connectome_val = (minus_log_Connectome_val-min_abs)/(max_abs-min_abs)
fit.beta.log.connectome <- fitdist(norm_minus_log_Connectome_val, "beta", method="mme")
png("beta_normalogConnectome_val.png")
plot(fit.beta.log.connectome)
dev.off()

#->gamma
fit.gamma.log.connectome <- fitdist(abs(minus_log_Connectome_val), "gamma")
png("gamma_abslogConnectome_val.png")
plot(fit.gamma.log.connectome)
dev.off()


# Non log
## Gamma
fit.gamma<- fitdist(Connectome_val, "gamma", method="mme")
# Saving as PNG
png("gamma_Connectome_val.png")
# Code
plot(fit.gamma)
# Close device
dev.off()

maxi_Connectome_val <- max(Connectome_val)
mini_Connectome_val <- min(Connectome_val)
norm_Connectome_val <- (Connectome_val-mini_Connectome_val)/(maxi_Connectome_val-mini_Connectome_val )
norm_Connectome_val

## Beta
fit.beta<- fitdist(norm_Connectome_val, "beta", method="mme")
# Saving as PNG
png("beta_Connectome_val.png")
# Code
plot(fit.beta)
# Close device
dev.off()

