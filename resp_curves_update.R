## response_curves
# # data simulation
# set.seed(10)
# n_root <- 50
# n <- n_root * n_root
# clim1 <- rnorm(n, 10)
# clim2 <- rnorm(n)
# coors <- expand.grid(1:n_root,1:n_root)
# occ1 <- rbinom(n, 1, invLogit(-17 + 1.5 * clim1 - 2 * clim2))
# dat1 <- data.frame(coors, clim1, clim2, occ1)
# coordinates(dat1) = c("Var1", "Var2")
# gridded(dat1) <- TRUE
# dat1 <- as(dat1, "SpatialGridDataFrame")
# 
# clim1 <- rnorm(n, 8)
# clim2 <- rnorm(n, 1)
# occ2 <- rbinom(n, 1, invLogit(-13 + 1 * clim1 - 1.5 * clim2))
# dat2 <- data.frame(coors, clim1, clim2, occ2)
# coordinates(dat2) = c("Var1", "Var2")
# gridded(dat2) <- TRUE
# dat2 <- as(dat2, "SpatialGridDataFrame")
# 
# # model running
# mod1 <- sdmINLASpatRandEff(dat1[3], dat1[1:2], meshParameters = list(max.n.strict = -1), createGeoTIFF = FALSE, outFolder = getwd(), myName = "Species_occ1_ModelPredictions")
# mod2 <- sdmINLASpatRandEff(dat2[3], dat2[1:2], meshParameters = list(max.n.strict = -1), createGeoTIFF = FALSE, outFolder = getwd(), myName = "Species_occ2_ModelPredictions")
# 
# modRDS1 <- readRDS(file = "Species_occ1_ModelPredictions.rds")
# modRDS2 <- readRDS(file = "Species_occ2_ModelPredictions.rds")
# 
# preds_1 <-modRDS1$responsePredictions
# preds_2 <- modRDS2$responsePredictions
# 
# plot_response_curve(preds_1$clim1, "climate variable 1")
# 
# ggplot(preds_1$clim1, aes(covarVal, meanEst))+
#   geom_line(linewidth = 0.6, col = "blue")+
#   geom_line(aes(covarVal, lowerEst), linetype = "dashed")+
#   geom_line(aes(covarVal, upperEst), linetype = "dashed")+
#   # ylim(y_limit)+
#   labs(y = "Probability of occurence", x = "climate variable 1", title = "Native model")+
#   theme_bw()
# plot_response_curve(preds_2$clim2, "climate variable 1")
# ggplot(preds_2$clim1, aes(covarVal, meanEst))+
#   geom_line(linewidth = 0.6, col = "red")+
#   geom_line(aes(covarVal, lowerEst), linetype = "dashed")+
#   geom_line(aes(covarVal, upperEst), linetype = "dashed")+
#   ylim(c(0, 1))+
#   labs(y = "Probability of occurence", x = "climate variable 1", title = "Invasive model")+
#   theme_bw()

# functions
invLogit <- function(x) exp(x)/(1+ exp(x))
st <- function(x, y = x) (x - mean(y, na.rm = TRUE)) / sd(y, na.rm = TRUE)
getPars <- function(modRDS, climateData){
  pars <- c(modRDS$modelObject$summary.fixed$mean,
            tail(modRDS$modelObject$summary.random$climIndeces$mean, length(climateData@data) * 2))
  climNames <- names(climateData@data)
  names(pars) <- c("Intercept", climNames, paste0(climNames, "_quadratic"))
  pars
}
respCurve <- function(x, datSel, pars, climVar, integral = 1){ 
  x <- st(x, datSel)
  invLogit(pars["Intercept"] + x * pars[climVar] + x * x * pars[paste0(climVar, "_quadratic")]) / integral
}
jointCurves <- function(x, datSel1, datSel2, pars1, pars2, integral1, integral2, climVar, fn){
  fn(respCurve(x, datSel1, pars1, climVar, integral1), respCurve(x, datSel2, pars2, climVar, integral2))  
}
calcOverlap <- function(climVar, modRDS1, climDat1, modRDS2, climDat2){
  getAux <- function(modRDS, climDat, climVar){
    xSeq <- modRDS$responsePredictions[[climVar]]$covarVal
    pars <- getPars(modRDS, climDat)
    datSel <- climDat[[climVar]]
    integral <- integrate(respCurve, min(xSeq), max(xSeq), datSel, pars, climVar)$value
    list(xSeq = xSeq, pars = pars, datSel = datSel, integral = integral)
  }
  A <- getAux(modRDS1, climDat1, climVar)
  B <- getAux(modRDS2, climDat2, climVar)
  cRange <- c(max(min(A$xSeq), min(B$xSeq)), min(max(A$xSeq), max(B$xSeq)))
  argsL <- list(jointCurves, cRange[1], cRange[2], A$datSel, B$datSel, A$pars, B$pars, A$integral, B$integral, climVar)
  intMin <- do.call(integrate, c(argsL, pmin))$value
  intMax <- do.call(integrate, c(argsL, pmax))$value
  list(native = A, invasive = B, cRange = cRange, overlap = intMin / intMax)
}



# # calculation of overlap and figures
# overClim1 <- calcOverlap("clim1", modRDS1, dat1[1:2], modRDS2, dat2[1:2])
# overClim2 <- calcOverlap("clim2", modRDS1, dat1[1:2], modRDS2, dat2[1:2])
# 
# with(overClim1$native, plot(xSeq, respCurve(xSeq, datSel, pars, "clim1", integral), type = "l", col = "blue",
#   xlab = "climate variable 1", ylab = "Occurrence density", main = "Overlap plot"))
# with(overClim1$invasive, lines(xSeq, respCurve(xSeq, datSel, pars, "clim1", integral), col = "red"))
# abline(v = overClim1$cRange, lwd = 2, lty = 2)
# 
# with(overClim2$native, plot(xSeq, respCurve(xSeq, datSel, pars, "clim2", integral), type = "l",
#   xlab = "clim2", ylab = "Occurrence density"))
# with(overClim2$invasive, lines(xSeq, respCurve(xSeq, datSel, pars, "clim2", integral), col = "blue"))
# abline(v = overClim2$cRange, lwd = 2, lty = 2)
# 
# # all variables
# climVariables <- c("clim1", "clim2")
# overlap <- lapply(climVariables, calcOverlap, modRDS1, dat1[1:2], modRDS2, dat2[1:2])
# 
# _
