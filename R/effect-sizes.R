# # Effect sizes
# 
# 
# library(lavaan); library(semTools)
# 
# conf <- "
# f1 =~ 1*x1 + x2 + x3
# f2 =~ 1*x4 + x5 + x6
# 
# "
# 
# weak <- "
# f1 =~ NA*x1 + x2 + x3
# f2 =~ NA*x4 + x5 + x6
# f1 ~~ c(1, NA)*f1
# f2 ~~ c(1, NA)*f2
# "
# 
# configural <- cfa(conf, data = HolzingerSwineford1939, group="school")
# weak <- cfa(weak, data = HolzingerSwineford1939, group="school", group.equal="loadings")
# 
# dmacs <- function(lav.model) {
#   dmacs::lavaan_dmacs(configural, RefGroup=lavInspect(configural, "group.label")[1], "pooled" )$DMACS
# }
# 
# 
# dmacs::lavaan_dmacs(weak, RefGroup=lavInspect(weak, "group.label")[1], "pooled" )$DMACS
# 
# # Dmacs
# library(lavaan)
# Dmacs <- function(pairwise.lavaan.model.fit) {
#   
#   weak = pairwise.lavaan.model.fit
#   summary(weak)
#   fit = weak
#   
#   mean.dmacs = mean(as.matrix(dmacs::lavaan_dmacs(weak, RefGroup=lavInspect(weak, "group.label")[1], "pooled" )$DMACS), na.rm=T)
#   
# }
# 
# 
# dmacs_lavaan(weak, signed = F)
# 
# pf.signed.true = pairwiseFit("f1 =~ v1 + v2 + v3 + v4;
#                f2 =~ v11 + v12 + v13 + v14;",  data = four.clusters, group="group", signed = T)
# pf.unsigned.true = pairwiseFit("f1 =~ v1 + v2 + v3 + v4;
#                f2 =~ v11 + v12 + v13 + v14;",  data = four.clusters, group="group", signed = F)
# 
# 
# plotDistances(measures = pf.signed.true, fit.index="signed.dmacs" )
# plotDistances(measures = pf.unsigned.true, fit.index="average.dmacs" )
# 
# # 
# # 
# # SDI2 <- function(LambdaR, ThreshR, LambdaF, ThreshF, MeanF, VarF, SD) {
# #   
# #   expected_value <- function(Lambda, Thresh, Theta) { Thresh + Lambda * Theta }
# #   z <- seq(-5, 5, .001)
# #   
# #   Y_j1 = expected_value(LambdaF, ThreshF, MeanF + z * sqrt(VarF))
# #   Y_j2 = expected_value(LambdaR, ThreshR, MeanF + z * sqrt(VarF))
# #   
# #   sum((Y_j1- Y_j2) * dnorm(z) * .001 * sqrt(VarF))/SD
# # }
# # 
# # 
# # 
# # 
# 
# 
# 
# 
# 
# 
