# 
# library(LittleHelpers)
# library(MIE)
# 
# ess6 <- haven::read_sav("~/Dropbox/STAT/European Social Survey/Data/ESS6/ESS6e02_2.sav")
# 
# ess_trimmed <- ess6[, c( values$items, "cntry")]
# class(ess_trimmed) <- "data.frame"
# for(x in names(ess_trimmed)) attributes(ess_trimmed[,x])<-NULL
# ess_trimmed$cntry <- as.factor(ess_trimmed$cntry)
# 
# ess_trimmed$country <- ess_trimmed$cntry
# ess_trimmed$cntry <- NULL
# 
# 
# runMIE("F =~ impfree + iphlppl + ipsuces + ipstrgv", ess_trimmed[1:5000,], group="country")
# runMIE(model = )
# runMIE(data=ess_trimmed[1:5000,], group = "country")
# runMIE(data=ess6[, c()], group = "country")
# 
# runMIE()
# 
# 
# covars <- compute_covariance(ess_trimmed, group="country")
# crs <-compute_correlation(ess_trimmed, "country")
# 
# intrcps <- MGCFAparameters(model  = "F =~ impfree + iphlppl + ipsuces + ipstrgv", data=ess_trimmed,
#                            group="country", parameters = "intercepts")
# 
# 
# 
# mt1 <- pairwiseFit(model="F =~ impfree + iphlppl + ipsuces + ipstrgv",
#                    data = ess_trimmed,
#                    group="country")
# 
# mt2 <- pairwiseFit(data = ess_trimmed,  model="F =~ impfree + iphlppl + ipsuces + ipstrgv", shiny=F, group="country", constraints = "intercepts" )
# 
# 
# a <- incrementalFit("F =~ impfree + iphlppl + ipsuces + ipstrgv", ess_trimmed, group="country")
# 
# 
# plotDistances(covars)
# plotDistances(crs)
# plotDistances(intrcps)
# plotDistances(a, n.clusters = 3, drop = c("IE", "LT", "CY", "BG"), fit.index = "fit")
# 
# ddd <- read.csv("inst/application/simulated2.csv")
# m1 <- pairwiseFit(ddd, pairs_of_countries(ddd$cntry[1:5000]),
#     model = "#By default a simulated data - model mimics Schwartz values ESS scale
#   person=~ ipcrtiv +impfree +impfun +ipgdtim +impdiff +ipadvnt+ imprich +iprspot +ipshabt +ipsuces;
#     social=~ impenv +ipeqopt +ipudrst +iplylfr +iphlppl +impsafe +ipstrgv +ipfrule +ipbhprp +ipmodst +imptrad;",
#    shiny=F, group="cntry", constraints = c("loadings", "intercepts"))
# str(m1)



