library(tseries)

#####################Chargement des données###################
milk_production <- read.table("monthly-milk-production.csv", header = TRUE);

#####################Dévision des données#####################
test_df<-milk_production[157:168,]
train_df<-milk_production[1:156,]

#####################Création de la série#####################
M<-train_df[,2]
M<- floor( as.numeric( gsub(",", "", M) ) )
train_ts <- ts(M, start=1962,deltat=1/12)

###################representation graphique###################
plot(train_ts, type="l", col="red",main="Production mensuelle du lait",xlab="Temps",ylab="Quantité de lait" )

############################ACF###############################
plot(acf(train_ts,lag.max=36,plot=FALSE),ylim=c(-1,1),main="ACF sur l'échantillon d'apprentissage")

##################test de stationnarité######################
adf.test(train_ts)
kpss.test(train_ts)


############différenciation (I-B) de la série##################
y_dif1=diff(train_ts,lag=1,differences=1)
plot(acf(y_dif1,lag.max=36,plot=FALSE),ylim=c(-1,1),main="ACF après différenciation")

###############test de stationnarité en tendance###############
adf.test(y_dif1)
kpss.test(y_dif1)


######différenciation (I-B^12) saisonnière+ ACF et PACF#########
y_dif_1_12=diff(y_dif1,lag=12,differences=1)
plot(pacf(y_dif_1_12,lag.max=48,plot=FALSE),ylim=c(-1,1),main="PACF après différenciation")
plot(acf(y_dif_1_12,lag.max=36,plot=FALSE),ylim=c(-1,1),main="ACF après différenciation")

####test de stationnarité de la série différenciée en saisonnalité############
adf.test(y_dif_1_12)
kpss.test(y_dif_1_12)
pp.test(y_dif_1_12)
PP.test(y_dif_1_12)


###################METHODE BOX & JENKINS########################

#1- IDENTIFICATION DU MODELE
plot(acf(y_dif_1_12,lag.max=36,plot=FALSE),ylim=c(-1,1),main="autocorrélogramme")
plot(pacf(y_dif_1_12,lag.max=36,plot=FALSE),ylim=c(-1,1),main="autocorrélogramme partiel")

#2- ESTIMATION
library(lmtest)

########fonction qui renvoie les modéles significatifs##########
significant_models <- function(series,
                               p_max , q_max,
                               P_max, Q_max,
                               d, D, s ,
                               z_thresh = 1.96,
                               method = "CSS-ML",
                               include.mean = FALSE) {
  
  keep_specs  <- character()
  keep_models <- list()
  
  for (p in 0:p_max)
    for (q in 0:q_max)
      for (P in 0:P_max)
        for (Q in 0:Q_max) {
          
          if (p == 0 && q == 0 && P == 0 && Q == 0 && !include.mean) next
          
          spec <- sprintf("(%d,%d,%d)(%d,%d,%d)[%d]", p,d,q, P,D,Q, s)
          
          fit <- try(
            arima(series,
                  order    = c(p, d, q),
                  seasonal = list(order = c(P, D, Q), period = s),
                  include.mean = include.mean,
                  method   = method),
            silent = TRUE)
          
          if (inherits(fit, "try-error")) next
          if (length(coef(fit)) == 0)      next
          
          z_vals <- try(coeftest(fit)[ , "z value"], silent = TRUE)
          if (inherits(z_vals, "try-error")) next
          if (anyNA(z_vals))               next
          
          if (all(abs(z_vals) > z_thresh)) {
            keep_specs  <- c(keep_specs, spec)
            keep_models[[spec]] <- fit
          }
        }
  
  if (length(keep_specs) == 0) {
    message("Aucun modèle dont tous les |z| dépassent ", z_thresh)
    return(invisible(NULL))
  } else {
    out <- data.frame(spec = keep_specs, row.names = NULL)
    attr(out, "models") <- keep_models
    return(out)
  }
}


#########Application de la fonction sur notre série##########
mods <- significant_models(train_ts,
                           p_max = 1, q_max = 2,
                           P_max = 2, Q_max = 1,
                           d = 1, D = 1, s = 12)
mods


#3- VALIDATION

####################Test de Box-Pierce#####################
modlist <- attr(mods, "models")

# fonction utilitaire--------------------------------------
bp_test <- function(fit, lag) {
  k <- length(coef(fit))
  Box.test(residuals(fit),
           lag   = lag,
           type  = "Box-Pierce",
           fitdf = k)$p.value
}


# tableau récapitulatif -----------------------------------
lag_choice <- 12
bp_pvals <- sapply(modlist, bp_test, lag = lag_choice)

specs   <- names(modlist)
p_vals  <- sapply(modlist, bp_test, lag = lag_choice)

bp_tbl  <- data.frame(
  spec = specs,
  p_BP = p_vals,
  row.names = NULL
)

print(bp_tbl, digits = 4)


## filtrage : résidus (p > 0,05) ----------------------------
ok_idx   <- which(bp_tbl$p_BP > 0.05 & !is.na(bp_tbl$p_BP))
ok_specs <- bp_tbl$spec[ok_idx]
ok_specs

############################AIC##############################
aic_vals <- sapply(modlist[ok_specs], AIC)

final_tbl <- data.frame(
  spec = ok_specs,
  p_BP = bp_tbl$p_BP[ok_idx],
  AIC  = aic_vals,
  row.names = NULL
)

## tri par AIC croissant +affichage
final_tbl <- final_tbl[order(final_tbl$AIC, decreasing = FALSE), ]
print(final_tbl, digits = 4)


#4- PREVISION
model=arima(train_ts,order=c(0,1,1),list(order=c(0,1,1),period=12),include.mean=FALSE,method="CSS-ML")

pred_model=predict(model,n.ahead=12)
plot(pred_model$pred)
real_model<-floor( as.numeric( gsub(",", "", test_df[,2]) ) )
plot(real_model,type="l", col="red" )
t_test    <- as.Date(test_df[,1])
plot(t_test, real_model, type = "l", col = "red", lwd = 2,
     ylim = range(c(real_model, pred_model$pred)),
     xlab = "Année", ylab = "Production",
     main = "Prévision des 12 derniers mois\nSARIMA (0,1,1)(0,1,1)[12]")
lines(t_test, pred_model$pred, col = "blue", lwd = 2)
legend("topleft", bty = "n",
       col   = c("red",  "blue"),
       lwd   = 2,
       legend = c("Observé", "Prévision"))

###################Evaluation des prévisions###############
rmse=sqrt(mean((real_model - pred_model$pred)^2))
rmse
mape=mean(abs(1-pred_model$pred/real_model))*100
mape

###########Vérification de la normalité des résidus###########
residus  <- residuals(model)                  # résidus bruts
std_residuals <- residus / sd(residus)        # résidus standardisés
shapiro.test(std_residuals)

##########AUTO ARIMA ############
library(forecast)
m<-milk_production[,2]
m<-floor( as.numeric( gsub(",", "", m) ) )
m.ts<-ts(m,start=1962, delta=1/12)
auto.arima(m.ts, trace=TRUE)
