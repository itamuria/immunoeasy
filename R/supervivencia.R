# setwd("G:/Mi unidad/LMontuenga/Data/DatosFirma/Datos firmas/Estadística Octubre 2016/")
# setwd("G:/Mi unidad/LMontuenga/Data/DatosFirma/Datos firmas/Papers/SCC/Estadística SCC/")

library(foreign)
library(dplyr)
library(survival)
library(survminer)
library(ggpubr)

# dir()

# data <- read.spss("MDA SCC FINAL.sav", to.data.frame=TRUE)
# data <- read.spss("MDA SCC FINAL todos los pacientes completos 020518.sav", to.data.frame=TRUE)

setwd("G:/Mi unidad/LMontuenga/Datos firmas/Datos firmas/Papers/ADC/Estadística ADC")
data <- read.spss("MDA ADC FINAL LOS PACIENTES COMPLETOS.sav", to.data.frame=TRUE)
grep("SL", names(data), value = T)

# setwd("G:/Mi unidad/LMontuenga/Datos firmas/Datos firmas/Papers/SCC/Estadística SCC/")
# data <- read.spss("MDA SCC FINAL todos los pacientes completos.sav", to.data.frame=TRUE)
# grep("SL", names(data), value = T)

head(data)

# library("survival")
# library("survminer")

# res.cox <- coxph(Surv(TimeSurvival60months, Survival) ~ Stage, data = data)
# res.cox
#
# summary(res.cox)
# cox.zph(res.cox)

#  Univariate -------------------------------------------------------------


covariates <- c("SNRPE","RRM2","GLUT1MB",
                "GLUT1CIT","RAD51CIT","RAD51NUC",
                "QKINUC","QKICIT","SIRT2CIT",
                "SIRT2NUC","LIG1","SRSF1",
                "RAE1CIT","CDC6NUC","CDC6CIT",
                "STC1NUC","STC1CIT","BRCA1")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(TimeSurvival60months, Survival) ~ Stage + ', x)))

form <- as.formula("Surv(TimeSurvival60months, Survival) ~ Stage + RRM2")
su <- coxph(form, data = data)
summary(su)
cox.zph(su)


univ_models <- lapply( univ_formulas, function(x){coxph(x, data = data)})




# Extract data
univ_results <- lapply(univ_models,
                       function(x){
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)

                         prot <- rownames(x$coefficients)[2]

                         pvalue<-signif(x$coefficients[1,5], digits=2);#coeficient beta
                         beta<-signif(x$coef[1,1], digits=2);#coeficient beta
                         HR <-signif(x$coef[1,2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[1,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[1,"upper .95"],2)
                         HR <- paste0(HR, " (",
                                      HR.confint.lower, "-", HR.confint.upper, ")")

                         pvalue2<-signif(x$coefficients[2,5], digits=2);#coeficient beta
                         beta2<-signif(x$coef[2,1], digits=2);#coeficient beta
                         HR2 <-signif(x$coef[2,2], digits=2);#exp(beta)
                         HR.confint.lower2 <- signif(x$conf.int[2,"lower .95"], 2)
                         HR.confint.upper2 <- signif(x$conf.int[2,"upper .95"],2)
                         HR2 <- paste0(HR2, " (",
                                       HR.confint.lower2, "-", HR.confint.upper2, ")")

                         res<-c(prot, beta, HR, pvalue, beta2, HR2, pvalue, wald.test, p.value)
                         names(res)<-c("Protein", "betaStage", "HR stage (95% CI for HR)", "pvalue Stage", "betaProt", "HR Prot (95% CI for HR)", "pvalue Prot","wald.test",
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
# res <- t(as.data.frame(univ_results, check.names = FALSE))
# as.data.frame(res)

dat <- data.frame(matrix(unname(unlist(univ_results)), byrow = TRUE, ncol = 9))
names(dat) <- c("Proteina","betaStage", "HR stage (95% CI for HR)", "pvalue-Stage", "betaProt", "HR Prot (95% CI for HR)", "pvalue-Prot","wald.test",
                "p.value")
dat <- dat[order(dat$p.value),]
# sapply(dat, class)

for(d in c(2,4,5,7:9))
{
  dat[,d] <- as.numeric(as.character(dat[,d]))
}


# Univariate sin ajustar por stage ----------------------------------------

covariates <- c("SNRPE","RRM2","GLUT1MB",
                "GLUT1CIT","RAD51CIT","RAD51NUC",
                "QKINUC","QKICIT","SIRT2CIT",
                "SIRT2NUC","LIG1","SRSF1",
                "RAE1CIT","CDC6NUC","CDC6CIT",
                "STC1NUC","STC1CIT","BRCA1")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(TimeSurvival60months, Survival) ~ ', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = data)})


# Extract data
univ_results <- lapply(univ_models,
                       function(x){
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (",
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         formula_extract <- paste( c(x$formula[[2]],
                                                     x$formula[[1]],
                                                     x$formula[[3]]), collapse='')

                         res<-c(beta, HR, wald.test, p.value, formula_extract)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test",
                                       "p.value", "formula_extract")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
# res <- t(as.data.frame(univ_results, check.names = FALSE))
# as.data.frame(res)


dat2 <- data.frame(names(univ_results), matrix(unname(unlist(univ_results)), byrow = TRUE, ncol = 5))
names(dat2) <- c("Proteina","betaStage", "HR stage (95% CI for HR)", "wald.test","p.value", "formula")
dat2 <- dat2[order(dat2$p.value),]

# sapply(dat, class)

for(d in c(2,4,5))
{
  dat2[,d] <- as.numeric(as.character(dat2[,d]))
}


# Multivariante Ibon ------------------------------------------------------

all <- c("SNRPE","RRM2","GLUT1MB",
         "GLUT1CIT","RAD51CIT","RAD51NUC",
         "QKINUC","QKICIT","SIRT2CIT",
         "SIRT2NUC","LIG1","SRSF1",
         "RAE1CIT","CDC6NUC","CDC6CIT",
         "STC1NUC","STC1CIT","BRCA1")

library(gtools)

# combinations(n=18, r=1, v=all)
# combinations(n=18, r=2, v=all)
# combinations(n=18, r=3, v=all)
# combinations(n=18, r=4, v=all)
# combinations(n=18, r=5, v=all)
# combinations(n=18, r=6, v=all)
# combinations(n=18, r=7, v=all)
#
# combinations(n=18, r=8, v=all)
# combinations(n=18, r=9, v=all)
# combinations(n=18, r=10, v=all)
#
# combinations(n=18, r=11, v=all)
# combinations(n=18, r=12, v=all)
# combinations(n=18, r=13, v=all)
# combinations(n=18, r=14, v=all)
# combinations(n=18, r=15, v=all)
# combinations(n=18, r=16, v=all)
# combinations(n=18, r=17, v=all)
# combinations(n=18, r=18, v=all)

data_fil <- data[,1:60]
mi_cindex <- function (formula, data)
{
  cox1 <- coxph(formula, data)
  y <- cox1[["y"]]
  p <- ncol(y)
  time <- y[, p - 1]
  status <- y[, p]
  x <- cox1$linear.predictors
  n <- length(time)
  ord <- order(time, -status)
  time <- time[ord]
  status <- status[ord]
  x <- x[ord]
  wh <- which(status == 1)
  if(wh[length(wh)]==length(time)) wh <- wh[-length(wh)]
  total <- concordant <- 0
  for (i in wh) {
    # print(paste0("i:",i))
    for (j in ((i + 1):n)) {
      # print(paste0("i:",i,"j:",j))
      if (time[j] > time[i]) {
        total <- total + 1
        if (x[j] < x[i])
          concordant <- concordant + 1
        if (x[j] == x[i])
          concordant <- concordant + 0.5
      }
    }
  }
  return(list(concordant = concordant, total = total, cindex = concordant/total))
}


whole_prot_analysis <- function(x, datuak){

  print(x$formula)
  kk <- summary(x)
  k3 <- c(kk$n,kk$nevent,
          kk$rsq,
          kk$concordance,
          kk$logtest,
          kk$waldtest ,
          kk$sctest,
          kk$robscore)

  # prot <- names(x)
  cc <- cox.zph(x)
  tt <- cc$table

  # PI general

  datuak$PI <- datuak$Stage * x$coefficients[1]

  for(v in 2:length(x$coefficients))
  {
    datuak$PI <- datuak$PI + datuak[,names(x$coefficients)[v]] * x$coefficients[v]
    # probar con suma scalada
  }

  # PI <- 0.335946 * data$Stage + 0.003686 * data$RRM2 - 0.001029 * data$GLUT1CIT + 0.004553 * data$QKICIT - 0.002510  * data$SIRT2CIT

  # data$PI <- PI

  library(dynpred)
  # coxFormula(univ_models)
  formule <- as.formula(paste0("Surv(TimeSurvival60months, Survival) ~ PI"))
  c_calc <- cindex(formula = formule, datuak)
  # ggcoxzph(cox.zph(k))

  cind_all <- c(c_calc$cindex,c_calc$cindex)

  # http://www.sthda.com/english/wiki/cox-model-assumptions

  # ggcoxdiagnostics(fit = k, type = , linear.predictions = TRUE)
  #
  # ggcoxdiagnostics(res.cox, type = "dfbeta",
  #                  linear.predictions = FALSE, ggtheme = theme_bw())
  #
  # ggcoxdiagnostics(res.cox, type = "deviance",
  #                  linear.predictions = FALSE, ggtheme = theme_bw())

  # Shrinkage

  # datuak <- data_fil
  print(Sys.time())
  for(bb in 1:100)
  {
    # print(paste0(bb, "/100"))
    tf_Sample <- sample(1:nrow(datuak), replace = TRUE)
    datuak_boot <- datuak[tf_Sample,]

    # Coefficienteak hemendik atera
    cox_model_run <- coxph(x$formula, data = datuak_boot)

    # print(paste0("##########################  after cox_model_run"))
    #### C-VALIDATED

    datuak_boot$PI <- datuak_boot$Stage * cox_model_run$coefficients[1]
    # print(mean(datuak_boot$PI, na.rm = TRUE))

    for(v in 2:length(cox_model_run$coefficients))
    {
      # print(v)
      datuak_boot$PI <- datuak_boot$PI + datuak_boot[,names(cox_model_run$coefficients)[v]] * cox_model_run$coefficients[v]
      # probar con suma scalada
    }
    # print(formule)
    # print(dim(datuak_boot))
    datuak_boot2 <- datuak_boot[,c("TimeSurvival60months","Survival","PI")]
    datuak_boot2 <- datuak_boot2[complete.cases(datuak_boot2),]
    # print(datuak_boot2)

    tryCatch(
      {
        c_calc2 <- cindex(formula = formule, datuak_boot)
      },error=function(cond) {
        message("Aqui hay un error:")
        # Choose a return value in case of error
        c_calc2 <- mi_cindex(formula = formule, datuak_boot)

      })

    # print(paste0("##########################  after c_calc2"))
    # print(c_calc2$cindex)
    #### C-VALIDATED

    datuak_real <- datuak

    datuak_real$PI2 <- datuak_real$Stage * cox_model_run$coefficients[1]

    for(v in 2:length(cox_model_run$coefficients))
    {
      # datuak_real[,names(cox_model_run$coefficients)[v]] <- datuak_real[,names(cox_model_run$coefficients)[v]] * cox_model_run$coefficients[v]
      datuak_real$PI2 <- datuak_real$PI2 + datuak_real[,names(cox_model_run$coefficients)[v]] * cox_model_run$coefficients[v]
      # probar con suma scalada
    }

    formule2 <- as.formula(paste0("Surv(TimeSurvival60months, Survival) ~ PI2"))
    c_calc3 <- cindex(formula = formule2, datuak_real)
    # print(paste0("##########################  after c_calc3"))
    # print(c_calc3$cindex)
    cind_all <- c(cind_all, c_calc2$cindex, c_calc3$cindex)
  }
  # print("100 rounds")
  # print(Sys.time())
  df_cind_all <- data.frame(matrix(cind_all, byrow = TRUE, ncol = 2))
  df_cind_all$Dif <- df_cind_all$X1 - df_cind_all$X2

  error <- mean(df_cind_all$Dif[-1], na.rm = TRUE)

  valor_final <- c_calc$cindex - error
  print(paste0("Final: ",valor_final))

  res<-c(as.character(x$formula)[3],k3,tt[rownames(tt) == "GLOBAL",],c_calc$concordant, c_calc$total, c_calc$cindex, valor_final)
  return(res)

}
library(riskRegression)

kanpos <- c("formula","n","n_events","rsq","maxrsq", "C-concord", "se-concord", "logtest", "logtest_df", "logtest_pvalue","waldtest",
            "wald_df","wald_pvalue","sctest","sctest_df","sctest_pvalue", "Global_chi2","Global_df","Global_pvalues",
            "cindex_concord","cindex_total","cindex_cindex","Cindex_ajustado","NumMix")

df_final_total <- data.frame(matrix(NA,ncol = length(kanpos)))
names(df_final_total) <- kanpos

konb3 <- NULL
for(ww in 2:9) #1:9
{
  print(paste0("###################################################################  ", ww, "      ##########"))
  print(Sys.time())
  # combinations(n=18, r=ww, v=all)

  konb <- data.frame(combinations(n=18, r=ww, v=all))
  konb2 <- apply(konb, 1, function(x) paste0(x,collapse = "+"))

  konb3 <- c(konb3, konb2)
}
# openxlsx::write.xlsx(konb2, file = "../../../Konbinations_4.xlsx")

# all_formulas_temp <- sapply(konb2,
#                        function(x) as.formula(paste('Surv(TimeSurvival60months, Survival) ~ Stage + ', x)))
#
# all_formulas <- list(all_formulas, all_formulas_temp)
#
# print(paste0("-------------------------------------------------------------------- all formulas done      --------------------------------------------------------------------"))
# }

t1 <- Sys.time()

p1 <- 1
stepp <- 50
tround <- Sys.time()
tm <- Sys.time()

for(sek in 1:300) #1:9
{

  ronda <- tround - tm
  tm <- Sys.time()

  tdif <- tm - t1
  hhours <- ((tdif * length(konb3) / sek) - tdif) / 60 / 60



  print(paste0("###################################################################  ", round(sek), " - faltan ", round(hhours),"   horas - Ronda: ", round(ronda), " ##########"))
  print(Sys.time())
  # combinations(n=18, r=ww, v=all)

  # konb <- data.frame(combinations(n=18, r=ww, v=all))
  # konb2 <- apply(konb, 1, function(x) paste0(x,collapse = "+"))

  # openxlsx::write.xlsx(konb2, file = "../../../Konbinations_4.xlsx")

  p2 <- p1 + stepp - 1

  all_formulas <- sapply(konb3[p1:p2],
                         function(x) as.formula(paste('Surv(TimeSurvival60months, Survival) ~ Stage + ', x)))

  # print(paste0("-------------------------------------------------------------------- all formulas done      --------------------------------------------------------------------"))
  # print(Sys.time())
  all_models <- lapply( all_formulas, function(x){coxph(x, data = data_fil)})

  # print(paste0("-------------------------------------------------------------------- all models done      --------------------------------------------------------------------"))
  # print(Sys.time())
  # Extract data
  all_results2 <- lapply(all_models,whole_prot_analysis, datuak = data_fil)

  print(paste0("-------------------------------------------------------------------- all extraction done  sek: ", sek, "    --------------------------------------------------------------------"))
  print(Sys.time())
  dat <- data.frame(matrix(unname(unlist(all_results2)), byrow = TRUE, ncol = 23))
  names(dat) <- c("formula","n","n_events","rsq","maxrsq", "C-concord", "se-concord", "logtest", "logtest_df", "logtest_pvalue","waldtest",
                  "wald_df","wald_pvalue","sctest","sctest_df","sctest_pvalue", "Global_chi2","Global_df","Global_pvalues",
                  "cindex_concord","cindex_total","cindex_cindex","Cindex_ajustado")

  dat <- dat[order(dat$Cindex_ajustado, decreasing = TRUE),]

  for(d in c(2:23))
  {
    dat[,d] <- as.numeric(as.character(dat[,d]))
  }
  dat$NumMix <- ww
  openxlsx::write.xlsx(dat, file = paste0("../../../Results_Ibon/scc_model_mix",sek,"_sep.xlsx"))

  # print(paste0("-------------------------------------------------------------------- all excel done      --------------------------------------------------------------------"))
  # print(Sys.time())
  df_final_total <- rbind(df_final_total, dat)

  # print(paste0("-------------------------------------------------------------------- all rbind done      --------------------------------------------------------------------"))
  # print(Sys.time())

  p1 <- p1 + stepp
  tround <- Sys.time()
}


t2 <- Sys.time()
print(t1)
print(t2)
print(t2-t1)

openxlsx::write.xlsx(df_final_total, file = paste0("../../../scc_model_df_final_total.xlsx"))



# form <- as.formula("Surv(TimeSurvival60months, Survival) ~ Stage + RRM2")
# su <- coxph(form, data = data)


all_formulas[1]


# res.cox <- coxph(Surv(TimeSurvival60months, Survival) ~ Stage + RRM2+GLUT1CIT+QKICIT+SIRT2CIT , data = data)
# res.cox
#
# summary(res.cox)
# cox.zph(res.cox)

h <- all_models[1]
names(h)
# data.frame(h)

# mod0 <- broom::tidy(h[[1]])

# pvalue.coxph(h[[1]])
# res.cox.sum <- summary(h[[1]])$coefficients
# summary(h[[1]])




# myfxn <- function(var1,var2,var3){
#   var1*var2*var3
#
# }
#
# lapply(1:3,myfxn,var2=2,var3=100)




# for1 <- as.formula(paste0("Surv(TimeSurvival60months, Survival) ~ Stage + CDC6CIT + GLUT1MB"))
# co <- coxph(for1, data = data_fil)
#
# gg <- openxlsx::read.xlsx("G:/Mi unidad/LMontuenga/Libro1.xlsx")
# gg2 <- gg[order(gg$TimeSurvival60months),]
# gg2$PI <- ifelse(gg2$PI < 0, abs(gg2$PI), gg2$PI)
# cindex(formula = formule, gg2)


# res <- t(as.data.frame(univ_results2, check.names = FALSE))
# as.data.frame(res)
















# Multivariate ------------------------------------------------------------

all <- c("SNRPE","RRM2","GLUT1MB",
         "GLUT1CIT","RAD51CIT","RAD51NUC",
         "QKINUC","QKICIT","SIRT2CIT",
         "SIRT2NUC","LIG1","SRSF1",
         "RAE1CIT","CDC6NUC","CDC6CIT",
         "STC1NUC","STC1CIT","BRCA1")

# All to include

res.cox <- coxph(Surv(TimeSurvival60months, Survival) ~ Stage + SNRPE + STC1CIT+ STC1NUC+RAE1CIT+SRSF1+RRM2+GLUT1CIT+GLUT1MB+RAD51CIT+RAD51NUC+QKINUC+QKICIT+SIRT2CIT+SIRT2NUC+LIG1+CDC6NUC+CDC6CIT+BRCA1, data =  data)
summary(res.cox)

ss <- survfit(res.cox)
cox.zph(res.cox)
ggsurvplot(data = data, ss, color = "#2E9FDF",
           ggtheme = theme_minimal())



# Kepler-meier

# http://www.sthda.com/english/wiki/cox-proportional-hazards-model

library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)

setwd("G:/Mi unidad/LMontuenga")
dato <- openxlsx::read.xlsx("Results/20210215_merged_PI.xlsx")

plot(dato$Cindex_ajustado)
summary(dato$Cindex_ajustado)
hist(dato$Cindex_ajustado, breaks = 50)

# For each quantity of predictors max adjusted C-index

library(stringr)
dato$num_var <- str_count(dato$formula, "\\+")

dborra <- data.frame(dato$formula, dato$num_var)

library(dplyr)
dato %>%
  group_by(num_var) %>% summarise(max_cindex = max(Cindex_ajustado),
                                  max_chi = max(Global_chi2)) -> resumen

plot(resumen$max_cindex)
plot(resumen$max_chi)
plot(resumen$max_chi, resumen$max_cindex)

# Continue

km <- with(data, Surv(TimeSurvival60months, Survival))
head(km,80)

km_fit <- survfit(Surv(TimeSurvival60months, Survival) ~ 1, data=data)
summary(km_fit, times = c(1,25,50,75,100,150))

autoplot(km_fit)

mPI <- median(data$PI, na.rm = TRUE)

data$PI2 <- ifelse(data$PI <= mPI, "Low","Up")

km_trt_fit <- survfit(Surv(TimeSurvival60months, Survival) ~ PI2, data=data)
autoplot(km_trt_fit, conf.int = FALSE)
# https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/

library(survminer)

ggsurvplot(data,
           fit = survfit(Surv(TimeSurvival60months, Survival) ~ 1, data),
           xlab = "Days",
           ylab = "Overall survival probability")

ggsurvplot(data,
           fit = survfit(Surv(TimeSurvival60months, Survival) ~ PI2, data),
           xlab = "Days",
           ylab = "Overall survival probability")

ggsurvplot(fit = survfit(Surv(TimeSurvival60months, Survival) ~ PI2, data),
           data = data,
           surv.median.line = "hv", # Add medians survival

           # Change legends: title & labels
           legend.title = "PI",
           legend.labs = c("Low", "High"),
           # Add p-value and tervals
           pval = TRUE,

           conf.int = TRUE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),

           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("#E7B800", "#2E9FDF"),
           ggtheme = theme_bw() # Change ggplot2 theme
)

summary(survfit(Surv(TimeSurvival60months, Survival) ~ 1, data = data), times = 10)
survdiff(Surv(TimeSurvival60months, Survival) ~ PI2, data = data)

library(gtsummary)
coxph(Surv(TimeSurvival60months, Survival) ~ PI2, data = data) %>%
  gtsummary::tbl_regression(exp = TRUE)




# res.cox <- coxph(Surv(TimeSurvival60months, Survival) ~ Stage + RRM2 +  GLUT1MB + SNRPE + RAD51CIT, data =  data)
res.cox <- coxph(Surv(TimeSurvival60months, Survival) ~ Stage + STC1CIT+ STC1NUC+RAE1CIT+SRSF1+RRM2+GLUT1CIT, data =  data)
summary(res.cox)

ss <- survfit(res.cox)
cox.zph(res.cox)
ggsurvplot(data = data, ss, color = "#2E9FDF",
           ggtheme = theme_minimal())

# Create the new data
new_df <- with(data,
               data.frame(Stage = c(1, 2,3),
                          RAD51CIT = rep(mean(RAD51CIT, na.rm = TRUE), 3)
               )
)
new_df

# Survival curves
fit <- survfit(res.cox, newdata = new_df)
ggsurvplot(data,fit, conf.int = TRUE, legend.labs=c("Stage=1", "Stage=2","Stage=3"),
           ggtheme = theme_minimal())



# Frogak ------------------------------------------------------------------

setwd("G:/Mi unidad/LMontuenga/Data/DatosFirma/Datos firmas/Papers/SCC/Estadística SCC/Bootstrapping/Después del análisis para OS/")

# Todos los bootstraps juntos (ona)
data <- read.spss("BootdataSCC.sav", to.data.frame=TRUE)
#
data <- read.spss("OriginalSCC.sav", to.data.frame=TRUE)
#
data <- read.spss("XBetas.sav", to.data.frame=TRUE)


# Todos los bootstraps juntos (ona)
data <- read.spss("../OriginalData.sav", to.data.frame=TRUE)
#
# data <- read.spss("OriginalSCC.sav", to.data.frame=TRUE)
# #
# data <- read.spss("XBetas.sav", to.data.frame=TRUE)


# Protocol - 1 Until 20 ---------------------------------------------------

setwd("G:/Mi unidad/LMontuenga/Data/DatosFirma/Datos firmas/Papers/SCC/Estadística SCC/")

library(foreign)
library(dplyr)
library(survival)
library(survminer)
library(ggpubr)

dir()

data <- read.spss("MDA SCC FINAL todos los pacientes completos 020518.sav", to.data.frame=TRUE)


res.cox <- coxph(Surv(TimeSurvival60months, Survival) ~ Stage + KRASMUT+EGFRMUT+SNRPE+RRM2+GLUT1MB+GLUT1CIT+RAD51CIT+RAD51NUC+QKINUC+QKICIT+SIRT2CIT+SIRT2NUC+LIG1+SRSF1+RAE1CIT+CDC6NUC+CDC6CIT+STC1NUC+STC1CIT+BRCA1 , data = data)
res.cox <- coxph(Surv(TimeSurvival60months, Survival) ~ Stage + RRM2+GLUT1CIT+QKICIT+SIRT2CIT , data = data)
res.cox

summary(res.cox)
cox.zph(res.cox)


# 2-Calculo PI ------------------------------------------------------------

PI <- 0.335946 * data$Stage + 0.003686 * data$RRM2 - 0.001029 * data$GLUT1CIT + 0.004553 * data$QKICIT - 0.002510  * data$SIRT2CIT

data$PI <- PI

library(dynpred)
cindex(formula =Surv(TimeSurvival60months, Survival) ~ Stage + RRM2+GLUT1CIT+QKICIT+SIRT2CIT, data)
ggcoxzph(cox.zph(res.cox))

# http://www.sthda.com/english/wiki/cox-model-assumptions

ggcoxdiagnostics(fit = res.cox, type = , linear.predictions = TRUE)

ggcoxdiagnostics(res.cox, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())

ggcoxdiagnostics(res.cox, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_bw())


# Kepler-meier

# http://www.sthda.com/english/wiki/cox-proportional-hazards-model

library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)

km <- with(data, Surv(TimeSurvival60months, Survival))
head(km,80)

km_fit <- survfit(Surv(TimeSurvival60months, Survival) ~ 1, data=data)
summary(km_fit, times = c(1,25,50,75,100,150))

autoplot(km_fit)

mPI <- median(data$PI, na.rm = TRUE)

data$PI2 <- ifelse(data$PI <= mPI, "Low","Up")

km_trt_fit <- survfit(Surv(TimeSurvival60months, Survival) ~ PI2, data=data)
autoplot(km_trt_fit, conf.int = FALSE)
# https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/

library(survminer)

ggsurvplot(data,
           fit = survfit(Surv(TimeSurvival60months, Survival) ~ 1, data),
           xlab = "Days",
           ylab = "Overall survival probability")

ggsurvplot(data,
           fit = survfit(Surv(TimeSurvival60months, Survival) ~ PI2, data),
           xlab = "Days",
           ylab = "Overall survival probability")

ggsurvplot(fit = survfit(Surv(TimeSurvival60months, Survival) ~ PI2, data),
           data = data,
           surv.median.line = "hv", # Add medians survival

           # Change legends: title & labels
           legend.title = "PI",
           legend.labs = c("Low", "High"),
           # Add p-value and tervals
           pval = TRUE,

           conf.int = TRUE,
           # Add risk table
           risk.table = TRUE,
           tables.height = 0.2,
           tables.theme = theme_cleantable(),

           # Color palettes. Use custom color: c("#E7B800", "#2E9FDF"),
           # or brewer color (e.g.: "Dark2"), or ggsci color (e.g.: "jco")
           palette = c("#E7B800", "#2E9FDF"),
           ggtheme = theme_bw() # Change ggplot2 theme
)

summary(survfit(Surv(TimeSurvival60months, Survival) ~ 1, data = data), times = 10)
survdiff(Surv(TimeSurvival60months, Survival) ~ PI2, data = data)

library(gtsummary)
coxph(Surv(TimeSurvival60months, Survival) ~ PI2, data = data) %>%
  gtsummary::tbl_regression(exp = TRUE)

# https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html

# https://daviddalpiaz.github.io/r4sl/shrinkage-methods.html
# https://bookdown.org/egarpor/PM-UC3M/lm-iii.html



# Comparando --------------------------------------------------------------

res.cox <- coxph(Surv(TimeSurvival60months, Survival) ~ RRM2+RAD51CIT+CDC6NUC+SIRT2CIT+RAD51NUC+BRCA1+GLUT1MB , data = data)
res.cox

summary(res.cox)
cox.zph(res.cox)
