## function for the model 1 binary ##
uni.fit1 <- function(model,var) {
  s <-  summary(model)
  rownam <- rownames(s$coefficients)
  s <-  s$coefficients[-1,]
  sand <- sqrt(diag(sandwich::vcovCL(model,cluster = var)))
  r <- as.data.frame(list(exp(s[1]), exp(s[1] - 1.96*sand[2]),exp(s[1] + 1.96*sand[2]),round(2*(1-pnorm(abs(s[1]/sand[2]))),digits =5)))
  
  colnames(r) <- c( "RR","Lower CI", "Upper CI","p-value")
  r$RR <- round(r$RR,digits =2)
  r$`Lower CI` <- round(r$`Lower CI`,digits=3)
  r$`Upper CI` <- round(r$`Upper CI`,digits =3)
  rownames(r) <- rownam[2]
  return(r)
}

uni.fit.conti1 <- function(model,var,incr) {
  s <-  summary(model)
  rownam <- rownames(s$coefficients)
  s <-  s$coefficients[-1,]
  sand <- sqrt(diag(sandwich::vcovCL(model,cluster = var)))
  r <- as.data.frame(list(exp(s[1]*incr), exp((s[1] - 1.96*sand[2])*incr),exp((s[1] + 1.96*sand[2])*incr),round(2*(1-pnorm(abs(s[1]/sand[2]))),digits =5)))
  colnames(r) <- c( "RR","Lower CI", "Upper CI","p-value")
  r$RR <- round(r$RR,digits =2)
  r$`Lower CI` <- round(r$`Lower CI`,digits=3)
  r$`Upper CI` <- round(r$`Upper CI`,digits =3)
  rownames(r) <- rownam[2]
  return(r)
}

## function for the model 1 binary ##
uni.fit <- function(model) {
  s <-  summary(model)
  rownam <- rownames(s$coefficients)
  s <-  s$coefficients[-1,]
  sand <- sqrt(diag(sandwich(model)))
  r <- as.data.frame(list(exp(s[1]), exp(s[1] - 1.96*sand[2]),exp(s[1] + 1.96*sand[2]),round(2*(1-pnorm(abs(s[1]/sand[2]))),digits =5)))
  
  colnames(r) <- c( "RR","Lower CI", "Upper CI","p-value")
  r$RR <- round(r$RR,digits =2)
  r$`Lower CI` <- round(r$`Lower CI`,digits=3)
  r$`Upper CI` <- round(r$`Upper CI`,digits =3)
  rownames(r) <- rownam[2]
  return(r)
}

uni.fit.conti <- function(model,incr) {
  s <-  summary(model)
  rownam <- rownames(s$coefficients)
  s <-  s$coefficients[-1,]
  sand <- sqrt(diag(sandwich(model)))
  r <- as.data.frame(list(exp(s[1]*incr), exp((s[1] - 1.96*sand[2])*incr),exp((s[1] + 1.96*sand[2])*incr),round(2*(1-pnorm(abs(s[1]/sand[2]))),digits =5)))
  
  colnames(r) <- c( "RR","Lower CI", "Upper CI","p-value")
  r$RR <- round(r$RR,digits =2)
  r$`Lower CI` <- round(r$`Lower CI`,digits=3)
  r$`Upper CI` <- round(r$`Upper CI`,digits =3)
  rownames(r) <- rownam[2]
  return(r)
}
multi.fit.conti <- function(model) {
  s <-  summary(model)
  rownam <- rownames(s$coefficients)
  s <-  s$coefficients[-1,]
  sand <- sqrt(diag(sandwich(model)))
  r <- rbind(c(exp(s[1]*10), exp((s[1] - 1.96*sand[2])*10),exp((s[1] + 1.96*sand[2])*10),round(2*(1-pnorm(abs(s[1]/sand[2]))),digits =5)),
             c(exp(s[2]*50), exp((s[2] - 1.96*sand[3])*50),exp((s[2] + 1.96*sand[3])*50),round(2*(1-pnorm(abs(s[2]/sand[3]))),digits =5)))
  r <- as.data.frame(r)
  colnames(r) <- c( "RR","Lower CI", "Upper CI","p-value")
  r$RR <- round(r$RR,digits =2)
  r$`Lower CI` <- round(r$`Lower CI`,digits=3)
  r$`Upper CI` <- round(r$`Upper CI`,digits =3)
  rownam <- rownames(s)
  rownames(r) <- rownam[1:2]
  return(r)
}
multi.fit.conti1 <- function(model,var) {
  s <-  summary(model)
  s <-  s$coefficients[-1,]
  sand <- sqrt(diag(sandwich::vcovHC(model,cluster=var)))
  r = matrix(0,nrow = (length(coef(model)) - 1), ncol = 4)
  r <- rbind(c(exp(s[1]*10), exp((s[1] - 1.96*sand[2])*10),exp((s[1] + 1.96*sand[2])*10),round(2*(1-pnorm(abs(s[1]/sand[2]))),digits =5)),
             c(exp(s[2]*50), exp((s[2] - 1.96*sand[3])*50),exp((s[2] + 1.96*sand[3])*50),round(2*(1-pnorm(abs(s[2]/sand[3]))),digits =5)),
             c(exp(s[3]*4), exp((s[3] - 1.96*sand[3])*4),exp((s[3] + 1.96*sand[4])*4),round(2*(1-pnorm(abs(s[3]/sand[4]))),digits =5)))
  r <- as.data.frame(r)
  colnames(r) <- c( "RR","Lower CI", "Upper CI","p-value")
  r$RR <- round(r$RR,digits =2)
  r$`Lower CI` <- round(r$`Lower CI`,digits=3)
  r$`Upper CI` <- round(r$`Upper CI`,digits =3)
  rownam <- rownames(s)
  rownames(r) <- rownam[1:3]
  return(r)
}
multi.fit.conti2 <- function(model) {
  s <-  summary(model)
  rownam <- rownames(s$coefficients)
  s <-  s$coefficients[-1,]
  sand <- sqrt(diag(sandwich(model)))
  r <- rbind(c(exp(s[1]*10), exp((s[1] - 1.96*sand[2])*10),exp((s[1] + 1.96*sand[2])*10),round(2*(1-pnorm(abs(s[1]/sand[2]))),digits =5)),
             c(exp(s[2]*50), exp((s[2] - 1.96*sand[3])*50),exp((s[2] + 1.96*sand[3])*50),round(2*(1-pnorm(abs(s[2]/sand[3]))),digits =5)),
             c(exp(s[3]*4), exp((s[3] - 1.96*sand[3])*4),exp((s[3] + 1.96*sand[4])*4),round(2*(1-pnorm(abs(s[3]/sand[4]))),digits =5)),
             c(exp(s[4]*10), exp((s[4] - 1.96*sand[5])*10),exp((s[4] + 1.96*sand[5])*10),round(2*(1-pnorm(abs(s[4]/sand[5]))),digits =5)))
  r <- as.data.frame(r)
  colnames(r) <- c( "RR","Lower CI", "Upper CI","p-value")
  r$RR <- round(r$RR,digits =2)
  r$`Lower CI` <- round(r$`Lower CI`,digits=3)
  r$`Upper CI` <- round(r$`Upper CI`,digits =3)
  rownam <- rownames(s)
  rownames(r) <- rownam[1:4]
  return(r)
}

uni.fit.gee <- function(model) {
  s <-  summary(model)
  s <-  s$coefficients[-1,]
  sand <- sqrt(diag(model$robust.variance))
  r = matrix(0,nrow = (length(coef(model)) - 1), ncol = 4)
  for (i in 1:(length(coef(model)) - 1)){
    r[i,] <- c(exp(s[i,1]), exp(s[i,1] - 1.96*sand[i+1]),exp(s[i,1] + 1.96*sand[i+1]),round((1 - pnorm(abs(s[i,1]/sand[i+1]))) * 2,digits=5))
  }
  
  r <- as.data.frame(r)
  colnames(r) <- c( "RR","Lower CI", "Upper CI","p-value") 
  rownam <- rownames(s)
  rownames(r) <- rownam
  return(r)
}



## function for the model 1 long ##
uni.fit2 <- function(model,var) {
  s <-  summary(model)
  s <-  s$coefficients[-1,]
  sand <- sqrt(diag(sandwich::vcovHC(model,cluster=var)))
  r = matrix(0,nrow = (length(coef(model)) - 1), ncol = 4)
  for (i in 1:(length(coef(model)) - 1)){
    r[i,] <- c(exp(s[i,1]), exp(s[i,1] - 1.96*sand[i+1]),exp(s[i,1] + 1.96*sand[i+1]),round((1 - pnorm(abs(s[i,1]/sand[i+1]))) * 2,digits=5))
  }
  
  r <- as.data.frame(r)
  colnames(r) <- c( "RR","Lower CI", "Upper CI","p-value") 
  rownam <- rownames(s)
  rownames(r) <- rownam
  return(r)
}


uni.fit3 <- function(model) {
  s <-  summary(model)
  s <-  s$coefficients[-1,]
  sand <- sqrt(diag(sandwich(model)))
  r = matrix(0,nrow = (length(coef(model)) - 1), ncol = 4)
  for (i in 1:(length(coef(model)) - 1)){
    r[i,] <- c(exp(s[i,1]), exp(s[i,1] - 1.96*sand[i+1]),exp(s[i,1] + 1.96*sand[i+1]),round((1 - pnorm(abs(s[i,1]/sand[i+1]))) * 2,digits=5))
  }
  
  r <- as.data.frame(r)
  colnames(r) <- c( "RR","Lower CI", "Upper CI","p-value") 
  rownam <- rownames(s)
  rownames(r) <- rownam
  return(r)
}


rr <- function(myage,model){
  exp(model$Coefs[model$Variables %in% "ns(age_b, df = 4)1"]*
        (sp[rownames(sp) %in% myage, 1]-sp[rownames(sp) %in% 35, 1]) +
        model$Coefs[model$Variables %in% "ns(age_b, df = 4)2"]*
        (sp[rownames(sp) %in% myage, 2]- sp[rownames(sp) %in% 35, 2]) +
        model$Coefs[model$Variables %in% "ns(age_b, df = 4)3"]*
        (sp[rownames(sp) %in% myage, 3]- sp[rownames(sp) %in% 35, 3]) +
        model$Coefs[model$Variables %in% "ns(age_b, df = 4)4"]*
        (sp[rownames(sp) %in% myage, 4]- sp[rownames(sp) %in% 35, 4]))
  
}


# my function to calculate variance and SE
spline_se_fun <- function(myage){
  
  tmp <- v_spline
  a <- sp$X1[rownames(sp) %in% myage]-sp$X1[rownames(sp) %in% 35]
  b <- sp$X2[rownames(sp) %in% myage]-sp$X2[rownames(sp) %in% 35]
  c <- sp$X3[rownames(sp) %in% myage]-sp$X3[rownames(sp) %in% 35]
  d <- sp$X4[rownames(sp) %in% myage]-sp$X4[rownames(sp) %in% 35]
  
  Var_x <- tmp["ns(age_b, df = 4)1", "ns(age_b, df = 4)1"]
  Var_y <- tmp["ns(age_b, df = 4)2", "ns(age_b, df = 4)2"]
  Var_z <- tmp["ns(age_b, df = 4)3", "ns(age_b, df = 4)3"]
  Var_w <- tmp["ns(age_b, df = 4)4", "ns(age_b, df = 4)4"]
  
  Cov_xy <- tmp["ns(age_b, df = 4)1", "ns(age_b, df = 4)2"]
  Cov_xz <- tmp["ns(age_b, df = 4)1", "ns(age_b, df = 4)3"]
  Cov_xw <- tmp["ns(age_b, df = 4)1", "ns(age_b, df = 4)4"]
  Cov_yz <- tmp["ns(age_b, df = 4)2", "ns(age_b, df = 4)3"]
  Cov_yw <- tmp["ns(age_b, df = 4)2", "ns(age_b, df = 4)4"]
  Cov_zw <- tmp["ns(age_b, df = 4)3", "ns(age_b, df = 4)4"]
  
  pooled_var_splines_yr <- a^2*Var_x + b^2*Var_y + c^2*Var_z + d^2*Var_w +
    2*a*b*Cov_xy + 2*a*c*Cov_xz + 2*a*d*Cov_xw + 2*b*c*Cov_yz + 2*b*d*Cov_yw + 2*c*d*Cov_zw
  pooled_se_splines_yr <- sqrt(pooled_var_splines_yr)
  
  # print(pooled_var_splines_yr)
  print(pooled_se_splines_yr)
  
}


abs_val <- function(myage,model){
  (model$Coefs[model$Variables %in% "ns(age_b, df = 4)1"]*
     (sp[rownames(sp) %in% myage, 1]-sp[rownames(sp) %in% 35, 1]) +
     model$Coefs[model$Variables %in% "ns(age_b, df = 4)2"]*
     (sp[rownames(sp) %in% myage, 2]- sp[rownames(sp) %in% 35, 2]) +
     model$Coefs[model$Variables %in% "ns(age_b, df = 4)3"]*
     (sp[rownames(sp) %in% myage, 3]- sp[rownames(sp) %in% 35, 3]) +
     model$Coefs[model$Variables %in% "ns(age_b, df = 4)4"]*
     (sp[rownames(sp) %in% myage, 4]- sp[rownames(sp) %in% 35, 4]))
  
}

CI <- function(rr,SE) 
{ exp(log(rr) + c(-1,1)*1.96*SE)}


rr_cd4 <- function(myage,model){
  exp(model$Coefs[model$Variables %in% "ns(cd4_b, df = 4)1"]*
        (sp_cd4[rownames(sp_cd4) %in% myage, 1]-sp_cd4[rownames(sp_cd4) %in% 1000, 1]) +
        model$Coefs[model$Variables %in% "ns(cd4_b, df = 4)2"]*
        (sp_cd4[rownames(sp_cd4) %in% myage, 2]- sp_cd4[rownames(sp_cd4) %in% 1000, 2]) +
        model$Coefs[model$Variables %in% "ns(cd4_b, df = 4)3"]*
        (sp_cd4[rownames(sp_cd4) %in% myage, 3]- sp_cd4[rownames(sp_cd4) %in% 1000, 3]) +
        model$Coefs[model$Variables %in% "ns(cd4_b, df = 4)4"]*
        (sp_cd4[rownames(sp_cd4) %in% myage, 4]- sp_cd4[rownames(sp_cd4) %in% 1000, 4]))
  
}


# my function to calculate variance and SE
spline_se_fun_cd <- function(myage){
  
  tmp <- v_spline
  a <- sp_cd4$X1[rownames(sp_cd4) %in% myage]-sp_cd4$X1[rownames(sp_cd4) %in% 1000]
  b <- sp_cd4$X2[rownames(sp_cd4) %in% myage]-sp_cd4$X2[rownames(sp_cd4) %in% 1000]
  c <- sp_cd4$X3[rownames(sp_cd4) %in% myage]-sp_cd4$X3[rownames(sp_cd4) %in% 1000]
  d <- sp_cd4$X4[rownames(sp_cd4) %in% myage]-sp_cd4$X4[rownames(sp_cd4) %in% 1000]
  
  Var_x <- tmp["ns(cd4_b, df = 4)1", "ns(cd4_b, df = 4)1"]
  Var_y <- tmp["ns(cd4_b, df = 4)2", "ns(cd4_b, df = 4)2"]
  Var_z <- tmp["ns(cd4_b, df = 4)3", "ns(cd4_b, df = 4)3"]
  Var_w <- tmp["ns(cd4_b, df = 4)4", "ns(cd4_b, df = 4)4"]
  
  Cov_xy <- tmp["ns(cd4_b, df = 4)1", "ns(cd4_b, df = 4)2"]
  Cov_xz <- tmp["ns(cd4_b, df = 4)1", "ns(cd4_b, df = 4)3"]
  Cov_xw <- tmp["ns(cd4_b, df = 4)1", "ns(cd4_b, df = 4)4"]
  Cov_yz <- tmp["ns(cd4_b, df = 4)2", "ns(cd4_b, df = 4)3"]
  Cov_yw <- tmp["ns(cd4_b, df = 4)2", "ns(cd4_b, df = 4)4"]
  Cov_zw <- tmp["ns(cd4_b, df = 4)3", "ns(cd4_b, df = 4)4"]
  
  pooled_var_splines_yr <- a^2*Var_x + b^2*Var_y + c^2*Var_z + d^2*Var_w +
    2*a*b*Cov_xy + 2*a*c*Cov_xz + 2*a*d*Cov_xw + 2*b*c*Cov_yz + 2*b*d*Cov_yw + 2*c*d*Cov_zw
  pooled_se_splines_yr <- sqrt(pooled_var_splines_yr)
  
  # print(pooled_var_splines_yr)
  print(pooled_se_splines_yr)
  
}


abs_val_cd <- function(myage,model){
  (model$Coefs[model$Variables %in% "ns(cd4_b, df = 4)1"]*
     (sp_cd4[rownames(sp_cd4) %in% myage, 1]-sp_cd4[rownames(sp_cd4) %in% 1000, 1]) +
     model$Coefs[model$Variables %in% "ns(cd4_b, df = 4)2"]*
     (sp_cd4[rownames(sp_cd4) %in% myage, 2]- sp_cd4[rownames(sp_cd4) %in% 1000, 2]) +
     model$Coefs[model$Variables %in% "ns(cd4_b, df = 4)3"]*
     (sp_cd4[rownames(sp_cd4) %in% myage, 3]- sp_cd4[rownames(sp_cd4) %in% 1000, 3]) +
     model$Coefs[model$Variables %in% "ns(cd4_b, df = 4)4"]*
     (sp_cd4[rownames(sp_cd4) %in% myage, 4]- sp_cd4[rownames(sp_cd4) %in% 1000, 4]))
  
}

rr_cd4_tn <- function(myage,model){
  exp(model$Coefs[model$Variables %in% "ns(cd4_b, df = 4)1"]*
        (sp_cd4[rownames(sp_cd4) %in% myage, 1]-sp_cd4[rownames(sp_cd4) %in% 500, 1]) +
        model$Coefs[model$Variables %in% "ns(cd4_b, df = 4)2"]*
        (sp_cd4[rownames(sp_cd4) %in% myage, 2]- sp_cd4[rownames(sp_cd4) %in% 500, 2]) +
        model$Coefs[model$Variables %in% "ns(cd4_b, df = 4)3"]*
        (sp_cd4[rownames(sp_cd4) %in% myage, 3]- sp_cd4[rownames(sp_cd4) %in% 500, 3]) +
        model$Coefs[model$Variables %in% "ns(cd4_b, df = 4)4"]*
        (sp_cd4[rownames(sp_cd4) %in% myage, 4]- sp_cd4[rownames(sp_cd4) %in% 500, 4]))
  
}


# my function to calculate variance and SE
spline_se_fun_cd_tn <- function(myage){
  
  tmp <- v_spline
  a <- sp_cd4$X1[rownames(sp_cd4) %in% myage]-sp_cd4$X1[rownames(sp_cd4) %in% 500]
  b <- sp_cd4$X2[rownames(sp_cd4) %in% myage]-sp_cd4$X2[rownames(sp_cd4) %in% 500]
  c <- sp_cd4$X3[rownames(sp_cd4) %in% myage]-sp_cd4$X3[rownames(sp_cd4) %in% 500]
  d <- sp_cd4$X4[rownames(sp_cd4) %in% myage]-sp_cd4$X4[rownames(sp_cd4) %in% 500]
  
  Var_x <- tmp["ns(cd4_b, df = 4)1", "ns(cd4_b, df = 4)1"]
  Var_y <- tmp["ns(cd4_b, df = 4)2", "ns(cd4_b, df = 4)2"]
  Var_z <- tmp["ns(cd4_b, df = 4)3", "ns(cd4_b, df = 4)3"]
  Var_w <- tmp["ns(cd4_b, df = 4)4", "ns(cd4_b, df = 4)4"]
  
  Cov_xy <- tmp["ns(cd4_b, df = 4)1", "ns(cd4_b, df = 4)2"]
  Cov_xz <- tmp["ns(cd4_b, df = 4)1", "ns(cd4_b, df = 4)3"]
  Cov_xw <- tmp["ns(cd4_b, df = 4)1", "ns(cd4_b, df = 4)4"]
  Cov_yz <- tmp["ns(cd4_b, df = 4)2", "ns(cd4_b, df = 4)3"]
  Cov_yw <- tmp["ns(cd4_b, df = 4)2", "ns(cd4_b, df = 4)4"]
  Cov_zw <- tmp["ns(cd4_b, df = 4)3", "ns(cd4_b, df = 4)4"]
  
  pooled_var_splines_yr <- a^2*Var_x + b^2*Var_y + c^2*Var_z + d^2*Var_w +
    2*a*b*Cov_xy + 2*a*c*Cov_xz + 2*a*d*Cov_xw + 2*b*c*Cov_yz + 2*b*d*Cov_yw + 2*c*d*Cov_zw
  pooled_se_splines_yr <- sqrt(pooled_var_splines_yr)
  
  # print(pooled_var_splines_yr)
  print(pooled_se_splines_yr)
  
}


abs_val_cd_tn <- function(myage,model){
  (model$Coefs[model$Variables %in% "ns(cd4_b, df = 4)1"]*
     (sp_cd4[rownames(sp_cd4) %in% myage, 1]-sp_cd4[rownames(sp_cd4) %in% 500, 1]) +
     model$Coefs[model$Variables %in% "ns(cd4_b, df = 4)2"]*
     (sp_cd4[rownames(sp_cd4) %in% myage, 2]- sp_cd4[rownames(sp_cd4) %in% 500, 2]) +
     model$Coefs[model$Variables %in% "ns(cd4_b, df = 4)3"]*
     (sp_cd4[rownames(sp_cd4) %in% myage, 3]- sp_cd4[rownames(sp_cd4) %in% 500, 3]) +
     model$Coefs[model$Variables %in% "ns(cd4_b, df = 4)4"]*
     (sp_cd4[rownames(sp_cd4) %in% myage, 4]- sp_cd4[rownames(sp_cd4) %in% 500, 4]))
  
}




abs_sq_cd <- function(cdref,model){
  (model$Coefs[model$Variables %in% "ns(sq_cd4, df = 4)1"]*
     (sp_cd4[rownames(sp_cd4) %in% cdref, 1]-sp_cd4[rownames(sp_cd4) %in% 14.14, 1]) +
     model$Coefs[model$Variables %in% "ns(sq_cd4, df = 4)2"]*
     (sp_cd4[rownames(sp_cd4) %in% cdref, 2]- sp_cd4[rownames(sp_cd4) %in% 14.14, 2]) +
     model$Coefs[model$Variables %in% "ns(sq_cd4, df = 4)3"]*
     (sp_cd4[rownames(sp_cd4) %in% cdref, 3]- sp_cd4[rownames(sp_cd4) %in% 14.14, 3]) +
     model$Coefs[model$Variables %in% "ns(sq_cd4, df = 4)4"]*
     (sp_cd4[rownames(sp_cd4) %in% cdref, 4]- sp_cd4[rownames(sp_cd4) %in% 14.14, 4]))
  
}

rr_sq <- function(cdref,model){
  exp(model$Coefs[model$Variables %in% "ns(sq_cd4, df = 4)1"]*
        (sp_cd4[rownames(sp_cd4) %in% cdref, 1]-sp_cd4[rownames(sp_cd4) %in% 14.14, 1]) +
        model$Coefs[model$Variables %in% "ns(sq_cd4, df = 4)2"]*
        (sp_cd4[rownames(sp_cd4) %in% cdref, 2]- sp_cd4[rownames(sp_cd4) %in% 14.14, 2]) +
        model$Coefs[model$Variables %in% "ns(sq_cd4, df = 4)3"]*
        (sp_cd4[rownames(sp_cd4) %in% cdref, 3]- sp_cd4[rownames(sp_cd4) %in% 14.14, 3]) +
        model$Coefs[model$Variables %in% "ns(sq_cd4, df = 4)4"]*
        (sp_cd4[rownames(sp_cd4) %in% cdref, 4]- sp_cd4[rownames(sp_cd4) %in% 14.14, 4]))
  
}


# my function to calculate variance and SE
spline_se_fun_sq_cd <- function(cdref){
  
  tmp <- v_spline
  a <- sp_cd4$X1[rownames(sp_cd4) %in% cdref]-sp_cd4$X1[rownames(sp_cd4) %in% 14.14]
  b <- sp_cd4$X2[rownames(sp_cd4) %in% cdref]-sp_cd4$X2[rownames(sp_cd4) %in% 14.14]
  c <- sp_cd4$X3[rownames(sp_cd4) %in% cdref]-sp_cd4$X3[rownames(sp_cd4) %in% 14.14]
  d <- sp_cd4$X4[rownames(sp_cd4) %in% cdref]-sp_cd4$X4[rownames(sp_cd4) %in% 14.14]
  
  Var_x <- tmp["ns(sq_cd4, df = 4)1", "ns(sq_cd4, df = 4)1"]
  Var_y <- tmp["ns(sq_cd4, df = 4)2", "ns(sq_cd4, df = 4)2"]
  Var_z <- tmp["ns(sq_cd4, df = 4)3", "ns(sq_cd4, df = 4)3"]
  Var_w <- tmp["ns(sq_cd4, df = 4)4", "ns(sq_cd4, df = 4)4"]
  
  Cov_xy <- tmp["ns(sq_cd4, df = 4)1", "ns(sq_cd4, df = 4)2"]
  Cov_xz <- tmp["ns(sq_cd4, df = 4)1", "ns(sq_cd4, df = 4)3"]
  Cov_xw <- tmp["ns(sq_cd4, df = 4)1", "ns(sq_cd4, df = 4)4"]
  Cov_yz <- tmp["ns(sq_cd4, df = 4)2", "ns(sq_cd4, df = 4)3"]
  Cov_yw <- tmp["ns(sq_cd4, df = 4)2", "ns(sq_cd4, df = 4)4"]
  Cov_zw <- tmp["ns(sq_cd4, df = 4)3", "ns(sq_cd4, df = 4)4"]
  
  pooled_var_splines_yr <- a^2*Var_x + b^2*Var_y + c^2*Var_z + d^2*Var_w +
    2*a*b*Cov_xy + 2*a*c*Cov_xz + 2*a*d*Cov_xw + 2*b*c*Cov_yz + 2*b*d*Cov_yw + 2*c*d*Cov_zw
  pooled_se_splines_yr <- sqrt(pooled_var_splines_yr)
  
  # print(pooled_var_splines_yr)
  print(pooled_se_splines_yr)
  
}


abs_val_yrs <- function(myyear,model){
  (model$Coefs[model$Variables %in% "ns(baseline_y, df = 4)1"]*
     (sp_yr[rownames(sp_yr) %in% myyear, 1]-sp_yr[rownames(sp_yr) %in% 2010, 1]) +
     model$Coefs[model$Variables %in% "ns(baseline_y, df = 4)2"]*
     (sp_yr[rownames(sp_yr) %in% myyear, 2]- sp_yr[rownames(sp_yr) %in% 2010, 2]) +
     model$Coefs[model$Variables %in% "ns(baseline_y, df = 4)3"]*
     (sp_yr[rownames(sp_yr) %in% myyear, 3]- sp_yr[rownames(sp_yr) %in% 2010, 3]) +
     model$Coefs[model$Variables %in% "ns(baseline_y, df = 4)4"]*
     (sp_yr[rownames(sp_yr) %in% myyear, 4]- sp_yr[rownames(sp_yr) %in% 2010, 4]))
  
}

abs_val_yrs2 <- function(myyear,model){
  (model$Coefs[model$Variables %in% "ns(year, df = 4)1"]*
     (sp_yr[rownames(sp_yr) %in% myyear, 1]-sp_yr[rownames(sp_yr) %in% 2010, 1]) +
     model$Coefs[model$Variables %in% "ns(year, df = 4)2"]*
     (sp_yr[rownames(sp_yr) %in% myyear, 2]- sp_yr[rownames(sp_yr) %in% 2010, 2]) +
     model$Coefs[model$Variables %in% "ns(year, df = 4)3"]*
     (sp_yr[rownames(sp_yr) %in% myyear, 3]- sp_yr[rownames(sp_yr) %in% 2010, 3]) +
     model$Coefs[model$Variables %in% "ns(year, df = 4)4"]*
     (sp_yr[rownames(sp_yr) %in% myyear, 4]- sp_yr[rownames(sp_yr) %in% 2010, 4]))
  
}

rr_yr<- function(myyear,model){
  exp(model$Coefs[model$Variables %in% "ns(baseline_y, df = 4)1"]*
        (sp_yr[rownames(sp_yr) %in% myyear, 1]-sp_yr[rownames(sp_yr) %in% 2010, 1]) +
        model$Coefs[model$Variables %in% "ns(baseline_y, df = 4)2"]*
        (sp_yr[rownames(sp_yr) %in% myyear, 2]- sp_yr[rownames(sp_yr) %in% 2010, 2]) +
        model$Coefs[model$Variables %in% "ns(baseline_y, df = 4)3"]*
        (sp_yr[rownames(sp_yr) %in% myyear, 3]- sp_yr[rownames(sp_yr) %in% 2010, 3]) +
        model$Coefs[model$Variables %in% "ns(baseline_y, df = 4)4"]*
        (sp_yr[rownames(sp_yr) %in% myyear, 4]- sp_yr[rownames(sp_yr) %in% 2010, 4]))
  
}

rr_yr<- function(myyear,model){
  exp(model$Coefs[model$Variables %in% "ns(baseline_y, df = 4)1"]*
        (sp_yr[rownames(sp_yr) %in% myyear, 1]-sp_yr[rownames(sp_yr) %in% 2010, 1]) +
        model$Coefs[model$Variables %in% "ns(baseline_y, df = 4)2"]*
        (sp_yr[rownames(sp_yr) %in% myyear, 2]- sp_yr[rownames(sp_yr) %in% 2010, 2]) +
        model$Coefs[model$Variables %in% "ns(baseline_y, df = 4)3"]*
        (sp_yr[rownames(sp_yr) %in% myyear, 3]- sp_yr[rownames(sp_yr) %in% 2010, 3]) +
        model$Coefs[model$Variables %in% "ns(baseline_y, df = 4)4"]*
        (sp_yr[rownames(sp_yr) %in% myyear, 4]- sp_yr[rownames(sp_yr) %in% 2010, 4]))
  
}

rr_yr2 <- function(myyear,model){
  exp(model$Coefs[model$Variables %in% "ns(year, df = 4)1"]*
        (sp_yr[rownames(sp_yr) %in% myyear, 1]-sp_yr[rownames(sp_yr) %in% 2010, 1]) +
        model$Coefs[model$Variables %in% "ns(year, df = 4)2"]*
        (sp_yr[rownames(sp_yr) %in% myyear, 2]- sp_yr[rownames(sp_yr) %in% 2010, 2]) +
        model$Coefs[model$Variables %in% "ns(year, df = 4)3"]*
        (sp_yr[rownames(sp_yr) %in% myyear, 3]- sp_yr[rownames(sp_yr) %in% 2010, 3]) +
        model$Coefs[model$Variables %in% "ns(year, df = 4)4"]*
        (sp_yr[rownames(sp_yr) %in% myyear, 4]- sp_yr[rownames(sp_yr) %in% 2010, 4]))
  
}
# my function to calculate variance and SE
spline_se_yr <- function(myyear){
  
  tmp <- v_spline
  a <- sp_yr$X1[rownames(sp_yr) %in% myyear]-sp_yr$X1[rownames(sp_yr) %in% 2010]
  b <- sp_yr$X2[rownames(sp_yr) %in% myyear]-sp_yr$X2[rownames(sp_yr) %in% 2010]
  c <- sp_yr$X3[rownames(sp_yr) %in% myyear]-sp_yr$X3[rownames(sp_yr) %in% 2010]
  d <- sp_yr$X4[rownames(sp_yr) %in% myyear]-sp_yr$X4[rownames(sp_yr) %in% 2010]
  
  Var_x <- tmp["ns(baseline_y, df = 4)1", "ns(baseline_y, df = 4)1"]
  Var_y <- tmp["ns(baseline_y, df = 4)2", "ns(baseline_y, df = 4)2"]
  Var_z <- tmp["ns(baseline_y, df = 4)3", "ns(baseline_y, df = 4)3"]
  Var_w <- tmp["ns(baseline_y, df = 4)4", "ns(baseline_y, df = 4)4"]
  
  Cov_xy <- tmp["ns(baseline_y, df = 4)1", "ns(baseline_y, df = 4)2"]
  Cov_xz <- tmp["ns(baseline_y, df = 4)1", "ns(baseline_y, df = 4)3"]
  Cov_xw <- tmp["ns(baseline_y, df = 4)1", "ns(baseline_y, df = 4)4"]
  Cov_yz <- tmp["ns(baseline_y, df = 4)2", "ns(baseline_y, df = 4)3"]
  Cov_yw <- tmp["ns(baseline_y, df = 4)2", "ns(baseline_y, df = 4)4"]
  Cov_zw <- tmp["ns(baseline_y, df = 4)3", "ns(baseline_y, df = 4)4"]
  
  pooled_var_splines_yr <- a^2*Var_x + b^2*Var_y + c^2*Var_z + d^2*Var_w +
    2*a*b*Cov_xy + 2*a*c*Cov_xz + 2*a*d*Cov_xw + 2*b*c*Cov_yz + 2*b*d*Cov_yw + 2*c*d*Cov_zw
  pooled_se_splines_yr <- sqrt(pooled_var_splines_yr)
  
  # print(pooled_var_splines_yr)
  print(pooled_se_splines_yr)
  
}

# my function to calculate variance and SE
spline_se_yr2 <- function(myyear){
  
  tmp <- v_spline
  a <- sp_yr$X1[rownames(sp_yr) %in% myyear]-sp_yr$X1[rownames(sp_yr) %in% 2010]
  b <- sp_yr$X2[rownames(sp_yr) %in% myyear]-sp_yr$X2[rownames(sp_yr) %in% 2010]
  c <- sp_yr$X3[rownames(sp_yr) %in% myyear]-sp_yr$X3[rownames(sp_yr) %in% 2010]
  d <- sp_yr$X4[rownames(sp_yr) %in% myyear]-sp_yr$X4[rownames(sp_yr) %in% 2010]
  
  Var_x <- tmp["ns(year, df = 4)1", "ns(year, df = 4)1"]
  Var_y <- tmp["ns(year, df = 4)2", "ns(year, df = 4)2"]
  Var_z <- tmp["ns(year, df = 4)3", "ns(year, df = 4)3"]
  Var_w <- tmp["ns(year, df = 4)4", "ns(year, df = 4)4"]
  
  Cov_xy <- tmp["ns(year, df = 4)1", "ns(year, df = 4)2"]
  Cov_xz <- tmp["ns(year, df = 4)1", "ns(year, df = 4)3"]
  Cov_xw <- tmp["ns(year, df = 4)1", "ns(year, df = 4)4"]
  Cov_yz <- tmp["ns(year, df = 4)2", "ns(year, df = 4)3"]
  Cov_yw <- tmp["ns(year, df = 4)2", "ns(year, df = 4)4"]
  Cov_zw <- tmp["ns(year, df = 4)3", "ns(year, df = 4)4"]
  
  pooled_var_splines_yr <- a^2*Var_x + b^2*Var_y + c^2*Var_z + d^2*Var_w +
    2*a*b*Cov_xy + 2*a*c*Cov_xz + 2*a*d*Cov_xw + 2*b*c*Cov_yz + 2*b*d*Cov_yw + 2*c*d*Cov_zw
  pooled_se_splines_yr <- sqrt(pooled_var_splines_yr)
  
  # print(pooled_var_splines_yr)
  print(pooled_se_splines_yr)
  
}

abs_val_rna <- function(myyear,model){
  (model$Coefs[model$Variables %in% "ns(log_rna, df = 4)1"]*
     (sp_rna[rownames(sp_rna) %in% myyear, 1]-sp_rna[rownames(sp_rna) %in% 5, 1]) +
     model$Coefs[model$Variables %in% "ns(log_rna, df = 4)2"]*
     (sp_rna[rownames(sp_rna) %in% myyear, 2]- sp_rna[rownames(sp_rna) %in% 5, 2]) +
     model$Coefs[model$Variables %in% "ns(log_rna, df = 4)3"]*
     (sp_rna[rownames(sp_rna) %in% myyear, 3]- sp_rna[rownames(sp_rna) %in% 5, 3]) +
     model$Coefs[model$Variables %in% "ns(log_rna, df = 4)4"]*
     (sp_rna[rownames(sp_rna) %in% myyear, 4]- sp_rna[rownames(sp_rna) %in% 5, 4]))
  
}

rr_rna<- function(myyear,model){
  exp(model$Coefs[model$Variables %in% "ns(log_rna, df = 4)1"]*
        (sp_rna[rownames(sp_rna) %in% myyear, 1]-sp_rna[rownames(sp_rna) %in% 5, 1]) +
        model$Coefs[model$Variables %in% "ns(log_rna, df = 4)2"]*
        (sp_rna[rownames(sp_rna) %in% myyear, 2]- sp_rna[rownames(sp_rna) %in% 5, 2]) +
        model$Coefs[model$Variables %in% "ns(log_rna, df = 4)3"]*
        (sp_rna[rownames(sp_rna) %in% myyear, 3]- sp_rna[rownames(sp_rna) %in% 5, 3]) +
        model$Coefs[model$Variables %in% "ns(log_rna, df = 4)4"]*
        (sp_rna[rownames(sp_rna) %in% myyear, 4]- sp_rna[rownames(sp_rna) %in% 5, 4]))
  
}


# my function to calculate variance and SE
spline_se_rna <- function(myyear){
  
  tmp <- v_spline
  a <- sp_rna$X1[rownames(sp_rna) %in% myyear]-sp_rna$X1[rownames(sp_rna) %in% 5]
  b <- sp_rna$X2[rownames(sp_rna) %in% myyear]-sp_rna$X2[rownames(sp_rna) %in% 5]
  c <- sp_rna$X3[rownames(sp_rna) %in% myyear]-sp_rna$X3[rownames(sp_rna) %in% 5]
  d <- sp_rna$X4[rownames(sp_rna) %in% myyear]-sp_rna$X4[rownames(sp_rna) %in% 5]
  
  Var_x <- tmp["ns(log_rna, df = 4)1", "ns(log_rna, df = 4)1"]
  Var_y <- tmp["ns(log_rna, df = 4)2", "ns(log_rna, df = 4)2"]
  Var_z <- tmp["ns(log_rna, df = 4)3", "ns(log_rna, df = 4)3"]
  Var_w <- tmp["ns(log_rna, df = 4)4", "ns(log_rna, df = 4)4"]
  
  Cov_xy <- tmp["ns(log_rna, df = 4)1", "ns(log_rna, df = 4)2"]
  Cov_xz <- tmp["ns(log_rna, df = 4)1", "ns(log_rna, df = 4)3"]
  Cov_xw <- tmp["ns(log_rna, df = 4)1", "ns(log_rna, df = 4)4"]
  Cov_yz <- tmp["ns(log_rna, df = 4)2", "ns(log_rna, df = 4)3"]
  Cov_yw <- tmp["ns(log_rna, df = 4)2", "ns(log_rna, df = 4)4"]
  Cov_zw <- tmp["ns(log_rna, df = 4)3", "ns(log_rna, df = 4)4"]
  
  pooled_var_splines_yr <- a^2*Var_x + b^2*Var_y + c^2*Var_z + d^2*Var_w +
    2*a*b*Cov_xy + 2*a*c*Cov_xz + 2*a*d*Cov_xw + 2*b*c*Cov_yz + 2*b*d*Cov_yw + 2*c*d*Cov_zw
  pooled_se_splines_yr <- sqrt(pooled_var_splines_yr)
  
  # print(pooled_var_splines_yr)
  print(pooled_se_splines_yr)
  
}

uni.fit.gee1 <- function(model) {
  s <-  summary(model)
  rownam <- rownames(s$coefficients)
  s <-  s$coefficients[-1,]
  sand <- sqrt(diag(model$robust.variance))
  r <- as.data.frame(list(exp(s[1]), exp(s[1] - 1.96*sand[2]),exp(s[1] + 1.96*sand[2]),round(1 - pnorm(abs(s[1]/sand[2]))) * 2),digits =4)
  
  colnames(r) <- c( "RR","Lower CI", "Upper CI","p-value")
  r$RR <- round(r$RR,digits =2)
  r$`Lower CI` <- round(r$`Lower CI`,digits=3)
  r$`Upper CI` <- round(r$`Upper CI`,digits =3)
  rownames(r) <- rownam[2]
  return(r)
}


rr_sq_sh <- function(cdref,model){
  exp(model$Coefs[model$Variables %in% "ns(sq_cd4, df = 2)1"]*
        (sp_cd4[rownames(sp_cd4) %in% cdref, 1]-sp_cd4[rownames(sp_cd4) %in% 15, 1]) +
        model$Coefs[model$Variables %in% "ns(sq_cd4, df = 2)2"]*
        (sp_cd4[rownames(sp_cd4) %in% cdref, 2]- sp_cd4[rownames(sp_cd4) %in% 15, 2]))
  
}


# my function to calculate variance and SE
spline_se_sh <- function(cdref){
  
  tmp <- v_spline
  a <- sp_cd4$X1[rownames(sp_cd4) %in% cdref]-sp_cd4$X1[rownames(sp_cd4) %in% 14.14]
  b <- sp_cd4$X2[rownames(sp_cd4) %in% cdref]-sp_cd4$X2[rownames(sp_cd4) %in% 14.14]  
  Var_x <- tmp["ns(sq_cd4, df = 2)1", "ns(sq_cd4, df = 2)1"]
  Var_y <- tmp["ns(sq_cd4, df = 2)2", "ns(sq_cd4, df = 2)2"]
  
  Cov_xy <- tmp["ns(sq_cd4, df = 2)1", "ns(sq_cd4, df = 2)2"]
  
  
  pooled_var_splines_yr <- a^2*Var_x + b^2*Var_y + 2*a*b*Cov_xy 
  pooled_se_splines_yr <- sqrt(pooled_var_splines_yr)
  
  # print(pooled_var_splines_yr)
  print(pooled_se_splines_yr)
  
}


abs_sh <- function(myyear,model){
  (model$Coefs[model$Variables %in% "ns(sq_cd4, df = 2)1"]*
     (sp_cd4[rownames(sp_cd4) %in% myyear, 1]-sp_cd4[rownames(sp_cd4) %in% 15, 1]) +
     model$Coefs[model$Variables %in% "ns(sq_cd4, df = 2)2"]*
     (sp_cd4[rownames(sp_cd4) %in% myyear, 2]- sp_cd4[rownames(sp_cd4) %in% 15, 2]))
  
}

#rm(list=ls())


library(readr)
library(haven)
library(Hmisc)
library(dplyr)
library(tangram)
library(survival)
library(reshape2)
library(kableExtra)
library(survminer)
library(cmprsk)
library(rms)
library(psych)
library(gtsummary)
library(lattice)
library(survival)
library(purrr)
library(lubridate)
library(writexl)
library(countrycode)
library(splines)
library(sandwich)
library(lme4)
library(readxl)
library(gee)
library(broom)
library(mice)
library(lmtest)


setwd("C:/Users/ranadip/Desktop/Projects/Histo Project")

histo_demo <- read_csv("VCCC data/histo_demographics_16jun2023.csv",show_col_types = FALSE) 
histo_dx <- read_csv("VCCC data/histo_dx_17may2023.csv",show_col_types = FALSE)
histo_labs <- read_csv("VCCC data/histo_labs_21mar2023.csv",show_col_types = FALSE)

# Some IDS have repeated Histo Diagnosis dates so taking minimun of these ages
#histo_dx <- histo_dx %>% group_by(CFAR_PID) %>%  slice(which.min(AGE_AT_DX_ONSET))

ids <- read_csv("Corrected data/basic.csv",show_col_types = FALSE)
ids$peru <- with(ids,ifelse(site == "peru" & enrol_d >= "2019-12-20",0,1))
ids <- ids %>% filter(peru == 1)
basic <- read_csv("ccasanet_database_20230828/basic.csv",show_col_types = FALSE)
basic$peru <- with(basic,ifelse(site == "peru" & enrol_d >= "2019-12-20",0,1))
basic <- basic %>% filter(peru == 1)
ce <- read_csv("ccasanet_database_20230828/ce.csv",show_col_types = FALSE)
ce$peru <- with(ce,ifelse(site == "peru" & ce_d >= "2019-12-20",0,1))
ce <- ce %>% filter(peru == 1)
center <- read_csv("ccasanet_database_20230828/center.csv",show_col_types = FALSE)
ce_tb <- read_csv("ccasanet_database_20230828/ce_tb.csv",show_col_types = FALSE)
ce_cancer <- read_csv("ccasanet_database_20230828/ce_cancer.csv",show_col_types = FALSE)
origin <- read_csv("ccasanet_database_20230828/origin.csv",show_col_types = FALSE)
art <- read_csv("ccasanet_database_20230828/art.csv",show_col_types = FALSE)
art$peru <- with(art,ifelse(site == "peru" & art_sd >= "2019-12-20",0,1))
art <- art %>% filter(peru == 1)
follow <- read_csv("ccasanet_database_20230828/follow.csv",show_col_types = FALSE)
follow$l_alive_d <- with(follow,if_else(site == "peru" & l_alive_d >= "2019-12-20",as.Date("2019-12-20"),as.Date(l_alive_d)))
visit <- read_csv("ccasanet_database_20230828/visit.csv",show_col_types = FALSE)
visit$peru <- with(visit,ifelse(site == "peru" & visit_d >= "2019-12-20",0,1))
visit <- visit %>% filter(peru == 1)
lab <- read_csv("ccasanet_database_20230828/lab.csv",show_col_types = FALSE)
lab$peru <- with(lab,ifelse(site == "peru" & lab_d >= "2019-12-20",0,1))
lab <- lab %>% filter(peru == 1)
lab_cd4 <- read_csv("ccasanet_database_20230828/lab_cd4.csv",show_col_types = FALSE)
lab_cd4$peru <- with(lab_cd4,ifelse(site == "peru" & cd4_d >= "2019-12-20",0,1))
lab_cd4 <- lab_cd4 %>% filter(peru == 1)
lab_rna <- read_csv("ccasanet_database_20230828/lab_rna.csv",show_col_types = FALSE)
lab_rna$peru <- with(lab_rna,ifelse(site == "peru" & rna_d >= "2019-12-20",0,1))
lab_rna <- lab_rna %>% filter(peru == 1)
#visit <- read_csv("ccasanet_database_20230828/origin.csv",show_col_types = FALSE)

histo_la1 <- read_excel("ccasanet_database_20230828/countries.xlsx", sheet = "CCASAnet")
histo_tn1 <- read_excel("ccasanet_database_20230828/countries.xlsx", sheet = "VCCC")

ids <- ids %>% select(patient_id)
ce <- left_join(ids,ce, by ="patient_id")
baisc <- left_join(ids,basic, by ="patient_id")
#center <- left_join(ids,center, by ="patient_id")
lab <- left_join(ids,lab, by ="patient_id")
lab_rna <- left_join(ids,lab_rna, by ="patient_id")
lab_cd4 <- left_join(ids,lab_cd4, by ="patient_id")
follow <- left_join(ids,follow, by ="patient_id")
visit <- left_join(ids,visit, by ="patient_id")
art <- left_join(ids,art, by ="patient_id")
ce_cancer <- left_join(ids,ce_cancer, by ="patient_id")
ce_tb <- left_join(ids,ce_tb, by ="patient_id")
origin <- left_join(ids,origin, by ="patient_id")




#Since Haiti has no histo cases removing Haiti from Basic table 60265
la <- basic %>% filter(site == "brazil" | site == "chile" | site == "honduras" |site == "mexico" | site == "peru")  #32939
ce <- ce %>%  filter(site == "brazil" | site == "chile" | site == "honduras" |site == "mexico" | site == "peru")
visit <- visit %>% filter(site == "brazil" | site == "chile" | site == "honduras" |site == "mexico" | site == "peru")
la <- la %>% filter(site == "brazil" | site == "chile" | site == "honduras" |site == "mexico" | site == "peru")

#Filtering enrollment date between 1 st January 2000 to 31st December 2021 and age >= 18 #28554 TN #6870
la <- la %>%  filter(between(enrol_d, as.Date("2000-01-01"), as.Date("2021-12-31"))) %>%  mutate(age_b = difftime(enrol_d,birth_d, units = "days")/ 365.25) %>% filter(age_b >=18) 
tn <- histo_demo %>%   filter(age_at_first_visit>=18) %>% filter(between(year_of_enrollment,2000,2021))

#First occurrence of histo
histo_la <- ce %>% filter(grepl("histoplasmosis", ce_id,ignore.case = TRUE))
histo_la <- histo_la %>% group_by(patient_id) %>% slice_min(ce_d)
histo_la$histo_d <- as.Date(histo_la$ce_d,origin="1970-01-01")
histo_la$histo <- "Yes"
histo_la <- histo_la %>% select(-site)
histo_tn <- histo_dx %>% filter(grepl("Histoplasmosis",diagnosis,ignore.case = TRUE))
histo_tn <- histo_tn %>% group_by(CFAR_PID) %>% slice_min(AGE_AT_DX_ONSET)
histo_tn$histo_age <- histo_tn$AGE_AT_DX_ONSET
histo_tn$histo <- "Yes"


#Diagnosis of Histo before Enrollment (DX_BEFORE_ENROLL) to be excluded(- 30 days till date of clinic entry) LA #28496, TN #6837
tn <- left_join(tn,histo_tn,by="CFAR_PID") 
tn <- tn %>% group_by(CFAR_PID) %>% mutate(diff = sum(-age_at_first_visit,histo_age,na.rm=FALSE)) %>% filter((diff*365.25) >= -30 | is.na(diff))
la <- left_join(la,histo_la, by = "patient_id")
la <- la %>% dplyr::select(birth_d,enrol_d,patient_id,ce_d,mode,male_y,histo_d,aids_first_d,histo,histo_d,age_b) %>% mutate(diff = (difftime(histo_d, enrol_d, units = "days"))) %>% filter(diff >= -30 | is.na(diff)) 

#Histo Count(Yes/No)

la$histo <- with(la,ifelse(is.na(histo),"No","Yes"))
tn$histo <- with(tn,ifelse(is.na(histo),"No","Yes"))
la$histo_num <- with(la,ifelse(histo == "No",0,1))
tn$histo_num <- with(tn,ifelse(histo == "No",0,1))



#Applying the VCCC condition on the CCASAnet data
#(2 visits within a year of 1st visit) using labs and art table 23571

follow <- follow %>% dplyr::select(-center)
la <- left_join(la,follow,by = "patient_id") %>% filter(l_alive_d != enrol_d)
v <- left_join(la,lab_rna, by = "patient_id") %>% dplyr::select(patient_id,rna_d,l_alive_d,enrol_d)
v <- left_join(v,lab_cd4, by = "patient_id") %>% dplyr::select(patient_id,cd4_d,rna_d,l_alive_d,enrol_d)
v <- left_join(v,art,by = "patient_id") %>% dplyr::select(patient_id,cd4_d,rna_d,l_alive_d,enrol_d,art_sd)
library(data.table)
setDT(v)
v <- v[ , .(date = c(cd4_d, rna_d,art_sd)), by = .(patient_id)]
v <- v %>% group_by(patient_id) %>% distinct(date) %>% filter(between(date, as.Date("2000-01-01"), as.Date("2021-12-31")))
v <- v %>% group_by(patient_id) %>%
  arrange(date)
v <- v %>% group_by(patient_id) %>% mutate(diff = (date - first(date)))
v <- v %>% filter(diff <= 365.25)
v <- v %>% group_by(patient_id) %>%  add_count() %>% mutate(v = ifelse(n>=2,1,0))
v <- v %>% group_by(patient_id) %>%  slice_min(date)
v <- v %>% dplyr::select(patient_id,date,v)
la <- left_join(la,v,by="patient_id")  
la <- la %>% filter(v == 1)

## Excliding 56 ppl with male_y ==8
la <- la %>% filter(male_y != 8)
tn <- tn %>% filter(birthsex != "Intersexed")


# visit <- visit %>% dplyr::select(patient_id,visit_d)
# visit <- visit %>% group_by(patient_id) %>% mutate(diff =  visit_d - first(visit_d))
# visit <- visit %>% filter(diff <= 365.25)
# visit <- visit %>%  add_count() %>% mutate(visit = ifelse(n>=2,1,0))
# visit <- visit %>% group_by(patient_id) %>%  slice_min(visit_d)
# la <- left_join(la,visit,by="patient_id")  
# la <- la %>% filter(visit == 1)


#Renaming variables in TN tables
tn$site <- "Tennessee"
tn$baseline_y <- tn$year_of_enrollment


#Creating baseline variables CD4, rna

la$baseline <- la$enrol_d
tn$baseline_y <- tn$year_of_enrollment
tn$age_b <- tn$age_at_first_visit
la$baseline_y <- format(as.Date(la$baseline, format="%d/%m/%Y"),"%Y")
cd4_b <- left_join(la, lab_cd4, by = "patient_id")
cd4_b$upper_limit <- as.Date(cd4_b$baseline) + 30
cd4_b$lower_limit <- as.Date(cd4_b$baseline) %m-% months(6)
cd4_b <- cd4_b %>% filter(cd4_d <= upper_limit & cd4_d >= lower_limit)
cd4_b <- cd4_b %>% mutate(difference = as.numeric(difftime(as.Date(baseline),as.Date(cd4_d),units ="days")))
cd4_b <- cd4_b %>% group_by(patient_id) %>% slice_min(abs(difference))
cd4_b <- cd4_b %>% group_by(patient_id) %>% slice_min(abs(cd4_v)) 
cd4_b <- cd4_b %>% group_by(patient_id) %>% slice_min(difference)
cd4_b$cd4_b <- cd4_b$cd4_v 
cd4_b$cd4_base_date <- cd4_b$cd4_d
cd4_b <- cd4_b %>% dplyr::select(patient_id,cd4_b,cd4_base_date)

rna_b <- left_join(la, lab_rna, by = "patient_id")
rna_b$upper_limit <- as.Date(rna_b$baseline) + 30
rna_b$lower_limit <- as.Date(rna_b$baseline) %m-% months(6)
rna_b <- rna_b %>% filter(rna_d <= upper_limit & rna_d >= lower_limit)
rna_b <- rna_b %>% mutate(difference = as.numeric(difftime(as.Date(baseline),as.Date(rna_d),units ="days")))
rna_b <- rna_b %>% group_by(patient_id) %>% slice_min(abs(difference))
rna_b <- rna_b %>% group_by(patient_id) %>% slice_min(abs(rna_v))
rna_b$rna_b <- rna_b$rna_v 
rna_b$rna_base_date <- rna_b$rna_d
rna_b <- rna_b %>% dplyr::select(patient_id,rna_b,rna_base_date)

#Creating baseline and diagnosis date variables for CD4, rna for TN

tn_cd4_b <- left_join(tn,histo_labs, by = "CFAR_PID") %>% filter(testName == "CD4 COUNT")
tn_cd4_b$upper_limit <- tn_cd4_b$age_at_first_visit + 0.0821
tn_cd4_b$lower_limit <- tn_cd4_b$age_at_first_visit - 0.4928
tn_cd4_b <- tn_cd4_b %>% filter(AGE_AT_RESULT_DATE <= upper_limit & AGE_AT_RESULT_DATE >= lower_limit)
tn_cd4_b <- tn_cd4_b %>% mutate(difference = age_at_first_visit - AGE_AT_RESULT_DATE)
tn_cd4_b <- tn_cd4_b %>% group_by(CFAR_PID) %>% slice_min(abs(difference))
tn_cd4_b$cd4_b <- tn_cd4_b$RESULT_NUMERIC 
tn_cd4_b <- tn_cd4_b %>% group_by(CFAR_PID) %>% slice_min(cd4_b)
tn_cd4_b <- tn_cd4_b %>% dplyr::select(CFAR_PID,cd4_b)

tn_rna_b <- left_join(tn, histo_labs, by = "CFAR_PID") %>% filter(testName == "HIV-1 RNA")
tn_rna_b$upper_limit <- tn_rna_b$age_at_first_visit + 0.0821
tn_rna_b$lower_limit <- tn_rna_b$age_at_first_visit - 0.4928
tn_rna_b <- tn_rna_b %>% filter(AGE_AT_RESULT_DATE <= upper_limit & AGE_AT_RESULT_DATE >= lower_limit)
tn_rna_b <- tn_rna_b %>% mutate(difference = age_at_first_visit - AGE_AT_RESULT_DATE)
tn_rna_b <- tn_rna_b %>% group_by(CFAR_PID) %>% slice_min(difference)
tn_rna_b$rna_b <- tn_rna_b$RESULT_NUMERIC 
tn_rna_b <- tn_rna_b %>% group_by(CFAR_PID) %>% slice_max(rna_b)
tn_rna_b <- tn_rna_b %>% dplyr::select(CFAR_PID,rna_b)

tn_cd4_d <- left_join(tn,histo_labs, by = "CFAR_PID") %>% filter(testName == "CD4 COUNT")
tn_cd4_d$upper_limit <- tn_cd4_d$AGE_AT_DX_ONSET + 0.0821
tn_cd4_d$lower_limit <- tn_cd4_d$AGE_AT_DX_ONSET - 0.4928
tn_cd4_d <- tn_cd4_d %>% filter(AGE_AT_RESULT_DATE <= upper_limit & AGE_AT_RESULT_DATE >= lower_limit)
tn_cd4_d <- tn_cd4_d %>% mutate(difference = AGE_AT_DX_ONSET - AGE_AT_RESULT_DATE)
tn_cd4_d <- tn_cd4_d %>% group_by(CFAR_PID) %>% slice_min(difference)
tn_cd4_d$cd4_diag <- tn_cd4_d$RESULT_NUMERIC 
tn_cd4_d <- tn_cd4_d %>% group_by(CFAR_PID) %>% slice_min(cd4_diag)
tn_cd4_d <- tn_cd4_d %>% dplyr::select(CFAR_PID,cd4_diag)

tn_rna_d <- left_join(tn, histo_labs, by = "CFAR_PID") %>% filter(testName == "HIV-1 RNA")
tn_rna_d$upper_limit <- tn_rna_d$AGE_AT_DX_ONSET + 0.0821
tn_rna_d$lower_limit <- tn_rna_d$AGE_AT_DX_ONSET - 0.4928
tn_rna_d <- tn_rna_d %>% filter(AGE_AT_RESULT_DATE <= upper_limit & AGE_AT_RESULT_DATE >= lower_limit)
tn_rna_d <- tn_rna_d %>% mutate(difference = AGE_AT_DX_ONSET - AGE_AT_RESULT_DATE)
tn_rna_d <- tn_rna_d %>% group_by(CFAR_PID) %>% slice_min(difference)
tn_rna_d$rna_diag <- tn_rna_d$RESULT_NUMERIC 
tn_rna_d <- tn_rna_d %>% group_by(CFAR_PID) %>% slice_max(rna_diag)
tn_rna_d <- tn_rna_d %>% dplyr::select(CFAR_PID,rna_diag)


#Risk
la$mode <- as.factor(la$mode)
la$risk <- as.numeric(la$mode)
la$risk <- ifelse(la$risk %in% c(1, 2,3,5), 0,
                  ifelse(la$risk %in% c(4), 1,
                         ifelse(la$risk %in% c(6), 2,
                                ifelse(la$risk %in% c(7,8,9,10), 3, NA))))
la$risk1 <- with(la,ifelse(is.na(risk),3,risk))
la$risk1 <- factor(la$risk1, levels=0:3, labels=c("MSM", "Hetero", "IDU", "Other/Unknown"))

tn$r <- as.factor(tn$risk1)
tn$risk <- as.numeric(tn$r)
tn$risk <- ifelse(tn$risk %in% c(4), 0,
                  ifelse(tn$risk %in% c(2), 1,
                         ifelse(tn$risk %in% c(3), 2,
                                ifelse(tn$risk %in% c(1,5,6,7,8,9,10), 3, NA))))
tn$risk1 <- with(tn,ifelse(is.na(risk),3,risk))
tn$risk1 <- factor(tn$risk1, levels=0:3, labels=c("MSM", "Hetero", "IDU", "Other/Unknown"))

#getting CD$,rna at diagnosis

cd4_d <- left_join(la, lab_cd4, by = "patient_id")
cd4_d$upper_limit <- as.Date(cd4_d$ce_d) + 30
cd4_d$lower_limit <- as.Date(cd4_d$ce_d) %m-% months(6)
cd4_d <- cd4_d %>% filter(cd4_d <= upper_limit & cd4_d >= lower_limit)
cd4_d <- cd4_d %>% mutate(difference = as.numeric(difftime(as.Date(ce_d),as.Date(cd4_d),units ="days")))
cd4_d <- cd4_d %>% group_by(patient_id) %>% slice_min(difference)
cd4_d$cd4_diag <- cd4_d$cd4_v 
cd4_d <- cd4_d %>% dplyr::select(patient_id,cd4_diag)

rna_d <- left_join(la, lab_rna, by = "patient_id")
rna_d$upper_limit <- as.Date(rna_d$ce_d) + 30
rna_d$lower_limit <- as.Date(rna_d$ce_d) %m-% months(6)
rna_d <- rna_d %>% filter(rna_d <= upper_limit & rna_d >= lower_limit)
rna_d <- rna_d %>% mutate(difference = as.numeric(difftime(as.Date(ce_d),as.Date(rna_d),units ="days")))
rna_d <- rna_d %>% group_by(patient_id) %>% slice_min(difference)
rna_d$rna_diag <- rna_d$rna_v 
rna_d <- rna_d %>% dplyr::select(patient_id,rna_diag)



#ART Naive(Yes/No) at enrollment date or baseline
art <- art %>% dplyr::select(-site,-center)
art_1 <- art %>% group_by(patient_id) %>% slice_min(art_sd)
la <- left_join(la,art_1,by = "patient_id")
la$naive <- with(la,ifelse(art_sd >= enrol_d,"Yes","No"))
#art_1$naive <- factor(art_1$naive, levels=0:1, labels=c("No","Yes"))
#art_1 <- art_1 %>% dplyr::select(naive,art_sd,patient_id,ii1,ii2,nnrti1,pi,nrti,nnrti2,t20,ccr5)
tn$naive <- with(tn,ifelse(AGE_AT_ART_START >= age_at_first_visit,"Yes","No"))
#tn$naive <- factor(tn$naive, levels=0:1, labels=c("No","Yes"))

#Gender at Birth
tn$sex <- tn$birthsex
la$sex <- with(la,ifelse(male_y == 1, "Male",
                         ifelse(male_y == 0, "Female", "Transgender")))

#ART regiment at baseline
art_1$init_reg_class <- with(art_1,  ifelse(ii1 > 0 | ii2 > 0, "IINSTI-based",
                                            ifelse(nnrti1 > 0 | nnrti2 > 0 & pi==0, "NNRTI-based",
                                                   ifelse(pi > 0, "PI-based", "Other"))))

tn$init_reg_class <- with(tn,  ifelse(II > 0, "IINSTI-based",
                                      ifelse(NNRTI > 0 & PI ==0, "NNRTI-based",
                                             ifelse(PI > 0, "PI-based", "Other"))))


#AIDS Defining illness
la <- la %>% dplyr::select(-ce_d)
ade <- left_join(la,ce, by = "patient_id") 
ade$upper_limit <- as.Date(ade$baseline) + 30
ade$lower_limit <- as.Date(ade$baseline) %m-% months(6)
ade <- ade %>% filter(ce_d  <= upper_limit & ce_d>= lower_limit)
ade <- ade %>% mutate(difference = as.numeric(difftime(as.Date(baseline),as.Date(ce_d),units ="days")))
ade <- ade %>% mutate(ade_type= if_else(grepl('^ade', ce_id) & ce_id != "ade_histoplasmosis" , 'Yes', 'No')) 
ade <- ade %>% filter(ade_type == "Yes") 
ade <- ade %>% group_by(patient_id) %>% distinct(patient_id,.keep_all = TRUE) %>% dplyr::select(patient_id,ade_type)

#For TN
tn$age_at_histo <- tn$AGE_AT_DX_ONSET
tn <- tn %>% dplyr::select(-AGE_AT_DX_ONSET,-diagnosis)
ade_tn <- left_join(tn,histo_dx, by = "CFAR_PID") 
ade_tn$upper_limit <- ade_tn$age_at_first_visit + 0.0821
ade_tn$lower_limit <- ade_tn$age_at_first_visit - 0.4928
ade_tn <- ade_tn %>% filter(AGE_AT_DX_ONSET  <= upper_limit & AGE_AT_DX_ONSET>= lower_limit)
ade_tn <- ade_tn %>% mutate(difference = as.numeric(age_at_first_visit- AGE_AT_DX_ONSET))
ade_tn <- ade_tn %>% mutate(ade_type= if_else( diagnosis != "Histoplasmosis disseminated or extrapulmonary" , 'Yes', 'No'))
ade_tn <- ade_tn %>% filter(ade_type == "Yes") 
ade_tn <- ade_tn %>% group_by(CFAR_PID) %>% distinct(CFAR_PID,.keep_all = TRUE) %>% dplyr::select(CFAR_PID,ade_type,AGE_AT_DX_ONSET)

#TB Diagnosis
#AIDS Defining illness

tb <- left_join(la,ce_tb, by = "patient_id") 
tb$upper_limit <- as.Date(tb$baseline) + 30
tb$lower_limit <- as.Date(tb$baseline) %m-% months(6)
tb <- tb %>% filter(tbdiagnosis_d<= upper_limit & tbdiagnosis_d>= lower_limit)
tb <- tb %>% mutate(difference = as.numeric(difftime(as.Date(baseline),as.Date(tbdiagnosis_d),units ="days")))
tb <- tb %>% group_by(patient_id) %>% slice_min(tbdiagnosis_d)
tb$tb <- "Yes" 
tb <- tb %>% group_by(patient_id) %>% distinct(patient_id,.keep_all = TRUE) %>% dplyr::select(patient_id,tb,cxpos_y,tbdiagnosis_d,afbpos_y,pcrpos_y,genexpertpos_y,modspos_y,otherpos_y,tbtx_init)


tb$tb_test <- with(tb,ifelse(cxpos_y == 1 | afbpos_y == 1 | pcrpos_y == 1 | genexpertpos_y == 1 | modspos_y == 1 |
                               otherpos_y == 1, "Yes", "No"))
tb$cxpos_y <- with(tb,ifelse(cxpos_y == 1, "Yes",
                             ifelse(cxpos_y == 2, "No",
                                    ifelse(cxpos_y == 8,  "culture not performed",
                                           ifelse(cxpos_y ==9, "Unknown", NA)))),levels =c (1:2,8,9))
tb$afbpos_y <- with(tb,ifelse(afbpos_y == 1, "Yes",
                              ifelse(afbpos_y == 2, "No",
                                     ifelse(afbpos_y == 8,  "culture not performed",
                                            ifelse(afbpos_y ==9, "Unknown", NA)))),levels =c (1:2,8,9))
tb$pcrpos_y <- with(tb,ifelse(pcrpos_y == 1, "Yes",
                              ifelse(pcrpos_y == 2, "No",
                                     ifelse(pcrpos_y == 8,  "culture not performed",
                                            ifelse(pcrpos_y ==9, "Unknown", NA)))),levels =c (1:2,8,9))
tb$genexpertpos_y <- with(tb,ifelse(genexpertpos_y == 1, "Yes",
                                    ifelse(genexpertpos_y == 2, "No",
                                           ifelse(genexpertpos_y == 8,  "culture not performed",
                                                  ifelse(genexpertpos_y ==9, "Unknown", NA)))),levels =c (1:2,8,9))
tb$modspos_y <- with(tb,ifelse(modspos_y == 1, "Yes",
                               ifelse(modspos_y == 2, "No",
                                      ifelse(modspos_y == 8,  "culture not performed",
                                             ifelse(modspos_y ==9, "Unknown", NA)))),levels =c (1:2,8,9))
tb$otherpos_y <- with(tb,ifelse(otherpos_y == 1, "Yes",
                                ifelse(otherpos_y == 2, "No",
                                       ifelse(otherpos_y == 8,  "culture not performed",
                                              ifelse(otherpos_y ==9, "Unknown", NA)))),levels =c (1:2,8,9))
tb$tb_trt <- with(tb,ifelse(is.na(tbtx_init),"No","Yes"))



label(tb$cxpos_y)  <- "Culture positive tuberculosis"
label(tb$afbpos_y) <- "AFB Smear-positive"
label(tb$pcrpos_y) <- "PCR positive"
label(tb$genexpertpos_y) <- "GeneXpert positive"
label(tb$modspos_y) <- "MODS positive"
label(tb$otherpos_y) <- "Other diagnostic test results positive"
label(tb$tb_trt) <- "TB treatment"
label(tb$tb_test) <- "TB test yes/no"

#Creating baseline dates and other variables
#follow <- follow %>% dplyr::select(-center)
la <- list(la,cd4_b,cd4_d,rna_b,rna_d,ade,tb) %>% reduce(left_join, by = "patient_id")

la$last_d <- with(la,ifelse(l_alive_d < as.Date("2021-12-31"),l_alive_d,as.Date("2021-12-31")))
la$last_d <- as.Date(la$last_d,origin="1970-01-01")
la$baseline <- as.Date(la$baseline,origin="1970-01-01")
la$follow <- with(la,ifelse(histo == "Yes",abs(difftime(as.Date(histo_d),as.Date(baseline),units ="days")/365.25),
                            abs(difftime(as.Date(last_d),as.Date(baseline),units ="days")/365.25)))




#la$age_at_onset <- with(la,difftime(as.Date(baseline),a.Date(ce_d),units ="days"))/365.25
tn <- list(tn,tn_cd4_b,tn_cd4_d,tn_rna_b,tn_rna_d,ade_tn) %>% reduce(left_join, by = "CFAR_PID")
tn <- tn %>% distinct(CFAR_PID, .keep_all = TRUE)
la$SITE <- with(la,ifelse(site == "argentina", "Argentina",
                          ifelse(site == "brazil", "Brazil",
                                 ifelse(site == "chile", "Chile",
                                        ifelse(site == "hondur", "Honduras", 
                                               ifelse(site == "mexico", "Mexico", 
                                                      ifelse(site == "peru", "Peru","Haiti")))))))

la$tb <- with(la,ifelse(is.na(tb),"No",tb))
la$cp_d <- with(la,ifelse(cxpos_y == 2,tbdiagnosis_d,NA))


la$aids_first_y <- format(as.Date(la$aids_first_d, format="%d/%m/%Y"),"%Y")
la$l_alive_y <- format(as.Date(la$l_alive_d, format="%d/%m/%Y"),"%Y")
la$histo_y <- format(as.Date(la$histo_d, format="%d/%m/%Y"),"%Y")
la$death_y <- format(as.Date(la$death_d, format="%d/%m/%Y"),"%Y")
la$aids_first_y <- as.numeric(la$aids_first_y)
la$l_alive_y <- as.numeric(la$l_alive_y)
la$histo_y <- as.numeric(la$histo_y)
la$death_y <- as.numeric(la$death_y)
la$baseline_y <- as.numeric(la$baseline_y)
la$ade_type <- with(la,ifelse(!is.na(ade_type), "Yes", "No"))


#Categories for CD4, RNA
la$rna_bb <- ifelse(la$rna_b < 400, 0, 1)
la$rna_bb <- factor(la$rna_bb, levels = 0:1, labels=c("Undetectable", "Detectable"))


tn$rna_bb <- ifelse(tn$rna_b < 400, 0, 1)
tn$rna_bb <- factor(tn$rna_bb, levels = 0:1, labels=c("Undetectable", "Detectable"))



la$rna_db <- ifelse(la$rna_diag < 400, 0, 1)
la$rna_db <- factor(la$rna_db, levels = 0:1, labels=c("Undetectable", "Detectable"))



tn$rna_db <- ifelse(tn$rna_diag < 400, 0, 1)
tn$rna_db <- factor(tn$rna_db, levels = 0:1, labels=c("Undetectable", "Detectable"))


la$cd4_bb <- ifelse(la$cd4_b <= 350, 0, 1)
la$cd4_bb <- factor(la$cd4_bb, levels = 0:1, labels=c("<=350 cells/ mm3", ">350 cells/ mm3"))


tn$cd4_bb <- ifelse(tn$cd4_b <= 350, 0, 1)
tn$cd4_bb <- factor(tn$cd4_bb, levels = 0:1, labels=c("<=350 cells/ mm3", ">350 cells/ mm3"))



la$cd4_db <- ifelse(la$cd4_diag <= 350, 0, 1)
la$cd4_db <- factor(la$cd4_db, levels = 0:1, labels=c("<=350 cells/ mm3", ">350 cells/ mm3"))



tn$cd4_db <- ifelse(tn$cd4_diag <= 350, 0, 1)
tn$cd4_db <- factor(tn$cd4_db, levels = 0:1, labels=c("<=350 cells/ mm3", ">350 cells/ mm3"))



#Culture Positive 
#la$cp <- with(la,ifelse(cxpos_y == 0,"No", "Yes"))

#Time from Baseline to histo diagnosis
la$histo_t <- with(la,difftime(as.Date(histo_d),as.Date(baseline),units ="days")/365.25)
tn$histo_t <- with(tn,age_at_histo - age_at_first_visit)

#follow up time in TN
tn$follow_up <- with(tn,ifelse(histo == "Yes", (age_at_histo - age_at_first_visit),
                               ifelse(CENSOR_TYPE == "Death" | CENSOR_TYPE == "12/31/2021", (min(age_at_last_visit,age_at_censor) - age_at_first_visit),(age_at_last_visit-age_at_first_visit))))

#Time till ART
la$art_t <- with(la,difftime(as.Date(art_sd),as.Date(baseline),units ="days"))
la$art_histo <- with(la,difftime(as.Date(art_sd),as.Date(histo_d),units ="days"))
la$pre_post_histo <- with(la,ifelse(histo =="Yes" & (art_histo <= 0 | is.na(art_sd)) ,"Histo diagnosis after ART start",
                                    ifelse(histo == "Yes" & art_histo >0,"Histo diagnosis before/at ART start",NA)))
la$art_pre <- with(la,ifelse(pre_post_histo == "Histo diagnosis after ART start",abs(art_histo),NA))
la$art_post <- with(la,ifelse(pre_post_histo == "Histo diagnosis before/at ART start",abs(art_histo),NA))
###TIME updated table
# rna <- lab_rna %>% dplyr::select(patient_id,rna_v,rna_d)
# cd <- lab_cd4 %>% dplyr::select(patient_id,cd4_d,cd4_v)
# ll <- list(la,rna,cd) %>% reduce(left_join, by = "patient_id")
#  ll <- ll %>% dplyr::select(patient_id,cd4_v,cd4_d,rna_v,rna_d,baseline,l_alive_d,)
# library(data.table)
# setDT(ll)
# ll_n <- ll[ , .(Date = c(cd4_d, rna_d)), by = .(patient_id, cd4_d, rna_d,baseline,l_alive_d,histo_d)]
# l <- ll %>% dplyr::select(patient_id,art_sd,baseline,cd4_v,rna_v,l_alive_d,SITE,histo_d)
# t <- ll_n %>% group_by(patient_id) %>% distinct(Date,  .keep_all =  TRUE)
#l <- left_join(t,by="patient_id")

#t <- t %>% group_by(patient_id) %>% mutate(time.diff = difftime(Date, lag(Date), unit = 'days'))

la$cp_d <- as.Date(la$cp_d, origin="1970-01-01")
la$cp_y <- as.numeric(substr(la$cp_d, 1, 4))
la$rna_b <- abs(la$rna_b)
la$rna_diag <-  abs(la$rna_diag)

# Time till Histo
tn$time_h <- with(tn,(age_at_histo - age_at_first_visit)*365.25)
la$time_h <- with(la,difftime(histo_d,enrol_d, units = "days"))
tn$art_histo <- with(tn,(AGE_AT_ART_START-age_at_histo)*365.25)
tn$pre_post_histo <- with(tn,ifelse(histo =="Yes" & (art_histo <= 0 | is.na(AGE_AT_ART_START)) ,"Histo diagnosis after ART start",
                                    ifelse(histo == "Yes" & art_histo >0,"Histo diagnosis before/at ART start",NA)))
tn$art_pre <- with(tn,ifelse(pre_post_histo == "Histo diagnosis after ART start",abs(art_histo),NA))
tn$art_post <- with(tn,ifelse(pre_post_histo == "Histo diagnosis before/at ART start",abs(art_histo),NA))

# Country of origin
origin <- origin %>% dplyr::select(patient_id,origin)
la <- left_join(la,origin,by = "patient_id")

la$Country <- with(la, ifelse(origin == "OUTSIDE", NA, origin))
la$Country <- countrycode(la$Country, origin = "iso3c", destination = "country.name")
la$country_cat <- with(la,ifelse(!is.na(la$Country) & tolower(Country) != tolower(site), "Yes",
                                 ifelse(!is.na(la$Country) & tolower(Country) == tolower(site), "No",NA)))

tn$country_cat <- with(tn,ifelse(birthcountry == "UNITED STATES", "No", 
                                 ifelse(!is.na(birthcountry) & birthcountry != "UNITED STATES","Yes", NA)))

## creating excel for unique countries in r 
ccasanet_countries <- as.data.frame(unique(la$Country))
colnames(ccasanet_countries) <- c("Countries")
vccc_countries <- as.data.frame(unique(tn$birthcountry))
colnames(vccc_countries) <- c("Countries")
#write_xlsx(ccasanet_countries, 'C:\\Users\\ranadip\\Desktop\\Projects\\Histo Project\\countries.xlsx')
#write_xlsx(vccc_countries,'C:\\Users\\ranadip\\Desktop\\Projects\\Histo Project\\countries2.xlsx')

## NA for outliers in CD$ at baseline
la$cd4_b <- ifelse(la$cd4_b > 4000, NA,la$cd4_b)


la$cd4_b_sq <- sqrt(la$cd4_b)
la$cd4_diag_sq <- sqrt(la$cd4_diag)
la$rna_b_log <- log10(abs(la$rna_b))
la$rna_diag_log <- log10(abs(la$rna_diag))

tn$cd4_b_sq <- sqrt(tn$cd4_b)
tn$cd4_diag_sq <- sqrt(tn$cd4_diag)
tn$rna_b_log <- log10(abs(tn$rna_b))
tn$rna_diag_log <- log10(abs(tn$rna_diag))

## Histo-endemic countries
histo_la_yes <- histo_la1 %>% filter(`Endemic(Yes/No)` == 1)
histo_la_no <- histo_la1 %>% filter(`Endemic(Yes/No)` == 0)
histo_tn_yes <- histo_tn1 %>% filter(`Endemic(Yes/No)` == 1)
histo_tn_no <- histo_tn1 %>% filter(`Endemic(Yes/No)` == 0)
la$ende_yes <- ifelse(la$Country %in% histo_la_yes$Countries,1,NA)
la$ende_no <- ifelse(la$Country %in% histo_la_no$Countries,1,NA)
tn$ende_yes <- ifelse(tn$birthcountry %in% histo_tn_yes$Countries,1,NA)
tn$ende_no <- ifelse(tn$birthcountry %in% histo_tn_no$Countries,1,NA)

la$birth_ende <- with(la,ifelse(!is.na(ende_yes) & ende_yes == 1,"Yes", 
                                ifelse(is.na(ende_yes) & ende_no==1 ,"No",NA)))
tn$birth_ende <- with(tn,ifelse(!is.na(ende_yes) & ende_yes == 1,"Yes", 
                                ifelse(is.na(ende_yes) & ende_no==1 ,"No",NA)))



la$mig_ende <- with(la,ifelse(birth_ende == "Yes" & country_cat == "Yes","Yes","No"))
tn$mig_ende <- with(tn,ifelse(birth_ende == "Yes" & country_cat == "Yes","Yes","No"))

tn$ade_type <- with(tn,ifelse(ade_type == "Yes","Yes","No"))

tn$death <- with(tn,ifelse(!is.na(age_at_death), "Yes","No"))
la$death <- with(la,ifelse(!is.na(death_d), "Yes", "No"))
#Labels
label(tn$AGE_AT_ART_START) <- "Age at first ART regimen start"
label(tn$risk1) <- "HIV acquisition risk factor"
label(tn$baseline_y) <- "Baseline Year"
label(tn$race) <- "Race"
label(tn$cd4_b) <- "CD4 at baseline"
label(tn$cd4_diag) <- "CD4 at Histoplasmosis diagnosis"
label(tn$rna_b) <- "HIV-RNA  at Baseline"
label(tn$rna_b_log) <- "log HIV-RNA  at Baseline"
label(tn$rna_diag) <- "HIV-RNA at Histoplasmosis diagnosis"
label(tn$birthsex) <- "Sex at birth"
label(tn$age_at_first_visit) <- "Age at Baseline"
label(tn$age_at_histo) <- "Age at Histoplasmosis diagnosis"
label(tn$init_reg_class) <- "First ART regimen (PI-, NNRTI-, or INSTI-based or other)"
label(tn$ade_type) <- "TB (Yes/No)"
label(tn$naive) <- "ART-naive at cohort entry (yes/no)"
label(tn$AGE_AT_ART_START) <- "Age at ART initiation"
label(tn$AGE_AT_FIRST_INFECTION) <- "AGE at HIV diagnosis"
label(tn$age_at_death) <- "Age at death"
label(tn$age_at_last_visit) <- "Age at last visit"
label(tn$year_of_enrollment) <- "Year of baseline"
label(tn$rna_bb) <- "VL at baseline binary"
label(tn$rna_db) <- "VL at diagnosis binary"
label(tn$cd4_bb) <- "CD4 at baseline binary"
label(tn$cd4_db) <- "CD4 at diagnosis binary"
label(tn$time_h) <- "Time till Histo diagnosis (days)"
label(tn$birthcountry) <- "Country of birth"
label(tn$country_cat) <- "Migration status (Yes/No)"
label(tn$presentsex) <- "Sex at present"
label(tn$rna_bb) <- "VL at baseline binary"
label(tn$rna_db) <- "VL at diagnosis binary"
label(tn$cd4_bb) <- "CD4 at baseline binary"
label(tn$cd4_db) <- "CD4 at diagnosis binary"
label(tn$birth_ende) <- "Birth country is Histo endemic country"
label(tn$mig_ende) <- "Migrated from a Histo endemic country"
label(tn$death) <- "Death count"
label(tn$art_pre) <- "Time from ART start date to Histo diagnosis"
label(tn$art_post) <- "Time from Histo diagnosis to ART start date"
label(tn$pre_post_histo) <- "ART diagnosis wrt ART start"


label(la$baseline) <- "Date of baseline"
label(la$SITE) <- "Site"
label(la$risk1) <- "HIV acquisition risk factor"
label(la$baseline_y) <- "Baseline Year"
label(la$cd4_b) <- "CD4 at Baseline"
label(la$cd4_diag) <- "CD4 at Histoplasmosis diagnosis"
label(la$rna_b) <- "HIV-RNA  at Baseline"
label(la$rna_b_log) <- "log HIV-RNA  at Baseline"
label(la$rna_diag) <- "HIV-RNA at Histoplasmosis diagnosis"
label(la$sex) <- "Sex at birth"
label(la$follow) <- "Follow-up time in years (from baseline to event/death)"
label(la$art_sd) <- "Date of ART initiation"
label(la$histo_d) <- "Date of Histoplasmosis diagnosis"
label(la$aids_first_d) <- "Date of HIV diagnosis"
label(la$l_alive_d) <- "Date of last visit"
label(la$death_d) <- "Date at death"
label(la$histo_y) <- "Year of Histoplasmosis diagnosis"
label(la$aids_first_y) <- "Year of HIV diagnosis"
label(la$l_alive_y) <- "Year of last visit"
label(la$death_y) <- "Year at death"
label(la$age_b) <- "Age at baseline"
label(la$naive) <- "ART-naive at cohort entry (yes/no)"
label(la$ade_type) <- "History of AIDS-defining event at clinic entry (other than histoplasmosis)"
label(la$rna_bb) <- "VL at baseline binary"
label(la$rna_db) <- "VL at diagnosis binary"
label(la$cd4_bb) <- "CD4 at baseline binary"
label(la$cd4_db) <- "CD4 at diagnosis binary"
label(la$tb) <- "Tb (Yes)"
label(la$cp_y) <- "Year at culture-negative TB diagnosis"
label(la$time_h) <- "Time till Histo diagnosis (days)"
label(la$Country) <- "Country of birth"
label(la$country_cat) <- "Migration status (Yes/No)"
label(la$rna_b_log) <- "Log of HIV-1 RNA at baseline"
label(la$rna_diag_log) <- "Log of HIV-1 RNA at diagnosis"
label(la$cd4_b_sq) <- "Square root transformation of CD4 at baseline"
label(la$cd4_diag_sq) <- "Square root transformation of CD4 at diagnosis"
label(la$birth_ende) <- "Birth country is Histo endemic country"
label(la$mig_ende) <- "Migrated from a Histo endemic country"
label(la$death) <- "Death count"
#label(la$mode_tb) <- "TB test type"
label(la$art_pre) <- "Time from ART start date to Histo diagnosis (days)"
label(la$art_post) <- "Time from Histo diagnosis to ART start date (days) "
label(la$pre_post_histo) <- "ART diagnosis wrt ART start"
save(la, file="la.Rdata")


#adding one day to histo date if the diagnosis date is same as enrollment date
la$histo_d <- with(la,ifelse(enrol_d == histo_d, histo_d+1,histo_d))
la$histo_d <- as.Date(la$histo_d, origin="1970-01-01")

#la <- la %>% filter(enrol_d < histo_d | is.na(histo_d))
la$end_d <-  ifelse(la$histo_num == 1, la$histo_d, la$last_d)
la$histo_d <- with(la,ifelse(enrol_d == histo_d, ymd(histo_d) + days(1), histo_d))
la$histo_d <- as.Date(la$histo_d, origin="1970-01-01")
la$follow <- with(la,ifelse(histo == "Yes",difftime(as.Date(histo_d),as.Date(baseline),units ="days")/365.25,
                            difftime(as.Date(last_d),as.Date(baseline),units ="days")/365.25))


la$end_d <- as.Date(la$end_d, origin="1970-01-01")
la$follow_time <- as.numeric(as.Date(la$end_d)- as.Date(la$enrol_d))
label(la$follow_time) <- "Follow-up time in days"

la$follow_time_yrs <- la$follow_time/ 365.25
label(la$follow_time_yrs) <- "Follow-up time in years"

incid <- la[, c("patient_id", "site", "enrol_d", "l_alive_d", "death_d", "histo_d", "follow_time_yrs","end_d","histo","histo_y")]
incid$enrol_y <- as.numeric(substr(incid$enrol_d, 1, 4))
incid$end_fu_y <- as.numeric(substr(incid$end_d, 1, 4))


yrs <- seq(2000, 2021, by=1)

for (i in 1:length(yrs)) {
  incid[[paste0("py_", yrs[i])]] <- ifelse(incid$enrol_y==yrs[i] & incid$end_fu_y==yrs[i], as.numeric(abs(incid$end_d - incid$enrol_d))/365.25,
                                           ifelse(incid$enrol_y==yrs[i] & incid$end_fu_y > yrs[i], 
                                                  as.numeric(abs(as.Date(paste0(yrs[i], "-12-31")) - as.Date(incid$enrol_d)))/365.25, 
                                                  ifelse(incid$enrol_y < yrs[i] & incid$end_fu_y==yrs[i], 
                                                         as.numeric(abs(as.Date(incid$end_d) - as.Date(paste0(yrs[i]-1, "-12-31"))))/365.25, 
                                                         ifelse(incid$enrol_y < yrs[i] & incid$end_fu_y> yrs[i], 1, 0))))
  
  incid[[paste0("num_persons_", yrs[i])]] <- ifelse(incid[[paste0("py_", yrs[i])]] >0, 1, 0)
}


incid_long <- incid[, c("patient_id", "site", grep("py_", names(incid), value=TRUE))]
incid_long <- melt(as.data.frame(incid_long), id.vars=c("patient_id", "site"), measure.vars= c(grep("py_", names(incid_long), value=TRUE)))

#value is person-years, so if value=0 this person should not have a record for that year
# since we are considering the window of diagnosis from -30 to date of enrollment, taking absolute value of "person years" and rounding them

incid_long <- subset(incid_long, value >0)

incid_long <- incid_long[order(incid_long$patient_id, incid_long$site), ]
incid_long$year <- as.numeric(substr(incid_long$variable, 4, 7))

incid2 <- incid[, c("patient_id", "enrol_d", "l_alive_d", "death_d", "end_d", "follow_time_yrs", "histo_d","site","histo","histo_y")]
incid_long <- merge(incid_long, incid2, all.x=TRUE)
incid_long$histo_y <- as.numeric(substr(incid_long$histo_d, 1, 4))
incid_long$histo <- ifelse(!is.na(incid_long$histo_d) & incid_long$histo_y == incid_long$year, 1, 0)

incid$histo_y <- as.numeric(substr(incid$histo_d, 1, 4))


#Person-years of follow-up & number of persons in database
#calculate incidence per year
incidence_histo <- incid %>% group_by(histo_y) %>%
  summarise(n=sum(!is.na(histo_d))) %>% ungroup %>% as.data.frame()


#number of persons is sum of 1 in data (1 if counted in year and 0 if not counted)
incidence_histo$num_persons <- NA
incidence_histo$py_follow_up <- NA

#removing the IDs with missing patient_id ID
incid <- incid %>% filter(!is.na(patient_id))

for (i in 1:length(yrs)) {
  incidence_histo[["num_persons"]][i] <- sum(incid[[paste0("num_persons_", yrs[i])]])
  incidence_histo[["py_follow_up"]][i] <- sum(incid[[paste0("py_", yrs[i])]])
}


incidence_histo$incidence_1000py <- incidence_histo$n/(incidence_histo$py_follow_up/1000)
incidence_histo <- incidence_histo %>% filter(!is.na(histo_y))
incidence_histo_total <- incidence_histo %>% summarise(n= sum(n),
                                                       num_persons= length(unique(la$patient_id)),
                                                       py_follow_up= sum(py_follow_up)) %>% ungroup()

incidence_histo_total$incidence_1000py <- incidence_histo_total$n/(incidence_histo_total$py_follow_up/1000)
incidence_histo_total$histo_y <- "Total"
incidence_histo_total <- incidence_histo_total[, c(5, 1:4)]

incidence_histo <- rbind(incidence_histo, incidence_histo_total)

## Incidents for Brazil with 
incid_b <- la[, c("patient_id", "site", "enrol_d", "l_alive_d", "death_d", "histo_d", "follow_time_yrs","end_d","histo","histo_y")]
incid_b <- incid_b %>% filter(site == "brazil")

incid_b$enrol_y <- as.numeric(substr(incid_b$enrol_d, 1, 4))
incid_b$end_fu_y <- as.numeric(substr(incid_b$end_d, 1, 4))


yrs <- seq(2000, 2021, by=1)

for (i in 1:length(yrs)) {
  incid_b[[paste0("py_", yrs[i])]] <- ifelse(incid_b$enrol_y==yrs[i] & incid_b$end_fu_y==yrs[i], as.numeric(abs(incid_b$end_d - incid_b$enrol_d))/365.25,
                                             ifelse(incid_b$enrol_y==yrs[i] & incid_b$end_fu_y > yrs[i], 
                                                    as.numeric(abs(as.Date(paste0(yrs[i], "-12-31")) - as.Date(incid_b$enrol_d)))/365.25, 
                                                    ifelse(incid_b$enrol_y < yrs[i] & incid_b$end_fu_y==yrs[i], 
                                                           as.numeric(abs(as.Date(incid_b$end_d) - as.Date(paste0(yrs[i]-1, "-12-31"))))/365.25, 
                                                           ifelse(incid_b$enrol_y < yrs[i] & incid_b$end_fu_y> yrs[i], 1, 0))))
  
  incid_b[[paste0("num_persons_", yrs[i])]] <- ifelse(incid_b[[paste0("py_", yrs[i])]] >0, 1, 0)
}


incid_b_long <- incid_b[, c("patient_id", "site", grep("py_", names(incid_b), value=TRUE))]
incid_b_long <- melt(as.data.frame(incid_b_long), id.vars=c("patient_id", "site"), measure.vars= c(grep("py_", names(incid_b_long), value=TRUE)))

#value is person-years, so if value=0 this person should not have a record for that year
# since we are considering the window of diagnosis from -30 to date of enrollment, taking absolute value of "person years" and rounding them

#incid_b_long$value <- round(abs(incid_b_long$value),digits =2)
#incid_b_long$value <- with(incid_b_long,ifelse(value<=0,))
incid_b_long <- subset(incid_b_long, value >0)

incid_b_long <- incid_b_long[order(incid_b_long$patient_id, incid_b_long$site), ]
incid_b_long$year <- as.numeric(substr(incid_b_long$variable, 4, 7))

incid_b2 <- incid_b[, c("patient_id", "enrol_d", "l_alive_d", "death_d", "end_d", "follow_time_yrs", "histo_d","site","histo","histo_y")]
incid_b_long <- merge(incid_b_long, incid_b2, all.x=TRUE)
incid_b_long$histo_y <- as.numeric(substr(incid_b_long$histo_d, 1, 4))
incid_b_long$histo <- ifelse(!is.na(incid_b_long$histo_d) & incid_b_long$histo_y == incid_b_long$year, 1, 0)

incid_b$histo_y <- as.numeric(substr(incid_b$histo_d, 1, 4))


#Person-years of follow-up & number of persons in database
#calculate incid_b per year
incid_b_histo <- incid_b %>% group_by(histo_y) %>%
  summarise(n=sum(!is.na(histo_d))) %>% ungroup %>% as.data.frame()


#number of persons is sum of 1 in data (1 if counted in year and 0 if not counted)
histo_y <- data.frame(c(2000:2021))
colnames(histo_y) <- "histo_y"
incid_b_histo <- left_join(histo_y,incid_b_histo,by="histo_y")
incid_b_histo$num_persons <- NA
incid_b_histo$py_follow_up <- NA

#removing the IDs with missing patient_id ID
incid_b <- incid_b %>% filter(!is.na(patient_id))

for (i in 1:length(yrs)) {
  incid_b_histo[["num_persons"]][i] <- sum(incid_b[[paste0("num_persons_", yrs[i])]])
  incid_b_histo[["py_follow_up"]][i] <- sum(incid_b[[paste0("py_", yrs[i])]])
}
la_b <- la %>% filter(site== "brazil")

incid_b_histo$incid_b_1000py <- incid_b_histo$n/(incid_b_histo$py_follow_up/1000)
incid_b_histo <- incid_b_histo %>% filter(!is.na(histo_y))
incid_b_histo_total <- incid_b_histo %>% summarise(n= sum(n,na.rm = T),
                                                   num_persons= length(unique(la_b$patient_id)),
                                                   py_follow_up= sum(py_follow_up,na.rm=T)) %>% ungroup()

incid_b_histo_total$incid_b_1000py <- incid_b_histo_total$n/(incid_b_histo_total$py_follow_up/1000)
incid_b_histo_total$histo_y <- "Total"
incid_b_histo_total <- incid_b_histo_total[, c(5, 1:4)]

incid_b_histo <- rbind(incid_b_histo, incid_b_histo_total)
incid_b_histo$n <- ifelse(is.na(incid_b_histo$n),0,incid_b_histo$n)
incid_b_histo$incid_b_1000py <- ifelse(is.na(incid_b_histo$incid_b_1000py),0,incid_b_histo$incid_b_1000py)
CI <- binconf(incid_b_histo$n,incid_b_histo$num_persons, method = "wilson")
CI <- as.data.frame(CI)
CI$histo_y <- c(2000:2021,"")
incid_b_histo <- left_join(incid_b_histo,CI, by="histo_y")
incid_b_histo$Lower <- ((incid_b_histo$Lower*incid_b_histo$num_persons)/incid_b_histo$py_follow_up)*1000
incid_b_histo$Upper <- ((incid_b_histo$Upper*incid_b_histo$num_persons)/incid_b_histo$py_follow_up)*1000


# For Mexico
incid_m <- la[, c("patient_id", "site", "enrol_d", "l_alive_d", "death_d", "histo_d", "follow_time_yrs","end_d","histo","histo_y")]
incid_m <- incid_m %>% filter(site == "mexico")

incid_m$enrol_y <- as.numeric(substr(incid_m$enrol_d, 1, 4))
incid_m$end_fu_y <- as.numeric(substr(incid_m$end_d, 1, 4))


yrs <- seq(2000, 2021, by=1)

for (i in 1:length(yrs)) {
  incid_m[[paste0("py_", yrs[i])]] <- ifelse(incid_m$enrol_y==yrs[i] & incid_m$end_fu_y==yrs[i], as.numeric(abs(incid_m$end_d - incid_m$enrol_d))/365.25,
                                             ifelse(incid_m$enrol_y==yrs[i] & incid_m$end_fu_y > yrs[i], 
                                                    as.numeric(abs(as.Date(paste0(yrs[i], "-12-31")) - as.Date(incid_m$enrol_d)))/365.25, 
                                                    ifelse(incid_m$enrol_y < yrs[i] & incid_m$end_fu_y==yrs[i], 
                                                           as.numeric(abs(as.Date(incid_m$end_d) - as.Date(paste0(yrs[i]-1, "-12-31"))))/365.25, 
                                                           ifelse(incid_m$enrol_y < yrs[i] & incid_m$end_fu_y> yrs[i], 1, 0))))
  
  incid_m[[paste0("num_persons_", yrs[i])]] <- ifelse(incid_m[[paste0("py_", yrs[i])]] >0, 1, 0)
}


incid_m_long <- incid_m[, c("patient_id", "site", grep("py_", names(incid_m), value=TRUE))]
incid_m_long <- melt(as.data.frame(incid_m_long), id.vars=c("patient_id", "site"), measure.vars= c(grep("py_", names(incid_m_long), value=TRUE)))

incid_m_long <- subset(incid_m_long, value >0)

incid_m_long <- incid_m_long[order(incid_m_long$patient_id, incid_m_long$site), ]
incid_m_long$year <- as.numeric(substr(incid_m_long$variable, 4, 7))

incid_m2 <- incid_m[, c("patient_id", "enrol_d", "l_alive_d", "death_d", "end_d", "follow_time_yrs", "histo_d","site","histo","histo_y")]
incid_m_long <- merge(incid_m_long, incid_m2, all.x=TRUE)
incid_m_long$histo_y <- as.numeric(substr(incid_m_long$histo_d, 1, 4))
incid_m_long$histo <- ifelse(!is.na(incid_m_long$histo_d) & incid_m_long$histo_y == incid_m_long$year, 1, 0)

incid_m$histo_y <- as.numeric(substr(incid_m$histo_d, 1, 4))


#Person-years of follow-up & number of persons in database
#calculate incid_m per year
incid_m_histo <- incid_m %>% group_by(histo_y) %>%
  summarise(n=sum(!is.na(histo_d))) %>% ungroup %>% as.data.frame()


#number of persons is sum of 1 in data (1 if counted in year and 0 if not counted)
histo_y <- data.frame(c(2000:2021))
colnames(histo_y) <- "histo_y"
incid_m_histo <- left_join(histo_y,incid_m_histo,by="histo_y")
incid_m_histo$num_persons <- NA
incid_m_histo$py_follow_up <- NA

#removing the IDs with missing patient_id ID
incid_m <- incid_m %>% filter(!is.na(patient_id))

for (i in 1:length(yrs)) {
  incid_m_histo[["num_persons"]][i] <- sum(incid_m[[paste0("num_persons_", yrs[i])]])
  incid_m_histo[["py_follow_up"]][i] <- sum(incid_m[[paste0("py_", yrs[i])]])
  
}

la_m <- la %>% filter(site== "mexico")
incid_m_histo$incid_m_1000py <- incid_m_histo$n/(incid_m_histo$py_follow_up/1000)
incid_m_histo <- incid_m_histo %>% filter(!is.na(histo_y))
incid_m_histo_total <- incid_m_histo %>% summarise(n= sum(n,na.rm = T),
                                                   num_persons= length(unique(la_m$patient_id)),
                                                   py_follow_up= sum(py_follow_up,na.rm = T)) %>% ungroup()

incid_m_histo_total$incid_m_1000py <- incid_m_histo_total$n/(incid_m_histo_total$py_follow_up/1000)
incid_m_histo_total$histo_y <- "Total"
incid_m_histo_total <- incid_m_histo_total[, c(5, 1:4)]

incid_m_histo <- rbind(incid_m_histo, incid_m_histo_total)
incid_m_histo$n <- ifelse(is.na(incid_m_histo$n),0,incid_m_histo$n)
incid_m_histo$incid_m_1000py <- ifelse(is.na(incid_m_histo$incid_m_1000py),0,incid_m_histo$incid_m_1000py)
CI <- binconf(incid_m_histo$n,incid_m_histo$num_persons, method = "wilson")
CI <- as.data.frame(CI)
CI$histo_y <- c(2000:2021,"")
incid_m_histo <- left_join(incid_m_histo,CI, by="histo_y")
incid_m_histo$Lower <- ((incid_m_histo$Lower*incid_m_histo$num_persons)/incid_m_histo$py_follow_up)*1000
incid_m_histo$Upper <- ((incid_m_histo$Upper*incid_m_histo$num_persons)/incid_m_histo$py_follow_up)*1000


save(incidence_histo, file="incidence_histo.Rdata")
save(incid_b_histo, file="incid_b_histo.Rdata")
save(incid_m_histo, file="incid_m_histo.Rdata")
save(incid_b_long, file="incid_b_long.Rdata")
save(incid_m_long, file="incid_m_long.Rdata")

# For Honduras 
incid_h <- la[, c("patient_id", "site", "enrol_d", "l_alive_d", "death_d", "histo_d", "follow_time_yrs","end_d","histo","histo_y")]
incid_h <- incid_h %>% filter(site == "honduras")

incid_h$enrol_y <- as.numeric(substr(incid_h$enrol_d, 1, 4))
incid_h$end_fu_y <- as.numeric(substr(incid_h$end_d, 1, 4))


yrs <- seq(2000, 2021, by=1)

for (i in 1:length(yrs)) {
  incid_h[[paste0("py_", yrs[i])]] <- ifelse(incid_h$enrol_y==yrs[i] & incid_h$end_fu_y==yrs[i], as.numeric(abs(incid_h$end_d - incid_h$enrol_d))/365.25,
                                             ifelse(incid_h$enrol_y==yrs[i] & incid_h$end_fu_y > yrs[i], 
                                                    as.numeric(abs(as.Date(paste0(yrs[i], "-12-31")) - as.Date(incid_h$enrol_d)))/365.25, 
                                                    ifelse(incid_h$enrol_y < yrs[i] & incid_h$end_fu_y==yrs[i], 
                                                           as.numeric(abs(as.Date(incid_h$end_d) - as.Date(paste0(yrs[i]-1, "-12-31"))))/365.25, 
                                                           ifelse(incid_h$enrol_y < yrs[i] & incid_h$end_fu_y> yrs[i], 1, 0))))
  
  incid_h[[paste0("num_persons_", yrs[i])]] <- ifelse(incid_h[[paste0("py_", yrs[i])]] >0, 1, 0)
}


incid_h_long <- incid_h[, c("patient_id", "site", grep("py_", names(incid_h), value=TRUE))]
incid_h_long <- melt(as.data.frame(incid_h_long), id.vars=c("patient_id", "site"), measure.vars= c(grep("py_", names(incid_h_long), value=TRUE)))

incid_h_long <- subset(incid_h_long, value >0)

incid_h_long <- incid_h_long[order(incid_h_long$patient_id, incid_h_long$site), ]
incid_h_long$year <- as.numeric(substr(incid_h_long$variable, 4, 7))

incid_h2 <- incid_h[, c("patient_id", "enrol_d", "l_alive_d", "death_d", "end_d", "follow_time_yrs", "histo_d","site","histo","histo_y")]
incid_h_long <- merge(incid_h_long, incid_h2, all.x=TRUE)
incid_h_long$histo_y <- as.numeric(substr(incid_h_long$histo_d, 1, 4))
incid_h_long$histo <- ifelse(!is.na(incid_h_long$histo_d) & incid_h_long$histo_y == incid_h_long$year, 1, 0)

incid_h$histo_y <- as.numeric(substr(incid_h$histo_d, 1, 4))


#Person-years of follow-up & number of persons in database
#calculate incid_h per year
incid_h_histo <- incid_h %>% group_by(histo_y) %>%
  summarise(n=sum(!is.na(histo_d))) %>% ungroup %>% as.data.frame()


#number of persons is sum of 1 in data (1 if counted in year and 0 if not counted)
histo_y <- data.frame(c(2000:2021))
colnames(histo_y) <- "histo_y"
incid_h_histo <- left_join(histo_y,incid_h_histo,by="histo_y")
incid_h_histo$num_persons <- NA
incid_h_histo$py_follow_up <- NA

#removing the IDs with missing patient_id ID
incid_h <- incid_h %>% filter(!is.na(patient_id))

for (i in 1:length(yrs)) {
  incid_h_histo[["num_persons"]][i] <- sum(incid_h[[paste0("num_persons_", yrs[i])]])
  incid_h_histo[["py_follow_up"]][i] <- sum(incid_h[[paste0("py_", yrs[i])]])
  
}

la_h <- la %>% filter(site== "honduras")

incid_h_histo$incid_h_1000py <- incid_h_histo$n/(incid_h_histo$py_follow_up/1000)
incid_h_histo <- incid_h_histo %>% filter(!is.na(histo_y))
incid_h_histo_total <- incid_h_histo %>% summarise(n= sum(n,na.rm = T),
                                                   num_persons= length(unique(la_h$patient_id)),
                                                   py_follow_up= sum(py_follow_up, na.rm = T)) %>% ungroup()

incid_h_histo_total$incid_h_1000py <- incid_h_histo_total$n/(incid_h_histo_total$py_follow_up/1000)
incid_h_histo_total$histo_y <- "Total"
incid_h_histo_total <- incid_h_histo_total[, c(5, 1:4)]

incid_h_histo <- rbind(incid_h_histo, incid_h_histo_total)
incid_h_histo$n <- ifelse(is.na(incid_h_histo$n),0,incid_h_histo$n)
incid_h_histo$incid_h_1000py <- ifelse(is.na(incid_h_histo$incid_h_1000py),0,incid_h_histo$incid_h_1000py)
CI <- binconf(incid_h_histo$n,incid_h_histo$num_persons, method = "wilson")
CI <- as.data.frame(CI)
CI$histo_y <- c(2000:2021,"")
incid_h_histo <- left_join(incid_h_histo,CI, by="histo_y")
incid_h_histo$Lower <- ((incid_h_histo$Lower*incid_h_histo$num_persons)/incid_h_histo$py_follow_up)*1000
incid_h_histo$Upper <- ((incid_h_histo$Upper*incid_h_histo$num_persons)/incid_h_histo$py_follow_up)*1000



# For Peru 
incid_p <- la[, c("patient_id", "site", "enrol_d", "l_alive_d", "death_d", "histo_d", "follow_time_yrs","end_d","histo","histo_y")]
incid_p <- incid_p %>% filter(site == "peru")

incid_p$enrol_y <- as.numeric(substr(incid_p$enrol_d, 1, 4))
incid_p$end_fu_y <- as.numeric(substr(incid_p$end_d, 1, 4))


yrs <- seq(2000, 2019, by=1)

for (i in 1:length(yrs)) {
  incid_p[[paste0("py_", yrs[i])]] <- ifelse(incid_p$enrol_y==yrs[i] & incid_p$end_fu_y==yrs[i], as.numeric(abs(incid_p$end_d - incid_p$enrol_d))/365.25,
                                             ifelse(incid_p$enrol_y==yrs[i] & incid_p$end_fu_y > yrs[i], 
                                                    as.numeric(abs(as.Date(paste0(yrs[i], "-12-31")) - as.Date(incid_p$enrol_d)))/365.25, 
                                                    ifelse(incid_p$enrol_y < yrs[i] & incid_p$end_fu_y==yrs[i], 
                                                           as.numeric(abs(as.Date(incid_p$end_d) - as.Date(paste0(yrs[i]-1, "-12-31"))))/365.25, 
                                                           ifelse(incid_p$enrol_y < yrs[i] & incid_p$end_fu_y> yrs[i], 1, 0))))
  
  incid_p[[paste0("num_persons_", yrs[i])]] <- ifelse(incid_p[[paste0("py_", yrs[i])]] >0, 1, 0)
}


incid_p_long <- incid_p[, c("patient_id", "site", grep("py_", names(incid_p), value=TRUE))]
incid_p_long <- melt(as.data.frame(incid_p_long), id.vars=c("patient_id", "site"), measure.vars= c(grep("py_", names(incid_p_long), value=TRUE)))

incid_p_long <- subset(incid_p_long, value >0)

incid_p_long <- incid_p_long[order(incid_p_long$patient_id, incid_p_long$site), ]
incid_p_long$year <- as.numeric(substr(incid_p_long$variable, 4, 7))

incid_p2 <- incid_p[, c("patient_id", "enrol_d", "l_alive_d", "death_d", "end_d", "follow_time_yrs", "histo_d","site","histo","histo_y")]
incid_p_long <- merge(incid_p_long, incid_p2, all.x=TRUE)
incid_p_long$histo_y <- as.numeric(substr(incid_p_long$histo_d, 1, 4))
incid_p_long$histo <- ifelse(!is.na(incid_p_long$histo_d) & incid_p_long$histo_y == incid_p_long$year, 1, 0)

incid_p$histo_y <- as.numeric(substr(incid_p$histo_d, 1, 4))


#Person-years of follow-up & number of persons in database
#calculate incid_p per year
incid_p_histo <- incid_p %>% group_by(histo_y) %>%
  summarise(n=sum(!is.na(histo_d))) %>% ungroup %>% as.data.frame()


#number of persons is sum of 1 in data (1 if counted in year and 0 if not counted)
histo_y <- data.frame(c(2000:2019))
colnames(histo_y) <- "histo_y"
incid_p_histo <- left_join(histo_y,incid_p_histo,by="histo_y")
incid_p_histo$num_persons <- NA
incid_p_histo$py_follow_up <- NA

#removing the IDs with missing patient_id ID
incid_p <- incid_p %>% filter(!is.na(patient_id))

for (i in 1:length(yrs)) {
  incid_p_histo[["num_persons"]][i] <- sum(incid_p[[paste0("num_persons_", yrs[i])]])
  incid_p_histo[["py_follow_up"]][i] <- sum(incid_p[[paste0("py_", yrs[i])]])
  
}

la_p <- la %>% filter(site== "peru")
incid_p_histo$incid_p_1000py <- incid_p_histo$n/(incid_p_histo$py_follow_up/1000)
incid_p_histo <- incid_p_histo %>% filter(!is.na(histo_y))
incid_p_histo_total <- incid_p_histo %>% summarise(n= sum(n,na.rm = T),
                                                   num_persons= length(unique(la_p$patient_id)),
                                                   py_follow_up= sum(py_follow_up,na.rm = T)) %>% ungroup()

incid_p_histo_total$incid_p_1000py <- incid_p_histo_total$n/(incid_p_histo_total$py_follow_up/1000)
incid_p_histo_total$histo_y <- "Total"
incid_p_histo_total <- incid_p_histo_total[, c(5, 1:4)]

incid_p_histo <- rbind(incid_p_histo, incid_p_histo_total)
incid_p_histo$n <- ifelse(is.na(incid_p_histo$n),0,incid_p_histo$n)
incid_p_histo$incid_p_1000py <- ifelse(is.na(incid_p_histo$incid_p_1000py),0,incid_p_histo$incid_p_1000py)
CI <- binconf(incid_p_histo$n,incid_p_histo$num_persons, method = "wilson")
CI <- as.data.frame(CI)
CI$histo_y <- c(2000:2019,"")
incid_p_histo <- left_join(incid_p_histo,CI, by="histo_y")
incid_p_histo$Lower <- ((incid_p_histo$Lower*incid_p_histo$num_persons)/incid_p_histo$py_follow_up)*1000
incid_p_histo$Upper <- ((incid_p_histo$Upper*incid_p_histo$num_persons)/incid_p_histo$py_follow_up)*1000


save(incid_h_histo, file="incid_h_histo.Rdata")
save(incid_p_histo, file="incid_p_histo.Rdata")
save(incid_h_long, file="incid_h_long.Rdata")
save(incid_p_long, file="incid_p_long.Rdata")

# For Chile 
incid_c <- la[, c("patient_id", "site", "enrol_d", "l_alive_d", "death_d", "histo_d", "follow_time_yrs","end_d","histo","histo_y")]
incid_c <- incid_c %>% filter(site == "chile")

incid_c$enrol_y <- as.numeric(substr(incid_c$enrol_d, 1, 4))
incid_c$end_fu_y <- as.numeric(substr(incid_c$end_d, 1, 4))


yrs <- seq(2000, 2021, by=1)

for (i in 1:length(yrs)) {
  incid_c[[paste0("py_", yrs[i])]] <- ifelse(incid_c$enrol_y==yrs[i] & incid_c$end_fu_y==yrs[i], as.numeric(abs(incid_c$end_d - incid_c$enrol_d))/365.25,
                                             ifelse(incid_c$enrol_y==yrs[i] & incid_c$end_fu_y > yrs[i], 
                                                    as.numeric(abs(as.Date(paste0(yrs[i], "-12-31")) - as.Date(incid_c$enrol_d)))/365.25, 
                                                    ifelse(incid_c$enrol_y < yrs[i] & incid_c$end_fu_y==yrs[i], 
                                                           as.numeric(abs(as.Date(incid_c$end_d) - as.Date(paste0(yrs[i]-1, "-12-31"))))/365.25, 
                                                           ifelse(incid_c$enrol_y < yrs[i] & incid_c$end_fu_y> yrs[i], 1, 0))))
  
  incid_c[[paste0("num_persons_", yrs[i])]] <- ifelse(incid_c[[paste0("py_", yrs[i])]] >0, 1, 0)
}


incid_c_long <- incid_c[, c("patient_id", "site", grep("py_", names(incid_c), value=TRUE))]
incid_c_long <- melt(as.data.frame(incid_c_long), id.vars=c("patient_id", "site"), measure.vars= c(grep("py_", names(incid_c_long), value=TRUE)))

incid_c_long <- subset(incid_c_long, value >0)

incid_c_long <- incid_c_long[order(incid_c_long$patient_id, incid_c_long$site), ]
incid_c_long$year <- as.numeric(substr(incid_c_long$variable, 4, 7))

incid_c2 <- incid_c[, c("patient_id", "enrol_d", "l_alive_d", "death_d", "end_d", "follow_time_yrs", "histo_d","site","histo","histo_y")]
incid_c_long <- merge(incid_c_long, incid_c2, all.x=TRUE)
incid_c_long$histo_y <- as.numeric(substr(incid_c_long$histo_d, 1, 4))
incid_c_long$histo <- ifelse(!is.na(incid_c_long$histo_d) & incid_c_long$histo_y == incid_c_long$year, 1, 0)

incid_c$histo_y <- as.numeric(substr(incid_c$histo_d, 1, 4))


#Person-years of follow-up & number of persons in database
#calculate incid_c per year
incid_c_histo <- incid_c %>% group_by(histo_y) %>%
  summarise(n=sum(!is.na(histo_d))) %>% ungroup %>% as.data.frame()


#number of persons is sum of 1 in data (1 if counted in year and 0 if not counted)
histo_y <- data.frame(c(2000:2021))
colnames(histo_y) <- "histo_y"
incid_c_histo <- left_join(histo_y,incid_c_histo,by="histo_y")
incid_c_histo$num_persons <- NA
incid_c_histo$py_follow_up <- NA

#removing the IDs with missing patient_id ID
incid_c <- incid_c %>% filter(!is.na(patient_id))

for (i in 1:length(yrs)) {
  incid_c_histo[["num_persons"]][i] <- sum(incid_c[[paste0("num_persons_", yrs[i])]])
  incid_c_histo[["py_follow_up"]][i] <- sum(incid_c[[paste0("py_", yrs[i])]])
  
}

la_c <- la %>% filter(site== "chile")

incid_c_histo$incid_c_1000py <- incid_c_histo$n/(incid_c_histo$py_follow_up/1000)
#incid_c_histo <- incid_c_histo %>% filter(!is.na(histo_y))
incid_c_histo_total <- incid_c_histo %>% summarise(n= sum(n, na.rm = T),
                                                   num_persons= length(unique(la_c$patient_id)),
                                                   py_follow_up= sum(py_follow_up, na.rm = T)) %>% ungroup()

incid_c_histo_total$incid_c_1000py <- incid_c_histo_total$n/(incid_c_histo_total$py_follow_up/1000)
incid_c_histo_total$histo_y <- "Total"
incid_c_histo_total <- incid_c_histo_total[, c(5, 1:4)]

incid_c_histo <- rbind(incid_c_histo, incid_c_histo_total)
incid_c_histo$n <- ifelse(is.na(incid_c_histo$n),0,incid_c_histo$n)
incid_c_histo$incid_c_1000py <- ifelse(is.na(incid_c_histo$incid_c_1000py),0,incid_c_histo$incid_c_1000py)
CI <- binconf(incid_c_histo$n,incid_c_histo$num_persons, method = "wilson")
CI <- as.data.frame(CI)
CI$histo_y <- c(2000:2021,"")
incid_c_histo <- left_join(incid_c_histo,CI, by="histo_y")
incid_c_histo$Lower <- ((incid_c_histo$Lower*incid_c_histo$num_persons)/incid_c_histo$py_follow_up)*1000
incid_c_histo$Upper <- ((incid_c_histo$Upper*incid_c_histo$num_persons)/incid_c_histo$py_follow_up)*1000


save(incid_c_histo, file="incid_c_histo.Rdata")
save(incid_c_long, file="incid_c_long.Rdata")
#TN ##
tn$enrol_d <- ifelse(!is.na(tn$year_of_enrollment),
                     paste0(tn$year_of_enrollment, "-07-01"),NA)
tn$enrol_d <- as.Date(tn$enrol_d)
tn$age_at_histo <- with(tn,ifelse(age_at_histo == age_at_first_visit,age_at_histo + 0.002737851,age_at_histo))
tn$diff_enrol_histo <- with(tn,(age_at_histo - age_at_first_visit)*365.25)
tn$histo_d <- with(tn,enrol_d + diff_enrol_histo)
tn$age_last <- ifelse(tn$CENSOR_TYPE == "LAST VISIT", tn$age_at_last_visit,tn$age_at_censor)
tn$diff_enrol_last <- with(tn,(age_last - age_at_first_visit)*365.25)
tn$last_d <- with(tn,enrol_d + diff_enrol_last)
# tn$last_d <- with(tn,ifelse(is.na(last_d),enrol_d + diff_enrol_last, last_d))
tn$last_d <- with(tn,ifelse(CENSOR_TYPE == "12/31/2021", as.Date("2021-12-31"),as.Date(last_d,origin="1970-01-01")))
tn$last_d <- as.Date(tn$last_d,origin="1970-01-01")
tn$last_d <- with(tn,ifelse(last_d > "2021-12-31",as.Date("2021-12-31"),last_d))
tn$last_d <- as.Date(tn$last_d,origin="1970-01-01")
tn$enrol_y <- tn$year_of_enrollment

# tn$end_d <- ifelse(tn$histo == "Yes" , tn$age_at_histo, 
#                    ifelse(tn$CENSOR_TYPE == "Death" | tn$CENSOR_TYPE == "12/31/2021", 
#                           min(tn$age_at_last_visit,tn$age_at_censor),tn$age_at_last_visit))

tn$end_d <-  ifelse(tn$histo_num == 1, tn$histo_d, tn$last_d)
tn$baseline <- tn$enrol_d
tn$follow <- with(tn,ifelse(histo == "Yes",difftime(as.Date(histo_d),as.Date(baseline),units ="days")/365.25,
                            difftime(as.Date(last_d),as.Date(baseline),units ="days")/365.25))


tn$end_d <- as.Date(tn$end_d, origin="1970-01-01")
tn$follow_time <- as.numeric(as.Date(tn$end_d)- as.Date(tn$enrol_d))
label(tn$follow_time) <- "Follow-up time in days"

tn$follow_time_yrs <- tn$follow_time/ 365.25
label(tn$follow_time_yrs) <- "Follow-up time in years"

incid_tn <- tn[, c("CFAR_PID", "age_at_first_visit", "age_at_last_visit", "age_at_death", "age_at_histo", 
                   "follow_time_yrs","end_d","year_of_enrollment","histo_d","enrol_d","histo","enrol_y")]
incid_tn$end_fu_y <- as.numeric(substr(incid_tn$end_d, 1, 4))

yrs <- seq(2000, 2021, by=1)

for (i in 1:length(yrs)) {
  incid_tn[[paste0("py_", yrs[i])]] <- ifelse(incid_tn$enrol_y==yrs[i] & incid_tn$end_fu_y==yrs[i], as.numeric(abs(incid_tn$end_d - incid_tn$enrol_d))/365.25,
                                              ifelse(incid_tn$enrol_y==yrs[i] & incid_tn$end_fu_y > yrs[i], 
                                                     as.numeric(abs(as.Date(paste0(yrs[i], "-12-31")) - as.Date(incid_tn$enrol_d)))/365.25, 
                                                     ifelse(incid_tn$enrol_y < yrs[i] & incid_tn$end_fu_y==yrs[i], 
                                                            as.numeric(abs(as.Date(incid_tn$end_d) - as.Date(paste0(yrs[i]-1, "-12-31"))))/365.25, 
                                                            ifelse(incid_tn$enrol_y < yrs[i] & incid_tn$end_fu_y> yrs[i], 1, 0))))
  
  incid_tn[[paste0("num_persons_", yrs[i])]] <- ifelse(incid_tn[[paste0("py_", yrs[i])]] >0, 1, 0)
}

incid_tn_long <- incid_tn[, c("CFAR_PID", grep("py_", names(incid_tn), value=TRUE))]
incid_tn_long <- melt(as.data.frame(incid_tn_long), id.vars=c("CFAR_PID"), measure.vars= c(grep("py_", names(incid_tn_long), value=TRUE)))

#value is person-years, so if value=0 this person should not have a record for that year
#since we are considering the window of diagnosis from -30 to date of enrollment, taking absolute value of "person years" and rounding them

#incid_tn_long$value <- round(abs(incid_tn_long$value),digits =2)
#incid_tn_long$value <- with(incid_tn_long,ifelse(value<=0,))
incid_tn_long <- subset(incid_tn_long, value > 0)

incid_tn_long <- incid_tn_long[order(incid_tn_long$CFAR_PID), ]
incid_tn_long$year <- as.numeric(substr(incid_tn_long$variable, 4, 7))

incid_tn2 <- incid_tn[, c("CFAR_PID", "age_at_first_visit", "age_at_last_visit", "age_at_death", "end_d", 
                          "follow_time_yrs", "age_at_histo","year_of_enrollment","histo_d")]
incid_tn_long <- merge(incid_tn_long, incid_tn2, all.x=TRUE)
incid_tn_long$histo_y <- as.numeric(substr(incid_tn_long$histo_d, 1, 4))
#incid_tn_long$histo_y <- incid_tn_long$year_of_enrollment + round((incid_tn_long$age_at_histo - incid_tn_long$age_at_first_visit),digits=0)
incid_tn_long$histo <- ifelse(!is.na(incid_tn_long$histo_d) & incid_tn_long$histo_y == incid_tn_long$year, 1, 0)

incid_tn$histo_y <- as.numeric(substr(incid_tn$histo_d, 1, 4))
#Person-years of follow-up & number of persons in database
#calculate incid_tn per year

incidence_tn_histo <- incid_tn %>% group_by(histo_y) %>%
  summarise(n=sum(!is.na(histo_d))) %>% ungroup %>% as.data.frame()

#number of persons is sum of 1 in data (1 if counted in year and 0 if not counted)
incidence_tn_histo$num_persons <- NA
incidence_tn_histo$py_follow_up <- NA

#removing the IDs with missing CFAR_PID ID
incid_tn <- incid_tn %>% filter(!is.na(CFAR_PID))

for (i in 1:length(yrs)) {
  incidence_tn_histo[["num_persons"]][i] <- sum(incid_tn[[paste0("num_persons_", yrs[i])]])
  incidence_tn_histo[["py_follow_up"]][i] <- sum(incid_tn[[paste0("py_", yrs[i])]])
}

incidence_tn_histo$incidence_tn_1000py <- incidence_tn_histo$n/(incidence_tn_histo$py_follow_up/1000)

## There was last row with missing so removing that row
incidence_tn_histo <- incidence_tn_histo %>% filter(!is.na(incidence_tn_histo$histo_y))
incidence_tn_histo_total <- incidence_tn_histo %>% summarise(n= sum(n),
                                                             num_persons= length(unique(tn$CFAR_PID)),
                                                             py_follow_up= sum(py_follow_up)) %>% ungroup()

incidence_tn_histo_total$incidence_tn_1000py <- incidence_tn_histo_total$n/(incidence_tn_histo_total$py_follow_up/1000)
incidence_tn_histo_total$histo_y<- "Total"
incidence_tn_histo_total <- incidence_tn_histo_total[, c(5, 1:4)]

incidence_tn_histo <- rbind(incidence_tn_histo, incidence_tn_histo_total)
save(incidence_tn_histo, file="incidence_tn_histo.Rdata")
save(incid_long, file="incid_long.Rdata")
save(incid_tn_long, file="incid_tn_long.Rdata")
label(tn$follow_time_yrs) <- "Follow-up time in years (from baseline to event/death)"
save(tn, file="tn.Rdata")

## ----------------- ###
## Regression Models ###
## ----------------- ###

la <- data.frame(la)
la$follow <- with(la,ifelse(follow ==0,0.002737851,la$follow ))
la$follow <- abs(la$follow)

# Relevel
la$sex <- relevel(as.factor(la$sex), ref = "Female")
la$SITE <- relevel(as.factor(la$SITE), ref = "Brazil")
la$risk1 <- relevel(la$risk1, ref = "MSM")
la$rna_bb <- relevel(la$rna_bb, ref = "Undetectable")
la$cd4_bb <- relevel(la$cd4_bb, ref = "<=350 cells/ mm3")

## Non spline univariate models
uni.gen <- glm(histo_num ~as.factor(sex) + offset(log(follow)),family = poisson(link = "log"), data = la)
uni.risk <- glm(histo_num ~as.factor(risk1) + offset(log(follow)),family = poisson(link = "log"), data = la)
uni.rna <- glm(histo_num ~as.factor(rna_bb) + offset(log(follow)),family = poisson(link = "log"), data = la)
uni.site <- glm(histo_num ~as.factor(SITE) + offset(log(follow)),family = poisson(link = "log"), data = la)
#uni.multi <- glm(histo_num ~ age_b +  as.factor(cd4_bb)+ as.factor(sex)  + as.factor(risk) + as.factor(rna_bb) +  offset(log(follow)) + as.factor(SITE),family = poisson(link = "log"), data = la)

## Univariate spline models

ages <- round(min(la$age_b)):round(max(la$age_b))
sp <- as.data.frame(ns(ages, df=4))
sp <- cbind(ages, sp) 
sp <- as.data.frame(sp)
rownames(sp) <- sp$age
sp <- sp[, -1]
colnames(sp) <- c("X1", "X2", "X3", "X4")

ms1 <- glm(histo_num ~ ns(age_b,df =4)+ offset(log(follow)), la, family= poisson(link="log"))
summ_spline  <- summary(ms1)
v_spline <- sandwich(ms1)
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Coefs", "SE", "Variables")

rr_age_20 <- rr(20,summ_spline)
rr_age_25 <- rr(25,summ_spline)
rr_age_30 <- rr(30,summ_spline)
rr_age_40 <- rr(40,summ_spline)
rr_age_45 <- rr(45,summ_spline)
rr_age_50 <- rr(50,summ_spline)
rr_age_55 <- rr(55,summ_spline)
rr_age_60 <- rr(60,summ_spline)


SE_20 <- spline_se_fun(20)
SE_25 <- spline_se_fun(25)
SE_30 <- spline_se_fun(30)
SE_40 <- spline_se_fun(40)
SE_45 <- spline_se_fun(45)
SE_50 <- spline_se_fun(50)
SE_55 <- spline_se_fun(55)
SE_60 <- spline_se_fun(60)

CI <- function(rr,SE) 
{ exp(log(rr) + c(-1,1)*1.96*SE)}


CI_20 <- CI(rr_age_20,SE_20)
CI_25 <- CI(rr_age_25,SE_25)
CI_30 <- CI(rr_age_30,SE_30)
CI_40 <- CI(rr_age_40,SE_40)
CI_45 <- CI(rr_age_45,SE_45)
CI_50 <- CI(rr_age_50,SE_50)
CI_55 <- CI(rr_age_55,SE_55)
CI_60 <- CI(rr_age_60,SE_60)




abs_20 <- abs_val(20,summ_spline)
abs_25 <- abs_val(25,summ_spline)
abs_30 <- abs_val(30,summ_spline)
abs_40 <- abs_val(40,summ_spline)
abs_45 <- abs_val(45,summ_spline)
abs_50 <- abs_val(50,summ_spline)
abs_55 <- abs_val(55,summ_spline)
abs_60 <- abs_val(60,summ_spline)



age_spline <- data.frame("RR" = c(rr_age_20,rr_age_25,rr_age_30,rr_age_40,rr_age_45,rr_age_50,rr_age_55,rr_age_60),
                         "Lower CI"= c(CI_20[1],CI_25[1],CI_30[1],CI_40[1],CI_45[1],CI_50[1],CI_55[1],CI_60[1]),
                         "Upper CI"= c(CI_20[2],CI_25[2],CI_30[2],CI_40[2],CI_45[2],CI_50[2],CI_55[2],CI_60[2]),
                         "p-value" = c(2*(1-pnorm(abs(abs_20/SE_20))), 2*(1-pnorm(abs(abs_25/SE_25))), 2*(1-pnorm(abs(abs_30/SE_30))),
                                       2*(1-pnorm(abs(abs_40/SE_40))), 2*(1-pnorm(abs(abs_45/SE_45))), 2*(1-pnorm(abs(abs_50/SE_50))),
                                       2*(1-pnorm(abs(abs_55/SE_55))), 2*(1-pnorm(abs(abs_60/SE_60)))))
colnames(age_spline) <- c("RR","Lower CI","Upper CI","p-value")



cd4s <- round(min(la$cd4_b,na.rm = TRUE)):round(max(la$cd4_b,na.rm = TRUE))
sp_cd4 <- as.data.frame(ns(cd4s, df=4))
sp_cd4 <- cbind(cd4s, sp_cd4) 
sp_cd4 <- as.data.frame(sp_cd4)
rownames(sp_cd4) <- sp_cd4$cd4
sp_cd4 <- sp_cd4[, -1]
colnames(sp_cd4) <- c("X1", "X2", "X3", "X4")

ms1 <- glm(histo_num ~ ns(cd4_b,df=4)+ offset(log(follow)),family = poisson(link = "log"), data = la)
summ_spline  <- summary(ms1)
v_spline <- sandwich(ms1)
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Coefs", "SE", "Variables")

rr_cd4_100 <- rr_cd4(100,summ_spline)
rr_cd4_500 <- rr_cd4(500,summ_spline)
rr_cd4_1500 <- rr_cd4(1500,summ_spline)
rr_cd4_2000 <- rr_cd4(2000,summ_spline)
rr_cd4_2500 <- rr_cd4(2500,summ_spline)



SE_100 <- spline_se_fun_cd(100)
SE_500 <- spline_se_fun_cd(500)
SE_1500 <- spline_se_fun_cd(1500)
SE_2000 <- spline_se_fun_cd(2000)
SE_2500 <- spline_se_fun_cd(2500)


CI_100 <- CI(rr_cd4_100,SE_100)
CI_500 <- CI(rr_cd4_500,SE_500)
CI_1500 <- CI(rr_cd4_1500,SE_1500)
CI_2000 <- CI(rr_cd4_2000,SE_2000)
CI_2500 <- CI(rr_cd4_2500,SE_2500)

abs_100 <- abs_val_cd(100,summ_spline)
abs_500 <- abs_val_cd(500,summ_spline)
abs_1500 <- abs_val_cd(1500,summ_spline)
abs_2000 <- abs_val_cd(2000,summ_spline)
abs_2500 <- abs_val_cd(2500,summ_spline)



cd4_spline <- data.frame("RR" = c(rr_cd4_100,rr_cd4_500,rr_cd4_1500,rr_cd4_2000,rr_cd4_2500),
                         "Lower CI"= c(CI_100[1],CI_500[1],CI_1500[1],CI_2000[1],CI_2500[1]),
                         "Upper CI"= c(CI_100[2],CI_500[2],CI_1500[2],CI_2000[2],CI_2500[2]),
                         "p-value" = c(2*(1-pnorm(abs(abs_100/SE_100))), 2*(1-pnorm(abs(abs_500/SE_500))), 2*(1-pnorm(abs(abs_1500/SE_1500))),
                                       2*(1-pnorm(abs(abs_2000/SE_2000))), 2*(1-pnorm(abs(abs_2500/SE_2500)))))
colnames(cd4_spline) <- c("RR","Lower CI","Upper CI","p-value")


uni.tab1 <- rbind(age_spline ,cd4_spline,uni.fit(uni.gen),uni.fit3(uni.risk),uni.fit(uni.rna),uni.fit3(uni.site))

ms1 <- glm(histo_num ~ ns(age_b,df =4) + ns(cd4_b,df=4)+ as.factor(sex)  + as.factor(risk) + as.factor(rna_bb) + offset(log(follow)) + as.factor(SITE) ,family = poisson(link = "log"), data = la)
summ_spline  <- summary(ms1)
v_spline <- sandwich(ms1)
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Coefs", "SE", "Variables")

rr_age_20 <- rr(20,summ_spline)
rr_age_25 <- rr(25,summ_spline)
rr_age_30 <- rr(30,summ_spline)
rr_age_40 <- rr(40,summ_spline)
rr_age_45 <- rr(45,summ_spline)
rr_age_50 <- rr(50,summ_spline)
rr_age_55 <- rr(55,summ_spline)
rr_age_60 <- rr(60,summ_spline)

SE_20 <- spline_se_fun(20)
SE_25 <- spline_se_fun(25)
SE_30 <- spline_se_fun(30)
SE_40 <- spline_se_fun(40)
SE_45 <- spline_se_fun(45)
SE_50 <- spline_se_fun(50)
SE_55 <- spline_se_fun(55)
SE_60 <- spline_se_fun(60)

CI_20 <- CI(rr_age_20,SE_20)
CI_25 <- CI(rr_age_25,SE_25)
CI_30 <- CI(rr_age_30,SE_30)
CI_40 <- CI(rr_age_40,SE_40)
CI_45 <- CI(rr_age_45,SE_45)
CI_50 <- CI(rr_age_50,SE_50)
CI_55 <- CI(rr_age_55,SE_55)
CI_60 <- CI(rr_age_60,SE_60)

abs_20 <- abs_val(20,summ_spline)
abs_25 <- abs_val(25,summ_spline)
abs_30 <- abs_val(30,summ_spline)
abs_40 <- abs_val(40,summ_spline)
abs_45 <- abs_val(45,summ_spline)
abs_50 <- abs_val(50,summ_spline)
abs_55 <- abs_val(55,summ_spline)
abs_60 <- abs_val(60,summ_spline)

rr_cd4_100 <- rr_cd4(100,summ_spline)
rr_cd4_500 <- rr_cd4(500,summ_spline)
rr_cd4_1500 <- rr_cd4(1500,summ_spline)
rr_cd4_2000 <- rr_cd4(2000,summ_spline)
rr_cd4_2500 <- rr_cd4(2500,summ_spline)

SE_100 <- spline_se_fun_cd(100)
SE_500 <- spline_se_fun_cd(500)
SE_1500 <- spline_se_fun_cd(1500)
SE_2000 <- spline_se_fun_cd(2000)
SE_2500 <- spline_se_fun_cd(2500)


CI_100 <- CI(rr_cd4_100,SE_100)
CI_500 <- CI(rr_cd4_500,SE_500)
CI_1500 <- CI(rr_cd4_1500,SE_1500)
CI_2000 <- CI(rr_cd4_2000,SE_2000)
CI_2500 <- CI(rr_cd4_2500,SE_2500)


abs_100 <- abs_val_cd(100,summ_spline)
abs_500 <- abs_val_cd(500,summ_spline)
abs_1500 <- abs_val_cd(1500,summ_spline)
abs_2000 <- abs_val_cd(2000,summ_spline)
abs_2500 <- abs_val_cd(2500,summ_spline)


age_spl_multi <- data.frame("RR" = c(rr_age_20,rr_age_25,rr_age_30,rr_age_40,rr_age_45,rr_age_50,rr_age_55,rr_age_60),
                            "Lower CI"= c(CI_20[1],CI_25[1],CI_30[1],CI_40[1],CI_45[1],CI_50[1],CI_55[1],CI_60[1]),
                            "Upper CI"= c(CI_20[2],CI_25[2],CI_30[2],CI_40[2],CI_45[2],CI_50[2],CI_55[2],CI_60[2]),
                            "p-value" = c(2*(1-pnorm(abs(abs_20/SE_20))), 2*(1-pnorm(abs(abs_25/SE_25))), 2*(1-pnorm(abs(abs_30/SE_30))),
                                          2*(1-pnorm(abs(abs_40/SE_40))), 2*(1-pnorm(abs(abs_45/SE_45))), 2*(1-pnorm(abs(abs_50/SE_50))),
                                          2*(1-pnorm(abs(abs_55/SE_55))), 2*(1-pnorm(abs(abs_60/SE_60)))))
colnames(age_spl_multi) <- c("RR","Lower CI","Upper CI","p-value")

cd4_spl_multi <- data.frame("RR" = c(rr_cd4_100,rr_cd4_500,rr_cd4_1500,rr_cd4_2000,rr_cd4_2500),
                            "Lower CI"= c(CI_100[1],CI_500[1],CI_1500[1],CI_2000[1],CI_2500[1]),
                            "Upper CI"= c(CI_100[2],CI_500[2],CI_1500[2],CI_2000[2],CI_2500[2]),
                            "p-value" = c(2*(1-pnorm(abs(abs_100/SE_100))), 2*(1-pnorm(abs(abs_500/SE_500))), 2*(1-pnorm(abs(abs_1500/SE_1500))),
                                          2*(1-pnorm(abs(abs_2000/SE_2000))), 2*(1-pnorm(abs(abs_2500/SE_2500)))))
colnames(cd4_spl_multi) <- c("RR","Lower CI","Upper CI","p-value")


multi.tab1 <- uni.fit3(ms1)
multi.tab1 <- multi.tab1[-c(1,2,3,4,5,6,7,8),]
multi.tab1 <- rbind(age_spl_multi,cd4_spl_multi,multi.tab1)


tab1 <- cbind(uni.tab1,multi.tab1)
rownames(tab1) <- c("Age 20","Age 25","Age 30","Age 40","Age 45","Age 50","Age 55","Age 60","CD4 100","CD4 500","CD4 1500","CD4 2000","CD4 2500",
                    "Male ", "Hetero ","IDU ", "Other/Unknown ",
                    "RNA Detectable", "Chile ", "Honduras ","Mexico ", "Peru ")

save(tab1, file="tab1.Rdata")



#### TN Data anaslysis ####
tn <- data.frame(tn)
tn$follow_time_yrs <- with(tn,ifelse(follow_time_yrs ==0,0.002737851, tn$follow_time_yrs))
tn$follow_time_yrs <- abs(tn$follow_time_yrs)
dd <- datadist(tn)
options(datadist= 'dd')
dd$limits$sex[2] <- "Female"
dd$limits$mig_ende[2] <- "No"
dd$limits$rna_bb[2] <- "Undetectable"
dd$limits$cd4_bb[2] <- "<=350 cells/ mm3"

# Relevel
tn$sex <- relevel(as.factor(tn$sex), ref = "Female")
tn$risk1 <- relevel(tn$risk1, ref = "MSM")
tn$rna_bb <- relevel(tn$rna_bb, ref = "Undetectable")
tn$cd4_bb <- relevel(tn$cd4_bb, ref = "<=350 cells/ mm3")

#uni.age <- glm(histo_num ~ age_at_first_visit + offset(log(follow_time_yrs)),family = poisson(link = "log"), data = tn)
uni.gen <- glm(histo_num ~as.factor(sex) + offset(log(follow_time_yrs)),family = poisson(link = "log"), data = tn)
uni.risk <- glm(histo_num ~as.factor(risk1) + offset(log(follow_time_yrs)),family = poisson(link = "log"), data = tn)
uni.rna <- glm(histo_num ~as.factor(rna_bb) + offset(log(follow_time_yrs)),family = poisson(link = "log"), data = tn)
uni.age <- glm(histo_num ~ age_b + offset(log(follow_time_yrs)), tn, family= poisson(link="log"))
uni.cd4 <- glm(histo_num ~ cd4_b + offset(log(follow_time_yrs)),family = poisson(link = "log"), data = tn)
uni.tab2 <- rbind(uni.fit.conti(uni.age,10),uni.fit.conti(uni.cd4,500),uni.fit(uni.gen),uni.fit3(uni.risk),uni.fit(uni.rna))

ms1 <- glm(histo_num ~ age_b + cd4_b + as.factor(sex)  + as.factor(risk) + as.factor(rna_bb) + offset(log(follow_time_yrs)),family = poisson(link = "log"), data = tn)
#MS1 <- Glm(histo_num ~ age_b + cd4_b + sex + risk + rna_bb + offset(log(follow_time_yrs)),family = poisson(link = "log"), data = tn)
# summ_spline  <- summary(ms1)
# v_spline <- sandwich(ms1)
# se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
# summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
# summ_spline$vars <- rownames(summ_spline)
# colnames(summ_spline) <- c("Coefs", "SE", "Variables")
# 
# rr_age_20 <- rr(20,summ_spline)
# rr_age_25 <- rr(25,summ_spline)
# rr_age_30 <- rr(30,summ_spline)
# rr_age_40 <- rr(40,summ_spline)
# rr_age_45 <- rr(45,summ_spline)
# rr_age_50 <- rr(50,summ_spline)
# rr_age_55 <- rr(55,summ_spline)
# rr_age_60 <- rr(60,summ_spline)
# 
# SE_20 <- spline_se_fun(20)
# SE_25 <- spline_se_fun(25)
# SE_30 <- spline_se_fun(30)
# SE_40 <- spline_se_fun(40)
# SE_45 <- spline_se_fun(45)
# SE_50 <- spline_se_fun(50)
# SE_55 <- spline_se_fun(55)
# SE_60 <- spline_se_fun(60)
# 
# CI_20 <- CI(rr_age_20,SE_20)
# CI_25 <- CI(rr_age_25,SE_25)
# CI_30 <- CI(rr_age_30,SE_30)
# CI_40 <- CI(rr_age_40,SE_40)
# CI_45 <- CI(rr_age_45,SE_45)
# CI_50 <- CI(rr_age_50,SE_50)
# CI_55 <- CI(rr_age_55,SE_55)
# CI_60 <- CI(rr_age_60,SE_60)
# 
# abs_20 <- abs_val(20,summ_spline)
# abs_25 <- abs_val(25,summ_spline)
# abs_30 <- abs_val(30,summ_spline)
# abs_40 <- abs_val(40,summ_spline)
# abs_45 <- abs_val(45,summ_spline)
# abs_50 <- abs_val(50,summ_spline)
# abs_55 <- abs_val(55,summ_spline)
# abs_60 <- abs_val(60,summ_spline)
# 
# rr_cd4_100 <- rr_cd4_tn(100,summ_spline)
# rr_cd4_200 <- rr_cd4_tn(200,summ_spline)
# rr_cd4_300 <- rr_cd4_tn(300,summ_spline)
# rr_cd4_400 <- rr_cd4_tn(400,summ_spline)
# rr_cd4_600 <- rr_cd4_tn(600,summ_spline)
# rr_cd4_700 <- rr_cd4_tn(700,summ_spline)
# rr_cd4_800 <- rr_cd4_tn(800,summ_spline)
# rr_cd4_1000 <- rr_cd4_tn(1000,summ_spline)
# rr_cd4_1500 <- rr_cd4_tn(1500,summ_spline)
# 
# SE_100 <- spline_se_fun_cd_tn(100)
# SE_200 <- spline_se_fun_cd_tn(200)
# SE_300 <- spline_se_fun_cd_tn(300)
# SE_400 <- spline_se_fun_cd_tn(400)
# SE_600 <- spline_se_fun_cd_tn(600)
# SE_700 <- spline_se_fun_cd_tn(700)
# SE_800 <- spline_se_fun_cd_tn(800)
# SE_1000 <- spline_se_fun_cd_tn(1000)
# SE_1500 <- spline_se_fun_cd_tn(1500)
# 
# CI_100 <- CI(rr_cd4_100,SE_100)
# CI_200 <- CI(rr_cd4_200,SE_200)
# CI_300 <- CI(rr_cd4_300,SE_300)
# CI_400 <- CI(rr_cd4_400,SE_400)
# CI_600 <- CI(rr_cd4_600,SE_600)
# CI_700 <- CI(rr_cd4_700,SE_700)
# CI_800 <- CI(rr_cd4_800,SE_800)
# CI_1000 <- CI(rr_cd4_1000,SE_1000)
# CI_1500 <- CI(rr_cd4_1500,SE_1500)
# 
# abs_100 <- abs_val_cd_tn(100,summ_spline)
# abs_200 <- abs_val_cd_tn(200,summ_spline)
# abs_300 <- abs_val_cd_tn(300,summ_spline)
# abs_400 <- abs_val_cd_tn(400,summ_spline)
# abs_600 <- abs_val_cd_tn(600,summ_spline)
# abs_700 <- abs_val_cd_tn(700,summ_spline)
# abs_800 <- abs_val_cd_tn(800,summ_spline)
# abs_1000 <- abs_val_cd_tn(1000,summ_spline)
# abs_1500 <- abs_val_cd_tn(1500,summ_spline)
# 
# 
# cd4_spl_multi <- data.frame("RR" = c(rr_cd4_100,rr_cd4_200,rr_cd4_300,rr_cd4_400,rr_cd4_600,rr_cd4_700,rr_cd4_800,rr_cd4_1000,rr_cd4_1500),
#                             "Lower CI"= c(CI_100[1],CI_200[1],CI_300[1],CI_400[1],CI_600[1],CI_700[1],CI_800[1],CI_1000[1],CI_1500[1]),
#                             "Upper CI"= c(CI_100[2],CI_200[2],CI_300[2],CI_400[2],CI_600[2],CI_700[2],CI_800[2],CI_1000[2],CI_1500[2]),
#                             "p-value" = c(2*(1-pnorm(abs(abs_100/SE_100))), 2*(1-pnorm(abs(abs_200/SE_200))), 2*(1-pnorm(abs(abs_300/SE_300))),
#                                           2*(1-pnorm(abs(abs_400/SE_400))),  2*(1-pnorm(abs(abs_600/SE_600))),
#                                           2*(1-pnorm(abs(abs_700/SE_700))), 2*(1-pnorm(abs(abs_800/SE_800))),2*(1-pnorm(abs(abs_1000/SE_1000))),2*(1-pnorm(abs(abs_1500/SE_1500)))))
# colnames(cd4_spl_multi) <- c("RR","Lower CI","Upper CI","p-value")
# age_spl_multi <- data.frame("RR" = c(rr_age_20,rr_age_25,rr_age_30,rr_age_40,rr_age_45,rr_age_50,rr_age_55,rr_age_60),
#                             "Lower CI"= c(CI_20[1],CI_25[1],CI_30[1],CI_40[1],CI_45[1],CI_50[1],CI_55[1],CI_60[1]),
#                             "Upper CI"= c(CI_20[2],CI_25[2],CI_30[2],CI_40[2],CI_45[2],CI_50[2],CI_55[2],CI_60[2]),
#                             "p-value" = c(2*(1-pnorm(abs(abs_20/SE_20))), 2*(1-pnorm(abs(abs_25/SE_25))), 2*(1-pnorm(abs(abs_30/SE_30))),
#                                           2*(1-pnorm(abs(abs_40/SE_40))), 2*(1-pnorm(abs(abs_45/SE_45))), 2*(1-pnorm(abs(abs_50/SE_50))),
#                                           2*(1-pnorm(abs(abs_55/SE_55))), 2*(1-pnorm(abs(abs_60/SE_60)))))
# colnames(age_spl_multi) <- c("RR","Lower CI","Upper CI","p-value")
# 
# 
# 

multi.tab2 <- uni.fit3(ms1)
multi.tab2 <- multi.tab2[-c(1,2),]
multi.tab2 <- rbind(multi.fit.conti(ms1),multi.tab2)


tab2 <- cbind(uni.tab2,multi.tab2)
rownames(tab2) <- c("Age 35 vs 25","CD4 1000 vs. 500" ,"Males ", "Hetero","IDU", "Other/Unknown",
                    "RNA Detectable ")

save(tab2, file="tab2.Rdata")


#Time updated ART
la$art_sdt <- with(la,ifelse(art_sd == enrol_d, art_sd +1,art_sd))
la_t <- la %>% select(patient_id,cd4_b,rna_b,histo_d,art_sdt,age_b,histo,sex,risk1,SITE,last_d,baseline,enrol_d,cd4_base_date,rna_base_date,baseline_y,mig_ende,art_sd)
la_t$histo_d <- as.Date(la_t$histo_d,origin="1970-01-01")
la_t <- la_t %>% tidyr::pivot_longer(-c(patient_id,age_b,histo,sex,risk1,SITE,cd4_b,rna_b,baseline,cd4_base_date,rna_base_date,baseline_y,mig_ende,art_sd,last_d,histo_d), 
                                     names_to = c("lab", ".value"), 
                                     names_sep="_" )
la_t$date <- with(la_t,ifelse(lab == "art",sdt,d))
la_t$date <- as.Date(la_t$date,origin="1970-01-01")
la_t <- la_t %>% group_by(patient_id) %>%  distinct(lab,date, .keep_all =  TRUE)
la_t <- la_t %>% filter(!is.na(date))
la_t <- la_t %>% group_by(patient_id) %>% arrange(date)
la_t <- la_t %>% filter(sdt >= baseline | is.na(sdt))
la_t <- la_t %>% filter(sdt <= histo_d | is.na(sdt) | is.na(histo_d))
la_t <- la_t %>% filter(sdt <= last_d | is.na(sdt))
# la_t <- la_t %>% mutate(lag_t = difftime(date,lag(date),units = "days"))
# la_t$lag_t <- with(la_t,ifelse(lag_t == 0,1,lag_t))
# la_t$lag_t <- with(la_t,ifelse(is.na(lag_t),0,lag_t))
# la_t <- la_t %>% mutate(follow = cumsum(lag_t))
# la_t$follow_yrs <- (la_t$follow/365.25)

la_t <- la_t %>% group_by(patient_id) %>% mutate(lag_t = difftime(date,lag(date),units = "days"))
la_t <- la_t %>% group_by(patient_id) %>%  tidyr::fill(lag_t, .direction = "up")
la_t$end_d <- with(la_t,ifelse(!is.na(histo_d), histo_d,last_d))
la_t <- la_t %>% group_by(patient_id) %>% add_count()
la_t$follow_up <- with(la_t,ifelse(lab == "art" | n ==1, abs(difftime(as.Date(end_d,origin="1970-01-01"),as.Date(date,origin="1970-01-01"),units = "days")),lag_t))
la_t$follow_up <- with(la_t,ifelse(follow_up ==0,1,follow_up))
la_t$rna <- abs(la_t$rna_b)
la_t$sq_cd4 <- sqrt(la_t$cd4_b)
la_t$log_rna <- log10(la_t$rna_b)
la_t$follow_up <- (la_t$follow_up)/365.25
la_t$end_d <- as.Date(la_t$end_d, origin = "1970-01-01")

# ## filter Dates beyond last date and histo date
# la_t <- la_t %>% group_by(patient_id)  %>% slice(seq_len(min(which(lab == "last" | lab == "histo"), n())))


## Histo updated
la_t <- la_t %>% group_by(patient_id) %>%  mutate(first = row_number() == 1)
la_t$histo_num <- with(la_t,ifelse(!is.na(histo_d),1,0))
la_t$histo_num <- with(la_t,ifelse(n ==2 & first == "TRUE",0,histo_num))
la_t$histo <- with(la_t,ifelse(lab == "histo","Yes","No"))

# ## Time variable for time updated analysis
# la_t$time <- with(la_t,ifelse(lab != "enrol" & follow_yrs == 0,0.002737851,follow_yrs))
# 
# ## Offset follow up calculation
# la_t <- la_t  %>% 
#   group_by(patient_id) %>% 
#   mutate(follow_up = max(follow_yrs))
# la_t$follow_up <- with(la_t,ifelse(follow_up == 0,0.0027,follow_up))
# la_t$time <- as.numeric(la_t$time)

## Pre ART and Post ART
#la_t <- la_t %>% group_by(patient_id) %>% tidyr::fill(sdt, .direction = "up")
la_t$art_status <- with(la_t,ifelse(date<= art_sd | is.na(art_sd), "pre-ART initiation", "post-ART initiation"))
la_t$art_status <- with(la_t,ifelse(first == "FALSE","post-ART initiation",art_status ))

## pre-ART division only for individuals with no ART at baseline
la_t_pre <- la_t %>% filter(art_status == "pre-ART initiation")

## Post-ART division only for individuals with no ART at baseline
la_t_post <- la_t %>% filter(art_status == "post-ART initiation")


### Regression Models ###


# Relevel
la_t$sex <- relevel(as.factor(la_t$sex), ref = "Female")
la_t$SITE <- relevel(as.factor(la_t$SITE), ref = "Brazil")
la_t$risk1 <- relevel(la_t$risk1, ref = "MSM")
la_t$art_status <- relevel(as.factor(la_t$art_status), ref = "post-ART initiation")
la_t$mig_ende <- relevel(as.factor(la_t$mig_ende), ref = "No")
check <- la_t %>% dplyr::select(patient_id,histo_num,baseline,date,art_sd,end_d,follow_up)
## Non spline univariate models
uni0 <- glm(histo_num ~ offset(log(follow_up)),family = poisson(link = "log"),data=la_t)
uni.gen <- glm(histo_num ~as.factor(sex) + offset(log(follow_up)),family = poisson(link = "log"), data = la_t)
uni.risk <- glm(histo_num ~as.factor(risk1) + offset(log(follow_up)),family = poisson(link = "log"), data = la_t)
uni.site <- glm(histo_num ~as.factor(SITE) + offset(log(follow_up)),family = poisson(link = "log"), data = la_t)
uni.art <- glm(histo_num ~as.factor(art_status) + offset(log(follow_up)),family = poisson(link = "log"), data = la_t)
uni.mig <- glm(histo_num ~as.factor(mig_ende) + offset(log(follow_up)),family = poisson(link = "log"), data = la_t)

## Splines ##
ages <- round(min(la_t$age_b)):round(max(la_t$age_b))
sp <- as.data.frame(ns(ages, df=4))
sp <- cbind(ages, sp) 
sp <- as.data.frame(sp)
rownames(sp) <- sp$age
sp <- sp[, -1]
colnames(sp) <- c("X1", "X2", "X3", "X4")

uni.age <- glm(histo_num ~ ns(age_b,df =4) + offset(log(follow_up)),family = poisson(link =log), data = la_t)
summ_spline  <- summary(uni.age)
v_spline <- sandwich::vcovHC(uni.age,cluster = la_t$patient_id)
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Coefs", "SE", "Variables")

rr_age_20 <- rr(20,summ_spline)
rr_age_25 <- rr(25,summ_spline)
rr_age_30 <- rr(30,summ_spline)
rr_age_40 <- rr(40,summ_spline)
rr_age_45 <- rr(45,summ_spline)
rr_age_50 <- rr(50,summ_spline)
rr_age_55 <- rr(55,summ_spline)
rr_age_60 <- rr(60,summ_spline)


SE_20 <- spline_se_fun(20)
SE_25 <- spline_se_fun(25)
SE_30 <- spline_se_fun(30)
SE_40 <- spline_se_fun(40)
SE_45 <- spline_se_fun(45)
SE_50 <- spline_se_fun(50)
SE_55 <- spline_se_fun(55)
SE_60 <- spline_se_fun(60)

CI <- function(rr,SE) 
{ exp(log(rr) + c(-1,1)*1.96*SE)}


CI_20 <- CI(rr_age_20,SE_20)
CI_25 <- CI(rr_age_25,SE_25)
CI_30 <- CI(rr_age_30,SE_30)
CI_40 <- CI(rr_age_40,SE_40)
CI_45 <- CI(rr_age_45,SE_45)
CI_50 <- CI(rr_age_50,SE_50)
CI_55 <- CI(rr_age_55,SE_55)
CI_60 <- CI(rr_age_60,SE_60)

abs_20 <- abs_val(20,summ_spline)
abs_25 <- abs_val(25,summ_spline)
abs_30 <- abs_val(30,summ_spline)
abs_40 <- abs_val(40,summ_spline)
abs_45 <- abs_val(45,summ_spline)
abs_50 <- abs_val(50,summ_spline)
abs_55 <- abs_val(55,summ_spline)
abs_60 <- abs_val(60,summ_spline)

age_spline <- data.frame("RR" = c(rr_age_20,rr_age_25,rr_age_30,rr_age_40,rr_age_45,rr_age_50,rr_age_55,rr_age_60),
                         "Lower CI"= c(CI_20[1],CI_25[1],CI_30[1],CI_40[1],CI_45[1],CI_50[1],CI_55[1],CI_60[1]),
                         "Upper CI"= c(CI_20[2],CI_25[2],CI_30[2],CI_40[2],CI_45[2],CI_50[2],CI_55[2],CI_60[2]),
                         "p-value" = c(2*(1-pnorm(abs(abs_20/SE_20))), 2*(1-pnorm(abs(abs_25/SE_25))), 2*(1-pnorm(abs(abs_30/SE_30))),
                                       2*(1-pnorm(abs(abs_40/SE_40))), 2*(1-pnorm(abs(abs_45/SE_45))), 2*(1-pnorm(abs(abs_50/SE_50))),
                                       2*(1-pnorm(abs(abs_55/SE_55))), 2*(1-pnorm(abs(abs_60/SE_60)))))
colnames(age_spline) <- c("RR","Lower CI","Upper CI","p-value")


cd4s <- seq(min(la_t$sq_cd4,na.rm = T),max(la_t$sq_cd4,na.rm=T), by = 0.01)
sp_cd4 <- as.data.frame(ns(cd4s, df=4))
sp_cd4 <- cbind(cd4s, sp_cd4) 
sp_cd4 <- as.data.frame(sp_cd4)
rownames(sp_cd4) <- sp_cd4$cd4
sp_cd4 <- sp_cd4[, -1]
colnames(sp_cd4) <- c("X1", "X2", "X3", "X4")

uni.cd <- glm(histo_num ~ ns(sq_cd4,df =4) + offset(log(follow_up)),family = poisson(link =log), data = la_t)
summ_spline  <- summary(uni.cd)
v_spline <-sandwich::vcovHC(uni.cd,cluster = la_t$patient_id)
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Coefs", "SE", "Variables")

rr_sq_5 <- rr_sq(5,summ_spline)
rr_sq_10 <- rr_sq(10,summ_spline)
rr_sq_18.70 <- rr_sq(18.70,summ_spline)
rr_sq_22.40 <- rr_sq(22.40,summ_spline)


SE_5 <- spline_se_fun_sq_cd(5)
SE_10 <- spline_se_fun_sq_cd(10)
SE_18.70 <- spline_se_fun_sq_cd(18.70)
SE_22.40 <- spline_se_fun_sq_cd(22.40)



CI_5 <- CI(rr_sq_5,SE_5)
CI_10 <- CI(rr_sq_10,SE_10)
CI_18.70 <- CI(rr_sq_18.70,SE_18.70)
CI_22.40 <- CI(rr_sq_22.40,SE_22.40)



abs_5 <- abs_sq_cd(5,summ_spline)
abs_10 <- abs_sq_cd(10,summ_spline)
abs_18.70 <- abs_sq_cd(18.70,summ_spline)
abs_22.40 <- abs_sq_cd(22.40,summ_spline)



cd4_spline <- data.frame("RR" = c(rr_sq_5,rr_sq_10,rr_sq_18.70,rr_sq_22.40),
                         "Lower CI"= c(CI_5[1],CI_10[1],CI_18.70[1],CI_20[1]),
                         "Upper CI"= c(CI_5[2],CI_10[2],CI_18.70[2],CI_20[2]),
                         "p-value" = c(2*(1-pnorm(abs(abs_5/SE_5))), 2*(1-pnorm(abs(abs_10/SE_10))), 2*(1-pnorm(abs(abs_18.70/SE_18.70))),
                                       2*(1-pnorm(abs(abs_25/SE_22.40)))))
colnames(cd4_spline) <- c("RR","Lower CI","Upper CI","p-value")

yr <- round(min(la_t$baseline_y,na.rm = TRUE)):round(max(la_t$baseline_y,na.rm = TRUE))
sp_yr <- as.data.frame(ns(yr, df=4))
sp_yr <- cbind(yr, sp_yr) 
sp_yr <- as.data.frame(sp_yr)
rownames(sp_yr) <- sp_yr$yr
sp_yr <- sp_yr[, -1]
colnames(sp_yr) <- c("X1", "X2", "X3", "X4")

uni.yrs <- glm(histo_num ~ ns(baseline_y,df =4) + offset(log(follow_up)),family = poisson(link =log), data = la_t)
summ_spline  <- summary(uni.yrs)
v_spline <- sandwich::vcovHC(uni.yrs, cluster =la_t$patient_id)
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Coefs", "SE", "Variables")


rr_yr_2005 <- rr_yr(2005,summ_spline)
rr_yr_2000 <- rr_yr(2000,summ_spline)
rr_yr_2015 <- rr_yr(2015,summ_spline)
rr_yr_2020 <- rr_yr(2020,summ_spline)

SE_2005 <- spline_se_yr(2005)
SE_2000 <- spline_se_yr(2000)
SE_2015 <- spline_se_yr(2015)
SE_2020 <- spline_se_yr(2020)


CI_2005 <- CI(rr_yr_2005,SE_2005)
CI_2000 <- CI(rr_yr_2000,SE_2000)
CI_2015 <- CI(rr_yr_2015,SE_2015)
CI_2020 <- CI(rr_yr_2020,SE_2020)



abs_2005 <- abs_val_yrs(2005,summ_spline)
abs_2000 <- abs_val_yrs(2000,summ_spline)
abs_2015 <- abs_val_yrs(2015,summ_spline)
abs_2020 <- abs_val_yrs(2020,summ_spline)



yr_spline <- data.frame("RR" = c(rr_yr_2000,rr_yr_2005,rr_yr_2015,rr_yr_2020),
                        "Lower CI"= c(CI_2000[1],CI_2005[1],CI_2015[1],CI_2020[1]),
                        "Upper CI"= c(CI_2000[2],CI_2005[2],CI_2015[2],CI_2020[2]),
                        "p-value" = c(2*(1-pnorm(abs(abs_2000/SE_2005))), 2*(1-pnorm(abs(abs_2005/SE_2000))),
                                      2*(1-pnorm(abs(abs_2015/SE_2015))), 2*(1-pnorm(abs(abs_2020/SE_2020)))))
colnames(yr_spline) <- c("RR","Lower CI","Upper CI","p-value")


rna <- round(min(la_t$log_rna,na.rm = TRUE)):round(max(la_t$log_rna,na.rm = TRUE))
sp_rna <- as.data.frame(ns(rna, df=4))
sp_rna <- cbind(rna,sp_rna) 
sp_rna <- as.data.frame(sp_rna)
rownames(sp_rna) <- sp_rna$rna
sp_rna <- sp_rna[, -1]
colnames(sp_rna) <- c("X1", "X2", "X3", "X4")

uni.rna <- glm(histo_num ~ ns(log_rna,df =4) + offset(log(follow_up)),family = poisson(link =log), data = la_t)
summ_spline  <- summary(uni.rna)
v_spline <- sandwich::vcovHC(uni.rna,cluster=la_t$patient_id)
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Coefs", "SE", "Variables")


rr_rna_1 <- rr_rna(1,summ_spline)
rr_rna_3 <- rr_rna(3,summ_spline)
rr_rna_7 <- rr_rna(7,summ_spline)


SE_1 <- spline_se_rna(1)
SE_3 <- spline_se_rna(3)
SE_7 <- spline_se_rna(7)

CI_1 <- CI(rr_rna_1,SE_1)
CI_3 <- CI(rr_rna_3,SE_3)
CI_7 <- CI(rr_rna_7,SE_7)

abs_1 <- abs_val_rna(1,summ_spline)
abs_3 <- abs_val_rna(3,summ_spline)
abs_7 <- abs_val_rna(7,summ_spline)


rna_spline <- data.frame("RR" = c(rr_rna_1,rr_rna_3,rr_rna_7),
                         "Lower CI"= c(CI_1[1],CI_3[1],CI_7[1]),
                         "Upper CI"= c(CI_1[2],CI_3[2],CI_7[2]),
                         "p-value" = c(2*(1-pnorm(abs(abs_1/SE_1))), 2*(1-pnorm(abs(abs_3/SE_3))),
                                       2*(1-pnorm(abs(abs_7/SE_7)))))
colnames(rna_spline) <- c("RR","Lower CI","Upper CI","p-value")



uni.tab1 <- rbind(age_spline ,cd4_spline,rna_spline,yr_spline,uni.fit1(uni.gen,la_t$patient_id),
                  uni.fit2(uni.risk,la_t$patient_id),uni.fit2(uni.site,la_t$patient_id),uni.fit1(uni.art,la_t$patient_id),uni.fit1(uni.mig,la_t$patient_id))

ms1 <- glm(histo_num ~  ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up))
           + as.factor(risk1) + as.factor(SITE) + as.factor(art_status) + as.factor(mig_ende),family = poisson(link = "log"), data = la_t)

summ_spline  <- summary(ms1)
v_spline <- sandwich::vcovHC(ms1, cluster=la_t$patient_id)
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Coefs", "SE", "Variables")

rr_age_20 <- rr(20,summ_spline)
rr_age_25 <- rr(25,summ_spline)
rr_age_30 <- rr(30,summ_spline)
rr_age_40 <- rr(40,summ_spline)
rr_age_45 <- rr(45,summ_spline)
rr_age_50 <- rr(50,summ_spline)
rr_age_55 <- rr(55,summ_spline)
rr_age_60 <- rr(60,summ_spline)

SE_20 <- spline_se_fun(20)
SE_25 <- spline_se_fun(25)
SE_30 <- spline_se_fun(30)
SE_40 <- spline_se_fun(40)
SE_45 <- spline_se_fun(45)
SE_50 <- spline_se_fun(50)
SE_55 <- spline_se_fun(55)
SE_60 <- spline_se_fun(60)

CI_20 <- CI(rr_age_20,SE_20)
CI_25 <- CI(rr_age_25,SE_25)
CI_30 <- CI(rr_age_30,SE_30)
CI_40 <- CI(rr_age_40,SE_40)
CI_45 <- CI(rr_age_45,SE_45)
CI_50 <- CI(rr_age_50,SE_50)
CI_55 <- CI(rr_age_55,SE_55)
CI_60 <- CI(rr_age_60,SE_60)

abs_20 <- abs_val(20,summ_spline)
abs_25 <- abs_val(25,summ_spline)
abs_30 <- abs_val(30,summ_spline)
abs_40 <- abs_val(40,summ_spline)
abs_45 <- abs_val(45,summ_spline)
abs_50 <- abs_val(50,summ_spline)
abs_55 <- abs_val(55,summ_spline)
abs_60 <- abs_val(60,summ_spline)

rr_sq_5 <- rr_sq_sh(5,summ_spline)
rr_sq_10 <- rr_sq_sh(10,summ_spline)
rr_sq_18.70 <- rr_sq_sh(18.70,summ_spline)
rr_sq_22.40 <- rr_sq_sh(22.40,summ_spline)


SE_sq5 <- spline_se_sh(5)
SE_sq10 <- spline_se_sh(10)
SE_sq18.70 <- spline_se_sh(18.70)
SE_sq22.40 <- spline_se_sh(22.40)


CI_sq5 <- CI(rr_sq_5,SE_sq5)
CI_sq10 <- CI(rr_sq_10,SE_sq10)
CI_sq18.70 <- CI(rr_sq_18.70,SE_sq18.70)
CI_sq22.40 <- CI(rr_sq_22.40,SE_sq22.40)


abs_sq5 <- abs_sh(5,summ_spline)
abs_sq10 <- abs_sh(10,summ_spline)
abs_sq18.70 <- abs_sh(18.70,summ_spline)
abs_sq22.40 <- abs_sh(22.40,summ_spline)



rr_yr_2005 <- rr_yr(2005,summ_spline)
rr_yr_2000 <- rr_yr(2000,summ_spline)
rr_yr_2015 <- rr_yr(2015,summ_spline)
rr_yr_2020 <- rr_yr(2020,summ_spline)

SE_2005 <- spline_se_yr(2005)
SE_2000 <- spline_se_yr(2000)
SE_2015 <- spline_se_yr(2015)
SE_2020 <- spline_se_yr(2020)


CI_2005 <- CI(rr_yr_2005,SE_2005)
CI_2000 <- CI(rr_yr_2000,SE_2000)
CI_2015 <- CI(rr_yr_2015,SE_2015)
CI_2020 <- CI(rr_yr_2020,SE_2020)



abs_2005 <- abs_val_yrs(2005,summ_spline)
abs_2000 <- abs_val_yrs(2000,summ_spline)
abs_2015 <- abs_val_yrs(2015,summ_spline)
abs_2020 <- abs_val_yrs(2020,summ_spline)

rr_rna_1 <- rr_rna(1,summ_spline)
rr_rna_3 <- rr_rna(3,summ_spline)
rr_rna_7 <- rr_rna(7,summ_spline)


SE_1 <- spline_se_rna(1)
SE_3 <- spline_se_rna(3)
SE_7 <- spline_se_rna(7)



CI_1 <- CI(rr_rna_1,SE_1)
CI_3 <- CI(rr_rna_3,SE_3)
CI_7 <- CI(rr_rna_7,SE_7)




abs_1 <- abs_val_rna(1,summ_spline)
abs_3 <- abs_val_rna(3,summ_spline)
abs_7 <- abs_val_rna(7,summ_spline)

age_spl_multi <- data.frame("RR" = c(rr_age_20,rr_age_25,rr_age_30,rr_age_40,rr_age_45,rr_age_50,rr_age_55,rr_age_60),
                            "Lower CI"= c(CI_20[1],CI_25[1],CI_30[1],CI_40[1],CI_45[1],CI_50[1],CI_55[1],CI_60[1]),
                            "Upper CI"= c(CI_20[2],CI_25[2],CI_30[2],CI_40[2],CI_45[2],CI_50[2],CI_55[2],CI_60[2]),
                            "p-value" = c(2*(1-pnorm(abs(abs_20/SE_20))), 2*(1-pnorm(abs(abs_25/SE_25))), 2*(1-pnorm(abs(abs_30/SE_30))),
                                          2*(1-pnorm(abs(abs_40/SE_40))), 2*(1-pnorm(abs(abs_45/SE_45))), 2*(1-pnorm(abs(abs_50/SE_50))),
                                          2*(1-pnorm(abs(abs_55/SE_55))), 2*(1-pnorm(abs(abs_60/SE_60)))))
colnames(age_spl_multi) <- c("RR","Lower CI","Upper CI","p-value")


cd4_spline_multi <- data.frame("RR" = c(rr_sq_5,rr_sq_10,rr_sq_18.70,rr_sq_22.40),
                               "Lower CI"= c(CI_sq5[1],CI_sq10[1],CI_sq18.70[1],CI_sq22.40[1]),
                               "Upper CI"= c(CI_sq5[2],CI_sq10[2],CI_sq18.70[2],CI_sq22.40[2]),
                               "p-value" = c(2*(1-pnorm(abs(abs_sq5/SE_sq5))), 2*(1-pnorm(abs(abs_sq10/SE_sq10))), 2*(1-pnorm(abs(abs_sq18.70/SE_sq18.70))),
                                             2*(1-pnorm(abs(abs_sq22.40/SE_sq22.40)))))

colnames(cd4_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")


yr_spline_multi <- data.frame("RR" = c(rr_yr_2000,rr_yr_2005,rr_yr_2015,rr_yr_2020),
                              "Lower CI"= c(CI_2000[1],CI_2005[1],CI_2015[1],CI_2020[1]),
                              "Upper CI"= c(CI_2000[2],CI_2005[2],CI_2015[2],CI_2020[2]),
                              "p-value" = c(2*(1-pnorm(abs(abs_2000/SE_2005))), 2*(1-pnorm(abs(abs_2005/SE_2000))),
                                            2*(1-pnorm(abs(abs_2015/SE_2015))), 2*(1-pnorm(abs(abs_2020/SE_2020)))))
colnames(yr_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")

rna_spline_multi <- data.frame("RR" = c(rr_rna_1,rr_rna_3,rr_rna_7),
                               "Lower CI"= c(CI_1[1],CI_3[1],CI_7[1]),
                               "Upper CI"= c(CI_1[2],CI_3[2],CI_7[2]),
                               "p-value" = c(2*(1-pnorm(abs(abs_1/SE_1))), 2*(1-pnorm(abs(abs_3/SE_3))),
                                             2*(1-pnorm(abs(abs_7/SE_7)))))
colnames(rna_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")

multi.tab1 <- uni.fit2(ms1)
multi.tab1 <- multi.tab1[-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14),]
multi.tab1 <- rbind(age_spl_multi,cd4_spline_multi,rna_spline_multi,yr_spline_multi,multi.tab1)


tab3 <- cbind(uni.tab1,multi.tab1)
rownames(tab3) <- c("Age 20","Age 25","Age 30","Age 40","Age 45","Age 50","Age 55","Age 60",
                    "CD4 25","CD4 100","CD4 350","CD4 500",
                    "RNA 1","RNA 3","RNA 7", "Year 2000","Year 2005"," Year 2015","Year 2020",
                    "Male ", "Hetero ","IDU ", "Other/Unknown ",  "Chile ", "Honduras ","Mexico ", 
                    "Peru ","pre-ART initiation", "Migration status")

save(tab3, file="tab3.Rdata")

incid_new <- incid_long %>% select(patient_id,value,year,follow_time_yrs)
ld <- full_join(la_t,incid_new, by = "patient_id") %>% select(patient_id,value,year,follow_time_yrs,art_sd,last_d,baseline,age_b,sq_cd4,log_rna,baseline_y,sex,follow_up,risk1,SITE,art_status,mig_ende,
                                                              histo,histo_d,histo_num)
#ld <- ld %>% distinct("patient_id", .keep_all = T)
ld$art_y <- year(ld$art_sd)
#ld <- ld %>% filter(!is.na(art_y))
ld$year_end_d <- as.Date(paste(ld$year,12,31, sep="-"))
ld$year_beg_d <- as.Date(paste(ld$year,01,01, sep="-"))
ld <- ld %>% group_by(patient_id) %>% mutate(d1 = pmin(baseline,art_sd, na.rm=T))
ld <- ld %>% group_by(patient_id) %>% mutate(d2 = pmax(baseline,art_sd, na.rm=T))
ld$t1 <- difftime(ld$d2,ld$d1,units = "days")/365.25
ld$t2 <- difftime(ld$year_end_d,ld$d2,units = "days")/365.25
ld$indi <- with(ld,ifelse(year == art_y & !is.na(art_y),1,0))

ld <- ld %>% group_by(patient_id) %>% dplyr::mutate(indi_row = cumsum(indi),
                                                    value1 = ifelse(indi ==1 & indi_row ==1 ,t1,
                                                                    ifelse(indi == 1 & indi_row == 2,t2,NA)))

ld$value1 <- with(ld,ifelse(is.na(value1),value,value1))

ld <- ld %>% distinct(patient_id,value1,year, .keep_all = T)

ld <- ld %>% group_by(patient_id) %>% arrange(year,value1)


## redefining art_atatus
ld$indi <- with(ld,ifelse(art_status == "post-ART initiation",1,0))
ld <- ld %>% group_by(patient_id) %>% dplyr::mutate(indi_row = cumsum(indi))
ld$art_status <- with(ld,ifelse(indi_row >= 1,"post-ART initiation","pre-ART initiation"))
ld$value1 <- with(ld,ifelse(value1 == 0,0.002737851,value1))
ms1 <- glm(histo_num ~  ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(year,df=4)+ as.factor(sex)  + offset(log(value1))
           + as.factor(risk1) + as.factor(SITE) + as.factor(art_status) + as.factor(mig_ende),family = poisson(link = "log"), data = ld)

summ_spline  <- summary(ms1)
v_spline <- sandwich::vcovHC(ms1, cluster=la_t$patient_id)
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Coefs", "SE", "Variables")

rr_age_20 <- rr(20,summ_spline)
rr_age_25 <- rr(25,summ_spline)
rr_age_30 <- rr(30,summ_spline)
rr_age_40 <- rr(40,summ_spline)
rr_age_45 <- rr(45,summ_spline)
rr_age_50 <- rr(50,summ_spline)
rr_age_55 <- rr(55,summ_spline)
rr_age_60 <- rr(60,summ_spline)

SE_20 <- spline_se_fun(20)
SE_25 <- spline_se_fun(25)
SE_30 <- spline_se_fun(30)
SE_40 <- spline_se_fun(40)
SE_45 <- spline_se_fun(45)
SE_50 <- spline_se_fun(50)
SE_55 <- spline_se_fun(55)
SE_60 <- spline_se_fun(60)

CI_20 <- CI(rr_age_20,SE_20)
CI_25 <- CI(rr_age_25,SE_25)
CI_30 <- CI(rr_age_30,SE_30)
CI_40 <- CI(rr_age_40,SE_40)
CI_45 <- CI(rr_age_45,SE_45)
CI_50 <- CI(rr_age_50,SE_50)
CI_55 <- CI(rr_age_55,SE_55)
CI_60 <- CI(rr_age_60,SE_60)

abs_20 <- abs_val(20,summ_spline)
abs_25 <- abs_val(25,summ_spline)
abs_30 <- abs_val(30,summ_spline)
abs_40 <- abs_val(40,summ_spline)
abs_45 <- abs_val(45,summ_spline)
abs_50 <- abs_val(50,summ_spline)
abs_55 <- abs_val(55,summ_spline)
abs_60 <- abs_val(60,summ_spline)

rr_sq_5 <- rr_sq_sh(5,summ_spline)
rr_sq_10 <- rr_sq_sh(10,summ_spline)
rr_sq_18.70 <- rr_sq_sh(18.70,summ_spline)
rr_sq_22.40 <- rr_sq_sh(22.40,summ_spline)


SE_sq5 <- spline_se_sh(5)
SE_sq10 <- spline_se_sh(10)
SE_sq18.70 <- spline_se_sh(18.70)
SE_sq22.40 <- spline_se_sh(22.40)


CI_sq5 <- CI(rr_sq_5,SE_sq5)
CI_sq10 <- CI(rr_sq_10,SE_sq10)
CI_sq18.70 <- CI(rr_sq_18.70,SE_sq18.70)
CI_sq22.40 <- CI(rr_sq_22.40,SE_sq22.40)


abs_sq5 <- abs_sh(5,summ_spline)
abs_sq10 <- abs_sh(10,summ_spline)
abs_sq18.70 <- abs_sh(18.70,summ_spline)
abs_sq22.40 <- abs_sh(22.40,summ_spline)



rr_yr_2005 <- rr_yr2(2005,summ_spline)
rr_yr_2000 <- rr_yr2(2000,summ_spline)
rr_yr_2015 <- rr_yr2(2015,summ_spline)
rr_yr_2020 <- rr_yr2(2020,summ_spline)

SE_2005 <- spline_se_yr2(2005)
SE_2000 <- spline_se_yr2(2000)
SE_2015 <- spline_se_yr2(2015)
SE_2020 <- spline_se_yr2(2020)


CI_2005 <- CI(rr_yr_2005,SE_2005)
CI_2000 <- CI(rr_yr_2000,SE_2000)
CI_2015 <- CI(rr_yr_2015,SE_2015)
CI_2020 <- CI(rr_yr_2020,SE_2020)



abs_2005 <- abs_val_yrs2(2005,summ_spline)
abs_2000 <- abs_val_yrs2(2000,summ_spline)
abs_2015 <- abs_val_yrs2(2015,summ_spline)
abs_2020 <- abs_val_yrs2(2020,summ_spline)

rr_rna_1 <- rr_rna(1,summ_spline)
rr_rna_3 <- rr_rna(3,summ_spline)
rr_rna_7 <- rr_rna(7,summ_spline)


SE_1 <- spline_se_rna(1)
SE_3 <- spline_se_rna(3)
SE_7 <- spline_se_rna(7)



CI_1 <- CI(rr_rna_1,SE_1)
CI_3 <- CI(rr_rna_3,SE_3)
CI_7 <- CI(rr_rna_7,SE_7)




abs_1 <- abs_val_rna(1,summ_spline)
abs_3 <- abs_val_rna(3,summ_spline)
abs_7 <- abs_val_rna(7,summ_spline)

age_spl_multi <- data.frame("RR" = c(rr_age_20,rr_age_25,rr_age_30,rr_age_40,rr_age_45,rr_age_50,rr_age_55,rr_age_60),
                            "Lower CI"= c(CI_20[1],CI_25[1],CI_30[1],CI_40[1],CI_45[1],CI_50[1],CI_55[1],CI_60[1]),
                            "Upper CI"= c(CI_20[2],CI_25[2],CI_30[2],CI_40[2],CI_45[2],CI_50[2],CI_55[2],CI_60[2]),
                            "p-value" = c(2*(1-pnorm(abs(abs_20/SE_20))), 2*(1-pnorm(abs(abs_25/SE_25))), 2*(1-pnorm(abs(abs_30/SE_30))),
                                          2*(1-pnorm(abs(abs_40/SE_40))), 2*(1-pnorm(abs(abs_45/SE_45))), 2*(1-pnorm(abs(abs_50/SE_50))),
                                          2*(1-pnorm(abs(abs_55/SE_55))), 2*(1-pnorm(abs(abs_60/SE_60)))))
colnames(age_spl_multi) <- c("RR","Lower CI","Upper CI","p-value")


cd4_spline_multi <- data.frame("RR" = c(rr_sq_5,rr_sq_10,rr_sq_18.70,rr_sq_22.40),
                               "Lower CI"= c(CI_sq5[1],CI_sq10[1],CI_sq18.70[1],CI_sq22.40[1]),
                               "Upper CI"= c(CI_sq5[2],CI_sq10[2],CI_sq18.70[2],CI_sq22.40[2]),
                               "p-value" = c(2*(1-pnorm(abs(abs_sq5/SE_sq5))), 2*(1-pnorm(abs(abs_sq10/SE_sq10))), 2*(1-pnorm(abs(abs_sq18.70/SE_sq18.70))),
                                             2*(1-pnorm(abs(abs_sq22.40/SE_sq22.40)))))

colnames(cd4_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")


yr_spline_multi <- data.frame("RR" = c(rr_yr_2000,rr_yr_2005,rr_yr_2015,rr_yr_2020),
                              "Lower CI"= c(CI_2000[1],CI_2005[1],CI_2015[1],CI_2020[1]),
                              "Upper CI"= c(CI_2000[2],CI_2005[2],CI_2015[2],CI_2020[2]),
                              "p-value" = c(2*(1-pnorm(abs(abs_2000/SE_2005))), 2*(1-pnorm(abs(abs_2005/SE_2000))),
                                            2*(1-pnorm(abs(abs_2015/SE_2015))), 2*(1-pnorm(abs(abs_2020/SE_2020)))))
colnames(yr_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")

rna_spline_multi <- data.frame("RR" = c(rr_rna_1,rr_rna_3,rr_rna_7),
                               "Lower CI"= c(CI_1[1],CI_3[1],CI_7[1]),
                               "Upper CI"= c(CI_1[2],CI_3[2],CI_7[2]),
                               "p-value" = c(2*(1-pnorm(abs(abs_1/SE_1))), 2*(1-pnorm(abs(abs_3/SE_3))),
                                             2*(1-pnorm(abs(abs_7/SE_7)))))
colnames(rna_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")

multi.tab1 <- uni.fit2(ms1)
multi.tab1 <- multi.tab1[-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14),]
tab3_extra <- rbind(age_spl_multi,cd4_spline_multi,rna_spline_multi,yr_spline_multi,multi.tab1)
rownames(tab3_extra) <- c("Age 20","Age 25","Age 30","Age 40","Age 45","Age 50","Age 55","Age 60",
                          "CD4 25","CD4 100","CD4 350","CD4 500",
                          "RNA 1","RNA 3","RNA 7", "Year 2000","Year 2005"," Year 2015","Year 2020",
                          "Male ", "Hetero ","IDU ", "Other/Unknown ",  "Chile ", "Honduras ","Mexico ", 
                          "Peru ","pre-ART initiation", "Migration status")


save(tab3_extra, file="tab3_extra.Rdata")
### TRT updated TN analysis ###
#Time updated ART
tn$age_b <- tn$age_at_first_visit
tn$art_sd <- tn$enrol_d + (tn$AGE_AT_ART_START - tn$age_at_first_visit)
tn$art_sdt <- with(tn,ifelse(art_sd == enrol_d, art_sd +1,art_sd))
tn_t <- tn %>% dplyr::select(CFAR_PID,birthsex,mig_ende,risk1,baseline,art_sdt,cd4_b,rna_b,last_d,histo_d,enrol_d,age_b,histo,baseline_y,art_sd)

tn_t <- tn_t %>% tidyr::pivot_longer(-c(CFAR_PID,age_b,histo,birthsex,risk1,cd4_b,rna_b,baseline,mig_ende,baseline_y,histo_d,last_d,art_sd), 
                                     names_to = c("lab", ".value"), 
                                     names_sep="_" )
tn_t$date <- with(tn_t,ifelse(lab == "art",sdt,d))
tn_t$date <- as.Date(tn_t$date,origin="1970-01-01")
tn_t <- tn_t %>% group_by(CFAR_PID) %>%  distinct(lab,date, .keep_all =  TRUE)
tn_t <- tn_t %>% filter(!is.na(date))
tn_t <- tn_t %>% group_by(CFAR_PID) %>% arrange(date)

tn_t <- tn_t %>% filter(sdt >= baseline | is.na(sdt))
tn_t <- tn_t %>% filter(sdt <= histo_d | is.na(sdt) | is.na(histo_d))
tn_t$last_d <- as.Date(tn_t$last_d,origin="1970-01-01")
tn_t <- tn_t %>% filter(sdt <=  last_d | is.na(sdt))

tn_t <- tn_t %>% group_by(CFAR_PID) %>% mutate(lag_t = difftime(date,lag(date),units = "days"))
tn_t <- tn_t %>% group_by(CFAR_PID) %>%  tidyr::fill(lag_t, .direction = "up")
tn_t$end_d <- with(tn_t,ifelse(!is.na(histo_d), histo_d,last_d))
tn_t <- tn_t %>% group_by(CFAR_PID) %>% add_count()
tn_t$follow_up <- with(tn_t,ifelse(lab == "art" | n ==1, abs(difftime(as.Date(end_d,origin="1970-01-01"),as.Date(date,origin="1970-01-01"),units = "days")),lag_t))
tn_t$follow_up <- with(tn_t,ifelse(follow_up ==0,1,follow_up))
tn_t$rna <- abs(tn_t$rna_b)
tn_t$sq_cd4 <- sqrt(tn_t$cd4_b)
tn_t$log_rna <- log10(tn_t$rna_b)

## Histo updated
tn_t <- tn_t %>% group_by(CFAR_PID) %>%  mutate(first = row_number() == 1)
tn_t$histo_num <- with(tn_t,ifelse(!is.na(histo_d),1,0))
tn_t$histo_num <- with(tn_t,ifelse(n ==2 & first == "TRUE",0,histo_num))
tn_t$histo <- with(tn_t,ifelse(lab == "histo","Yes","No"))


## Pre ART and Post ART
#tn_t <- tn_t %>% group_by(patient_id) %>% tidyr::fill(sdt, .direction = "up")
tn_t$art_status <- with(tn_t,ifelse(date< art_sd | is.na(art_sd), "pre-ART initiation", "post-ART initiation"))
tn_t$art_status <- with(tn_t,ifelse(first == "FALSE","post-ART initiation",art_status ))

## pre-ART division only for individuals with no ART at baseline
tn_t_pre <- tn_t %>% filter(art_status == "pre-ART initiation")

## Post-ART division only for individuals with no ART at baseline
#tn_t_post <- tn_t %>% filter(art_status == "post-ART initiation")
tn_t_post1 <- tn_t %>% filter(sdt > baseline)
tn_t_post <- tn_t %>% filter(art_status == "post-ART initiation")


### Regression Models ###

## removing tansgender cate. bacause of error 


# Relevel
tn_t$birthsex <- relevel(as.factor(tn_t$birthsex), ref = "Female")
tn_t$risk1 <- relevel(tn_t$risk1, ref = "MSM")
tn_t$art_status <- relevel(as.factor(tn_t$art_status), ref = "post-ART initiation")
tn_t$mig_ende <- relevel(as.factor(tn_t$mig_ende), ref = "No")

## Non spline univariate models
uni.gen <- glm(histo_num ~as.factor(birthsex) + offset(log(follow_up)),family = poisson(link = "log"), data = tn_t)
uni.risk <- glm(histo_num ~as.factor(risk1) + offset(log(follow_up)),family = poisson(link = "log"), data = tn_t)
uni.art <- glm(histo_num ~as.factor(art_status) + offset(log(follow_up)),family = poisson(link = "log"), data = tn_t)
uni.mig <- glm(histo_num ~as.factor(mig_ende) + offset(log(follow_up)),family = poisson(link = "log"), data = tn_t)
#uni0 <- glm(histo_num ~offset(log(follow_up)),family = poisson(link = "log"), data = tn_t)
## Splines ##
# ages <- round(min(tn_t$age_b)):round(max(tn_t$age_b))
# sp <- as.data.frame(ns(ages, df=4))
# sp <- cbind(ages, sp) 
# sp <- as.data.frame(sp)
# rownames(sp) <- sp$age
# sp <- sp[, -1]
# colnames(sp) <- c("X1", "X2", "X3", "X4")

uni.age <- glm(histo_num ~ age_b + offset(log(follow_up)),family = poisson(link =log), data = tn_t)
# summ_spline  <- summary(uni.age)
# v_spline <- sandwich::vcovHC(uni.age, cluster=tn_t$CFAR_PID)
# se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
# summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
# summ_spline$vars <- rownames(summ_spline)
# colnames(summ_spline) <- c("Coefs", "SE", "Variables")
# 
# rr_age_20 <- rr(20,summ_spline)
# rr_age_25 <- rr(25,summ_spline)
# rr_age_30 <- rr(30,summ_spline)
# rr_age_40 <- rr(40,summ_spline)
# rr_age_45 <- rr(45,summ_spline)
# rr_age_50 <- rr(50,summ_spline)
# rr_age_55 <- rr(55,summ_spline)
# rr_age_60 <- rr(60,summ_spline)
# 
# 
# SE_20 <- spline_se_fun(20)
# SE_25 <- spline_se_fun(25)
# SE_30 <- spline_se_fun(30)
# SE_40 <- spline_se_fun(40)
# SE_45 <- spline_se_fun(45)
# SE_50 <- spline_se_fun(50)
# SE_55 <- spline_se_fun(55)
# SE_60 <- spline_se_fun(60)
# 
# CI <- function(rr,SE) 
# { exp(log(rr) + c(-1,1)*1.96*SE)}
# 
# 
# CI_20 <- CI(rr_age_20,SE_20)
# CI_25 <- CI(rr_age_25,SE_25)
# CI_30 <- CI(rr_age_30,SE_30)
# CI_40 <- CI(rr_age_40,SE_40)
# CI_45 <- CI(rr_age_45,SE_45)
# CI_50 <- CI(rr_age_50,SE_50)
# CI_55 <- CI(rr_age_55,SE_55)
# CI_60 <- CI(rr_age_60,SE_60)
# 
# abs_20 <- abs_val(20,summ_spline)
# abs_25 <- abs_val(25,summ_spline)
# abs_30 <- abs_val(30,summ_spline)
# abs_40 <- abs_val(40,summ_spline)
# abs_45 <- abs_val(45,summ_spline)
# abs_50 <- abs_val(50,summ_spline)
# abs_55 <- abs_val(55,summ_spline)
# abs_60 <- abs_val(60,summ_spline)
# 
# age_spline <- data.frame("RR" = c(rr_age_20,rr_age_25,rr_age_30,rr_age_40,rr_age_45,rr_age_50,rr_age_55,rr_age_60),
#                          "Lower CI"= c(CI_20[1],CI_25[1],CI_30[1],CI_40[1],CI_45[1],CI_50[1],CI_55[1],CI_60[1]),
#                          "Upper CI"= c(CI_20[2],CI_25[2],CI_30[2],CI_40[2],CI_45[2],CI_50[2],CI_55[2],CI_60[2]),
#                          "p-value" = c(2*(1-pnorm(abs(abs_20/SE_20))), 2*(1-pnorm(abs(abs_25/SE_25))), 2*(1-pnorm(abs(abs_30/SE_30))),
#                                        2*(1-pnorm(abs(abs_40/SE_40))), 2*(1-pnorm(abs(abs_45/SE_45))), 2*(1-pnorm(abs(abs_50/SE_50))),
#                                        2*(1-pnorm(abs(abs_55/SE_55))), 2*(1-pnorm(abs(abs_60/SE_60)))))
# colnames(age_spline) <- c("RR","Lower CI","Upper CI","p-value")
# 
# 
# 
# cd4s <- seq(min(tn_t$sq_cd4,na.rm = T),max(tn_t$sq_cd4,na.rm=T), by = 0.01)
# sp_cd4 <- as.data.frame(ns(cd4s, df=4))
# sp_cd4 <- cbind(cd4s, sp_cd4) 
# sp_cd4 <- as.data.frame(sp_cd4)
# rownames(sp_cd4) <- sp_cd4$cd4
# sp_cd4 <- sp_cd4[, -1]
# colnames(sp_cd4) <- c("X1", "X2", "X3", "X4")
# 
# 

uni.cd <- glm(histo_num ~ sq_cd4 + offset(log(follow_up)),family = poisson(link =log), data = tn_t)
# summ_spline  <- summary(uni.cd)
# v_spline <- sandwich::vcovHC(uni.cd, cluster=tn_t$CFAR_PID)
# se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
# summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
# summ_spline$vars <- rownames(summ_spline)
# colnames(summ_spline) <- c("Coefs", "SE", "Variables")
# 
# rr_sq_5 <- rr_sq(5,summ_spline)
# rr_sq_10 <- rr_sq(10,summ_spline)
# rr_sq_18.70 <- rr_sq(18.70,summ_spline)
# rr_sq_22.40 <- rr_sq(22.40,summ_spline)
# 
# 
# SE_5 <- spline_se_fun_sq_cd(5)
# SE_10 <- spline_se_fun_sq_cd(10)
# SE_18.70 <- spline_se_fun_sq_cd(18.70)
# SE_22.40 <- spline_se_fun_sq_cd(22.40)
# 
# 
# 
# CI_5 <- CI(rr_sq_5,SE_5)
# CI_10 <- CI(rr_sq_10,SE_10)
# CI_18.70 <- CI(rr_sq_18.70,SE_18.70)
# CI_22.40 <- CI(rr_sq_22.40,SE_22.40)
# 
# 
# abs_5 <- abs_sq_cd(5,summ_spline)
# abs_10 <- abs_sq_cd(10,summ_spline)
# abs_18.70 <- abs_sq_cd(18.70,summ_spline)
# abs_22.40 <- abs_sq_cd(22.40,summ_spline)
# 
# 
# cd4_spline <- data.frame("RR" = c(rr_sq_5,rr_sq_10,rr_sq_18.70,rr_sq_22.40),
#                          "Lower CI"= c(CI_5[1],CI_10[1],CI_18.70[1],CI_22.40[1]),
#                          "Upper CI"= c(CI_5[2],CI_10[2],CI_18.70[2],CI_22.40[2]),
#                          "p-value" = c(2*(1-pnorm(abs(abs_5/SE_5))), 2*(1-pnorm(abs(abs_10/SE_10))), 2*(1-pnorm(abs(abs_18.70/SE_18.70))),
#                                        2*(1-pnorm(abs(abs_22.40/SE_22.40)))))
# colnames(cd4_spline) <- c("RR","Lower CI","Upper CI","p-value")
# 
# yr <- round(min(tn_t$baseline_y,na.rm = TRUE)):round(max(tn_t$baseline_y,na.rm = TRUE))
# sp_yr <- as.data.frame(ns(yr, df=4))
# sp_yr <- cbind(yr, sp_yr) 
# sp_yr <- as.data.frame(sp_yr)
# rownames(sp_yr) <- sp_yr$yr
# sp_yr <- sp_yr[, -1]
# colnames(sp_yr) <- c("X1", "X2", "X3", "X4")

uni.yrs <- glm(histo_num ~ baseline_y + offset(log(follow_up)),family = poisson(link =log), data = tn_t)
# summ_spline  <- summary(uni.yrs)
# v_spline <- sandwich::vcovHC(uni.yrs, cluster=tn_t$CFAR_PID)
# se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
# summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
# summ_spline$vars <- rownames(summ_spline)
# colnames(summ_spline) <- c("Coefs", "SE", "Variables")
# 
# 
# rr_yr_2005 <- rr_yr(2005,summ_spline)
# rr_yr_2000 <- rr_yr(2000,summ_spline)
# rr_yr_2015 <- rr_yr(2015,summ_spline)
# rr_yr_2020 <- rr_yr(2020,summ_spline)
# 
# SE_2005 <- spline_se_yr(2005)
# SE_2000 <- spline_se_yr(2000)
# SE_2015 <- spline_se_yr(2015)
# SE_2020 <- spline_se_yr(2020)
# 
# 
# CI_2005 <- CI(rr_yr_2005,SE_2005)
# CI_2000 <- CI(rr_yr_2000,SE_2000)
# CI_2015 <- CI(rr_yr_2015,SE_2015)
# CI_2020 <- CI(rr_yr_2020,SE_2020)
# 
# 
# 
# abs_2005 <- abs_val_yrs(2005,summ_spline)
# abs_2000 <- abs_val_yrs(2000,summ_spline)
# abs_2015 <- abs_val_yrs(2015,summ_spline)
# abs_2020 <- abs_val_yrs(2020,summ_spline)
# 
# 
# 
# yr_spline <- data.frame("RR" = c(rr_yr_2005,rr_yr_2000,rr_yr_2015,rr_yr_2020),
#                         "Lower CI"= c(CI_2005[1],CI_2000[1],CI_2015[1],CI_2020[1]),
#                         "Upper CI"= c(CI_2005[2],CI_2000[2],CI_2015[2],CI_2020[2]),
#                         "p-value" = c(2*(1-pnorm(abs(abs_2005/SE_2005))), 2*(1-pnorm(abs(abs_2000/SE_2000))),
#                                       2*(1-pnorm(abs(abs_2015/SE_2015))), 2*(1-pnorm(abs(abs_2020/SE_2020)))))
# colnames(yr_spline) <- c("RR","Lower CI","Upper CI","p-value")
# 
# 
# 
# rna <- round(min(tn_t$log_rna,na.rm = TRUE)):round(max(tn_t$log_rna,na.rm = TRUE))
# sp_rna <- as.data.frame(ns(rna, df=4))
# sp_rna <- cbind(rna,sp_rna) 
# sp_rna <- as.data.frame(sp_rna)
# rownames(sp_rna) <- sp_rna$rna
# sp_rna <- sp_rna[, -1]
# colnames(sp_rna) <- c("X1", "X2", "X3", "X4")

uni.rna <- glm(histo_num ~ log_rna + offset(log(follow_up)),family = poisson(link =log), data = tn_t)
# summ_spline  <- summary(uni.rna)
# v_spline <- sandwich::vcovHC(uni.rna, cluster=tn_t$CFAR_PID)
# se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
# summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
# summ_spline$vars <- rownames(summ_spline)
# colnames(summ_spline) <- c("Coefs", "SE", "Variables")
# 
# 
# rr_rna_1 <- rr_rna(1,summ_spline)
# rr_rna_3 <- rr_rna(3,summ_spline)
# rr_rna_7 <- rr_rna(7,summ_spline)
# 
# 
# SE_1 <- spline_se_rna(1)
# SE_3 <- spline_se_rna(3)
# SE_7 <- spline_se_rna(7)
# 
# 
# 
# CI_1 <- CI(rr_rna_1,SE_1)
# CI_3 <- CI(rr_rna_3,SE_3)
# CI_7 <- CI(rr_rna_7,SE_7)
# 
# 
# 
# 
# abs_1 <- abs_val_rna(1,summ_spline)
# abs_3 <- abs_val_rna(3,summ_spline)
# abs_7 <- abs_val_rna(7,summ_spline)
# 
# 
# 
# 
# rna_spline <- data.frame("RR" = c(rr_rna_1,rr_rna_3,rr_rna_7),
#                          "Lower CI"= c(CI_1[1],CI_3[1],CI_7[1]),
#                          "Upper CI"= c(CI_1[2],CI_3[2],CI_7[2]),
#                          "p-value" = c(2*(1-pnorm(abs(abs_1/SE_1))), 2*(1-pnorm(abs(abs_3/SE_3))),
#                                        2*(1-pnorm(abs(abs_7/SE_7)))))
# colnames(rna_spline) <- c("RR","Lower CI","Upper CI","p-value")
# 
# 

uni.tab1 <- rbind(uni.fit.conti1(uni.age,tn_t$CFAR_PID,10),uni.fit.conti1(uni.cd,tn_t$CFAR_PID,50),uni.fit.conti1(uni.rna,tn_t$CFAR_PID,4),uni.fit1(uni.yrs,tn_t$CFAR_PID),
                  uni.fit1(uni.gen,tn_t$CFAR_PID),uni.fit2(uni.risk,tn_t$CFAR_PID),uni.fit1(uni.art,tn_t$CFAR_PID),uni.fit1(uni.mig,tn_t$CFAR_PID))

ms1 <- glm(histo_num ~  age_b + sq_cd4 + log_rna + baseline_y + as.factor(birthsex)  + offset(log(follow_up)) + as.factor(risk1)  + as.factor(art_status) + as.factor(mig_ende),family = poisson(link = "log"), data = tn_t)
# summ_spline  <- summary(ms1)
# v_spline <- sandwich::vcovHC(ms1, cluster=tn_t$CFAR_PID)
# se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
# summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
# summ_spline$vars <- rownames(summ_spline)
# colnames(summ_spline) <- c("Coefs", "SE", "Variables")
# 
# rr_age_20 <- rr(20,summ_spline)
# rr_age_25 <- rr(25,summ_spline)
# rr_age_30 <- rr(30,summ_spline)
# rr_age_40 <- rr(40,summ_spline)
# rr_age_45 <- rr(45,summ_spline)
# rr_age_50 <- rr(50,summ_spline)
# rr_age_55 <- rr(55,summ_spline)
# rr_age_60 <- rr(60,summ_spline)
# 
# SE_20 <- spline_se_fun(20)
# SE_25 <- spline_se_fun(25)
# SE_30 <- spline_se_fun(30)
# SE_40 <- spline_se_fun(40)
# SE_45 <- spline_se_fun(45)
# SE_50 <- spline_se_fun(50)
# SE_55 <- spline_se_fun(55)
# SE_60 <- spline_se_fun(60)
# 
# CI_20 <- CI(rr_age_20,SE_20)
# CI_25 <- CI(rr_age_25,SE_25)
# CI_30 <- CI(rr_age_30,SE_30)
# CI_40 <- CI(rr_age_40,SE_40)
# CI_45 <- CI(rr_age_45,SE_45)
# CI_50 <- CI(rr_age_50,SE_50)
# CI_55 <- CI(rr_age_55,SE_55)
# CI_60 <- CI(rr_age_60,SE_60)
# 
# abs_20 <- abs_val(20,summ_spline)
# abs_25 <- abs_val(25,summ_spline)
# abs_30 <- abs_val(30,summ_spline)
# abs_40 <- abs_val(40,summ_spline)
# abs_45 <- abs_val(45,summ_spline)
# abs_50 <- abs_val(50,summ_spline)
# abs_55 <- abs_val(55,summ_spline)
# abs_60 <- abs_val(60,summ_spline)
# 
# rr_sq_5 <- rr_sq_sh(5,summ_spline)
# rr_sq_10 <- rr_sq_sh(10,summ_spline)
# rr_sq_18.70 <- rr_sq_sh(18.70,summ_spline)
# rr_sq_22.40 <- rr_sq_sh(22.40,summ_spline)
# 
# 
# SE_sq5 <- spline_se_sh(5)
# SE_sq10 <- spline_se_sh(10)
# SE_sq18.70 <- spline_se_sh(18.70)
# SE_sq22.40 <- spline_se_sh(22.40)
# 
# CI_sq5 <- CI(rr_sq_5,SE_5)
# CI_sq10 <- CI(rr_sq_10,SE_10)
# CI_sq18.70 <- CI(rr_sq_18.70,SE_18.70)
# CI_sq22.40 <- CI(rr_sq_22.40,SE_22.40)
# 
# abs_sq5 <- abs_sh(5,summ_spline)
# abs_sq10 <- abs_sh(10,summ_spline)
# abs_sq18.70 <- abs_sh(18.70,summ_spline)
# abs_sq22.40 <- abs_sh(22.40,summ_spline)
# 
# 
# rr_yr_2005 <- rr_yr(2005,summ_spline)
# rr_yr_2000 <- rr_yr(2000,summ_spline)
# rr_yr_2015 <- rr_yr(2015,summ_spline)
# rr_yr_2020 <- rr_yr(2020,summ_spline)
# 
# SE_2005 <- spline_se_yr(2005)
# SE_2000 <- spline_se_yr(2000)
# SE_2015 <- spline_se_yr(2015)
# SE_2020 <- spline_se_yr(2020)
# 
# 
# CI_2005 <- CI(rr_yr_2005,SE_2005)
# CI_2000 <- CI(rr_yr_2000,SE_2000)
# CI_2015 <- CI(rr_yr_2015,SE_2015)
# CI_2020 <- CI(rr_yr_2020,SE_2020)
# 
# 
# 
# abs_2005 <- abs_val_yrs(2005,summ_spline)
# abs_2000 <- abs_val_yrs(2000,summ_spline)
# abs_2015 <- abs_val_yrs(2015,summ_spline)
# abs_2020 <- abs_val_yrs(2020,summ_spline)
# 
# rr_rna_1 <- rr_rna(1,summ_spline)
# rr_rna_3 <- rr_rna(3,summ_spline)
# rr_rna_7 <- rr_rna(7,summ_spline)
# 
# 
# SE_1 <- spline_se_rna(1)
# SE_3 <- spline_se_rna(3)
# SE_7 <- spline_se_rna(7)
# 
# 
# 
# CI_1 <- CI(rr_rna_1,SE_1)
# CI_3 <- CI(rr_rna_3,SE_3)
# CI_7 <- CI(rr_rna_7,SE_7)
# 
# 
# 
# 
# abs_1 <- abs_val_rna(1,summ_spline)
# abs_3 <- abs_val_rna(3,summ_spline)
# abs_7 <- abs_val_rna(7,summ_spline)
# 
# age_spl_multi <- data.frame("RR" = c(rr_age_20,rr_age_25,rr_age_30,rr_age_40,rr_age_45,rr_age_50,rr_age_55,rr_age_60),
#                             "Lower CI"= c(CI_20[1],CI_25[1],CI_30[1],CI_40[1],CI_45[1],CI_50[1],CI_55[1],CI_60[1]),
#                             "Upper CI"= c(CI_20[2],CI_25[2],CI_30[2],CI_40[2],CI_45[2],CI_50[2],CI_55[2],CI_60[2]),
#                             "p-value" = c(2*(1-pnorm(abs(abs_20/SE_20))), 2*(1-pnorm(abs(abs_25/SE_25))), 2*(1-pnorm(abs(abs_30/SE_30))),
#                                           2*(1-pnorm(abs(abs_40/SE_40))), 2*(1-pnorm(abs(abs_45/SE_45))), 2*(1-pnorm(abs(abs_50/SE_50))),
#                                           2*(1-pnorm(abs(abs_55/SE_55))), 2*(1-pnorm(abs(abs_60/SE_60)))))
# colnames(age_spl_multi) <- c("RR","Lower CI","Upper CI","p-value")
# 
# 
# cd4_spline_multi <- data.frame("RR" = c(rr_sq_5,rr_sq_10,rr_sq_18.70,rr_sq_22.40),
#                                "Lower CI"= c(CI_sq5[1],CI_sq10[1],CI_sq18.70[1],CI_sq22.40[1]),
#                                "Upper CI"= c(CI_sq5[2],CI_sq10[2],CI_sq18.70[2],CI_sq22.40[2]),
#                                "p-value" = c(2*(1-pnorm(abs(abs_sq5/SE_sq5))), 2*(1-pnorm(abs(abs_sq10/SE_sq10))), 2*(1-pnorm(abs(abs_sq18.70/SE_sq18.70))),
#                                              2*(1-pnorm(abs(abs_sq22.40/SE_sq22.40)))))
# colnames(cd4_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")
# 
# 
# yr_spline_multi <- data.frame("RR" = c(rr_yr_2005,rr_yr_2000,rr_yr_2015,rr_yr_2020),
#                               "Lower CI"= c(CI_2005[1],CI_2000[1],CI_2015[1],CI_2020[1]),
#                               "Upper CI"= c(CI_2005[2],CI_2000[2],CI_2015[2],CI_2020[2]),
#                               "p-value" = c(2*(1-pnorm(abs(abs_2005/SE_2005))), 2*(1-pnorm(abs(abs_2000/SE_2000))),
#                                             2*(1-pnorm(abs(abs_2015/SE_2015))), 2*(1-pnorm(abs(abs_2020/SE_2020)))))
# colnames(yr_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")
# 
# rna_spline_multi <- data.frame("RR" = c(rr_rna_1,rr_rna_3,rr_rna_7),
#                                "Lower CI"= c(CI_1[1],CI_3[1],CI_7[1]),
#                                "Upper CI"= c(CI_1[2],CI_3[2],CI_7[2]),
#                                "p-value" = c(2*(1-pnorm(abs(abs_1/SE_1))), 2*(1-pnorm(abs(abs_3/SE_3))),
#                                              2*(1-pnorm(abs(abs_7/SE_7)))))
# colnames(rna_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")
# 

multi.tab1 <- uni.fit2(ms1,tn_t$CFAR_PID)
multi.tab1 <- multi.tab1[-c(1,2,3),]
multi.tab1 <- rbind(multi.fit.conti1(ms1,tn_t$CFAR_PID),multi.tab1)


tab4 <- cbind(uni.tab1,multi.tab1)
rownames(tab4) <- c("Age 35 vs 25","CD4 100 vs. 50","RNA 5 vs 1", "Year",
                    "Male ", "Hetero ","IDU ", "Other/Unknown ", "pre-ART initiation","Migration status")

save(tab4, file="tab4.Rdata")

## Imputation ##

# set.seed(7)
# la_imp <- aregImpute(~ age_b + sq_cd4+ log_rna + baseline_y+ as.factor(sex)  + as.factor(risk1)  + as.factor(art_status) + as.factor(mig_ende) + follow_up ,data = la_t, n.impute=20,x=TRUE, nk=0)
# 
# imp.la <- fit.mult.impute(histo_num ~  ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE) + as.factor(art_status) + as.factor(mig_ende), data=la_t, 
#                           fitter =glm, fit.reps=TRUE,xtrans=la_imp,n.impute=20,fitargs=list(x=TRUE, y=TRUE,family= poisson(link="log")))
# 
# 
# imp.la.age <- fit.mult.impute(histo_num ~  ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE) + as.factor(art_status) + as.factor(mig_ende), data=la_t,
#                               fitter = glm, fit.reps=TRUE,xtrans=la_imp,n.impute=20,fitargs=list(x=TRUE, y=TRUE,family= poisson(link="log")))


imp <- la_t %>% select(patient_id,age_b,sq_cd4,log_rna, baseline_y,sex,risk1,art_status,mig_ende,follow_up,histo_num,sex,SITE)
imp <- mice(imp,m=5,printFlag = FALSE, idvars = "patient_id")

## Running GLM on the imputed data
imp.la <- with(imp,
               glm(histo_num ~  ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex) 
                   + offset(log(follow_up)) + as.factor(risk1)  + as.factor(art_status) + as.factor(mig_ende) + as.factor(SITE),
                   family = poisson(link = "log")))

imp.la.age <- with(imp,
                   glm(histo_num ~  + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex) 
                       + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE) + as.factor(art_status) + as.factor(mig_ende),
                       family = poisson(link = "log")))

a_age <-D1(imp.la, imp.la.age)
p_age <- round(a_age$result[4],digits=3)

imp.la.cd <- with(imp,
                  glm(histo_num ~   ns(age_b,df =4) + ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE) + as.factor(art_status) + as.factor(mig_ende), 
                      family = poisson(link = "log")))

a_cd  <- D1(imp.la, imp.la.cd)
p_cd <- round(a_cd$result[4],digits=3)

imp.la.rna <- with(imp,
                   glm(histo_num ~ ns(age_b,df =4) + ns(sq_cd4,df=2) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE) + as.factor(art_status) + as.factor(mig_ende), 
                       family = poisson(link = "log")))

a_rna <- D1(imp.la, imp.la.rna)
p_rna <- round(a_rna$result[4],digits=3)

imp.la.yr <- with(imp,
                  glm(histo_num ~    ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE) + as.factor(art_status) + as.factor(mig_ende),
                      family = poisson(link = "log")))

a_yr <- D1(imp.la, imp.la.yr)
p_yr <- round(a_yr$result[4],digits=3)


imp.la.gen <- with(imp,
                   glm(histo_num ~   ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE) + as.factor(art_status) + as.factor(mig_ende),
                       family = poisson(link = "log")))
a_gen <- D1(imp.la, imp.la.gen)
p_gen <- round(a_gen$result[4],digits=3)

imp.la.risk <- with(imp,
                    glm(histo_num ~ ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up))  + as.factor(SITE) + as.factor(art_status) + as.factor(mig_ende),
                        family = poisson(link = "log")))

a_risk <- D1(imp.la, imp.la.risk)
p_risk <- round(a_risk$result[4],digits=4)

imp.la.site <- with(imp,
                    glm(histo_num ~ ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1)  + as.factor(art_status) + as.factor(mig_ende), 
                        family = poisson(link = "log")))

a_site <- D1(imp.la,imp.la.site)
p_site <- round(a_site$result[4],digits=3)


imp.la.art <- with(imp,
                   glm(histo_num ~   ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE)  + as.factor(mig_ende), 
                       family = poisson(link = "log")))


a_art <- D1(imp.la, imp.la.art)
p_art <- round(a_art$result[4],digits=3)


imp.la.mig <- with(imp,
                   glm(histo_num ~   ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE) + as.factor(art_status), 
                       family = poisson(link = "log")))


a_mig <- D1(imp.la, imp.la.mig)
p_mig <- round(a_mig$result[4],digits=3)

## Variance covariance
t2 <- pool(imp.la)
m <- t2$m

## Applying rubin's rule
vcov_list <- list()
for (i in 1:m) {
  vcov_list[[i]] <- vcovHC(imp.la$analyses[[i]], cluster = imp$data$patient_id)
} 

vw <- Reduce('+',vcov_list)/m
#vcov_array <- array(unlist(vcov_list), dim = c(dim(vcov_list[[1]]),length(vcov_list)))
## Computing within imputaion variance 
#var_vcov <- apply(vcov_array,c(1,2),var)
#vw <- var_vcov
## Between imputaion variance covariance
qbar <- getqbar(t2)
qhats <- sapply(imp.la$analyses, coef)
vb <- (1 / (m-1)) * (qhats - qbar) %*% t(qhats - qbar)
vt <- vw + (1 + 1 / (m)) * vb 
vt <- vt[-c(1), -c(1)]


impu.la <- function(model,var) {
  s <-  summary(pool(model))
  s <-  s$estimate[-1]
  sand <- sqrt(diag(var))
  r = matrix(0,nrow = length(s), ncol = 4)
  for (i in 1:(length(s))){
    r[i,] <- c(exp(s[i]), exp(s[i] - 1.96*sand[i]),exp(s[i] + 1.96*sand[i]),round((1 - pnorm(abs(s[i]/sand[i]))) * 2,digits=5))
  }
  
  r <- as.data.frame(r)
  colnames(r) <- c( "RR","Lower CI", "Upper CI","p-value") 
  rownam <- rownames(s)
  rownames(r) <- rownam
  r$ind <- seq_len(nrow(r))
  r <- rbind(r,c(NA,NA,NA,p_gen,ind=14.1),c(NA,NA,NA,p_risk,ind=15.1),c(NA,NA,NA,p_art,ind=18.1),c(NA,NA,NA,p_mig,ind=19.1),c(NA,NA,NA,p_site,ind=20.1))
  r <- r %>% arrange(ind) %>% dplyr::select(-ind)
  return(r)
}


impu.la <- impu.la(imp.la,vt)
cd4s <- seq(min(la_t$sq_cd4,na.rm = T),max(la_t$sq_cd4,na.rm=T), by = 0.01)
sp_cd4 <- as.data.frame(ns(cd4s, df=4))
sp_cd4 <- cbind(cd4s, sp_cd4) 
sp_cd4 <- as.data.frame(sp_cd4)
rownames(sp_cd4) <- sp_cd4$cd4
sp_cd4 <- sp_cd4[, -1]
colnames(sp_cd4) <- c("X1", "X2", "X3", "X4")

ages <- round(min(la_t$age_b)):round(max(la_t$age_b))
sp <- as.data.frame(ns(ages, df=4))
sp <- cbind(ages, sp) 
sp <- as.data.frame(sp)
rownames(sp) <- sp$age
sp <- sp[, -1]
colnames(sp) <- c("X1", "X2", "X3", "X4")

yr <- round(min(la_t$baseline_y,na.rm = TRUE)):round(max(la_t$baseline_y,na.rm = TRUE))
sp_yr <- as.data.frame(ns(yr, df=4))
sp_yr <- cbind(yr, sp_yr) 
sp_yr <- as.data.frame(sp_yr)
rownames(sp_yr) <- sp_yr$yr
sp_yr <- sp_yr[, -1]
colnames(sp_yr) <- c("X1", "X2", "X3", "X4")

rna <- round(min(la_t$log_rna,na.rm = TRUE)):round(max(la_t$log_rna,na.rm = TRUE))
sp_rna <- as.data.frame(ns(rna, df=4))
sp_rna <- cbind(rna,sp_rna) 
sp_rna <- as.data.frame(sp_rna)
rownames(sp_rna) <- sp_rna$rna
sp_rna <- sp_rna[, -1]
colnames(sp_rna) <- c("X1", "X2", "X3", "X4")

summ_spline  <- summary(pool(imp.la))
v_spline <- vt
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- summ_spline[-1,c(1,2)]
summ_spline <- as.data.frame(cbind(summ_spline, se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Variables","Coefs", "SE" )

rr_age_20 <- rr(20,summ_spline)
rr_age_25 <- rr(25,summ_spline)
rr_age_30 <- rr(30,summ_spline)
rr_age_40 <- rr(40,summ_spline)
rr_age_45 <- rr(45,summ_spline)
rr_age_50 <- rr(50,summ_spline)
rr_age_55 <- rr(55,summ_spline)
rr_age_60 <- rr(60,summ_spline)

SE_20 <- spline_se_fun(20)
SE_25 <- spline_se_fun(25)
SE_30 <- spline_se_fun(30)
SE_40 <- spline_se_fun(40)
SE_45 <- spline_se_fun(45)
SE_50 <- spline_se_fun(50)
SE_55 <- spline_se_fun(55)
SE_60 <- spline_se_fun(60)

CI_20 <- CI(rr_age_20,SE_20)
CI_25 <- CI(rr_age_25,SE_25)
CI_30 <- CI(rr_age_30,SE_30)
CI_40 <- CI(rr_age_40,SE_40)
CI_45 <- CI(rr_age_45,SE_45)
CI_50 <- CI(rr_age_50,SE_50)
CI_55 <- CI(rr_age_55,SE_55)
CI_60 <- CI(rr_age_60,SE_60)

abs_20 <- abs_val(20,summ_spline)
abs_25 <- abs_val(25,summ_spline)
abs_30 <- abs_val(30,summ_spline)
abs_40 <- abs_val(40,summ_spline)
abs_45 <- abs_val(45,summ_spline)
abs_50 <- abs_val(50,summ_spline)
abs_55 <- abs_val(55,summ_spline)
abs_60 <- abs_val(60,summ_spline)

rr_sq_5 <- rr_sq_sh(5,summ_spline)
rr_sq_10 <- rr_sq_sh(10,summ_spline)
rr_sq_18.70 <- rr_sq_sh(18.70,summ_spline)
rr_sq_22.40 <- rr_sq_sh(22.40,summ_spline)


SE_sq5 <- spline_se_sh(5)
SE_sq10 <- spline_se_sh(10)
SE_sq18.70 <- spline_se_sh(18.70)
SE_sq22.40 <- spline_se_sh(22.40)

CI_sq5 <- CI(rr_sq_5,SE_5)
CI_sq10 <- CI(rr_sq_10,SE_10)
CI_sq18.70 <- CI(rr_sq_18.70,SE_18.70)
CI_sq22.40 <- CI(rr_sq_22.40,SE_22.40)

abs_sq5 <- abs_sh(5,summ_spline)
abs_sq10 <- abs_sh(10,summ_spline)
abs_sq18.70 <- abs_sh(18.70,summ_spline)
abs_sq22.40 <- abs_sh(22.40,summ_spline)


rr_yr_2005 <- rr_yr(2005,summ_spline)
rr_yr_2000 <- rr_yr(2000,summ_spline)
rr_yr_2015 <- rr_yr(2015,summ_spline)
rr_yr_2020 <- rr_yr(2020,summ_spline)

SE_2005 <- spline_se_yr(2005)
SE_2000 <- spline_se_yr(2000)
SE_2015 <- spline_se_yr(2015)
SE_2020 <- spline_se_yr(2020)


CI_2005 <- CI(rr_yr_2005,SE_2005)
CI_2000 <- CI(rr_yr_2000,SE_2000)
CI_2015 <- CI(rr_yr_2015,SE_2015)
CI_2020 <- CI(rr_yr_2020,SE_2020)



abs_2005 <- abs_val_yrs(2005,summ_spline)
abs_2000 <- abs_val_yrs(2000,summ_spline)
abs_2015 <- abs_val_yrs(2015,summ_spline)
abs_2020 <- abs_val_yrs(2020,summ_spline)

rr_rna_1 <- rr_rna(1,summ_spline)
rr_rna_3 <- rr_rna(3,summ_spline)
rr_rna_7 <- rr_rna(7,summ_spline)


SE_1 <- spline_se_rna(1)
SE_3 <- spline_se_rna(3)
SE_7 <- spline_se_rna(7)



CI_1 <- CI(rr_rna_1,SE_1)
CI_3 <- CI(rr_rna_3,SE_3)
CI_7 <- CI(rr_rna_7,SE_7)




abs_1 <- abs_val_rna(1,summ_spline)
abs_3 <- abs_val_rna(3,summ_spline)
abs_7 <- abs_val_rna(7,summ_spline)

age_spl_multi <- data.frame("RR" = c(rr_age_20,rr_age_25,rr_age_30,rr_age_40,rr_age_45,rr_age_50,rr_age_55,rr_age_60),
                            "Lower CI"= c(CI_20[1],CI_25[1],CI_30[1],CI_40[1],CI_45[1],CI_50[1],CI_55[1],CI_60[1]),
                            "Upper CI"= c(CI_20[2],CI_25[2],CI_30[2],CI_40[2],CI_45[2],CI_50[2],CI_55[2],CI_60[2]),
                            "p-value" = c(2*(1-pnorm(abs(abs_20/SE_20))), 2*(1-pnorm(abs(abs_25/SE_25))), 2*(1-pnorm(abs(abs_30/SE_30))),
                                          2*(1-pnorm(abs(abs_40/SE_40))), 2*(1-pnorm(abs(abs_45/SE_45))), 2*(1-pnorm(abs(abs_50/SE_50))),
                                          2*(1-pnorm(abs(abs_55/SE_55))), 2*(1-pnorm(abs(abs_60/SE_60)))))
colnames(age_spl_multi) <- c("RR","Lower CI","Upper CI","p-value")


cd4_spline_multi <- data.frame("RR" = c(rr_sq_5,rr_sq_10,rr_sq_18.70,rr_sq_22.40),
                               "Lower CI"= c(CI_sq5[1],CI_sq10[1],CI_sq18.70[1],CI_sq22.40[1]),
                               "Upper CI"= c(CI_sq5[2],CI_sq10[2],CI_sq18.70[2],CI_sq22.40[2]),
                               "p-value" = c(2*(1-pnorm(abs(abs_sq5/SE_sq5))), 2*(1-pnorm(abs(abs_sq10/SE_sq10))), 2*(1-pnorm(abs(abs_sq18.70/SE_sq18.70))),
                                             2*(1-pnorm(abs(abs_sq22.40/SE_sq22.40)))))
colnames(cd4_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")


yr_spline_multi <- data.frame("RR" = c(rr_yr_2005,rr_yr_2000,rr_yr_2015,rr_yr_2020),
                              "Lower CI"= c(CI_2005[1],CI_2000[1],CI_2015[1],CI_2020[1]),
                              "Upper CI"= c(CI_2005[2],CI_2000[2],CI_2015[2],CI_2020[2]),
                              "p-value" = c(2*(1-pnorm(abs(abs_2005/SE_2005))), 2*(1-pnorm(abs(abs_2000/SE_2000))),
                                            2*(1-pnorm(abs(abs_2015/SE_2015))), 2*(1-pnorm(abs(abs_2020/SE_2020)))))
colnames(yr_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")

rna_spline_multi <- data.frame("RR" = c(rr_rna_1,rr_rna_3,rr_rna_7),
                               "Lower CI"= c(CI_1[1],CI_3[1],CI_7[1]),
                               "Upper CI"= c(CI_1[2],CI_3[2],CI_7[2]),
                               "p-value" = c(2*(1-pnorm(abs(abs_1/SE_1))), 2*(1-pnorm(abs(abs_3/SE_3))),
                                             2*(1-pnorm(abs(abs_7/SE_7)))))
colnames(rna_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")

impu.la <- impu.la[-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14),]
tab5 <- rbind(c(NA,NA,NA,p_age),age_spl_multi,c(NA,NA,NA,p_cd),cd4_spline_multi,c(NA,NA,NA,p_rna),rna_spline_multi,c(NA,NA,NA,p_yr),yr_spline_multi,impu.la)
rownames(tab5) <- c("Age (ref=35)","20","25","30","40","45","50","55","60",
                    "CD4 (ref=200)","25 ","100","350","500","log HIV RNA(ref=5)"," 1",  "3","7"
                    ,"Year (ref=2010)", "2000","2005","2015","2020","Sex","Male ", "Risk factor (ref =MSM)","Hetero ","IDU ", "Other/Unknown ",
                    "ART initiation status (ref=post)","pre-ART initiation", "Migration status (ref=No)","Yes",
                    "Country (ref=Brazil)", "Chile ", "Honduras ","Mexico ", "Peru ")

save(tab5, file="tab5.Rdata")

## Extra imputed model for CCASAnet

imp <- ld %>% select(patient_id,age_b,sq_cd4,log_rna, baseline_y,sex,risk1,art_status,mig_ende,follow_up,histo_num,sex,SITE,year,value1)
imp <- mice(imp,m=5,printFlag = FALSE, idvars = "patient_id")

## Running GLM on the impu2ted data
imp.la <- with(imp,
               glm(histo_num ~  ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(year,df=4)+ as.factor(sex) 
                   + offset(log(value1)) + as.factor(risk1)  + as.factor(art_status) + as.factor(mig_ende) + as.factor(SITE),
                   family = poisson(link = "log")))

imp.la.age <- with(imp,
                   glm(histo_num ~  + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(year,df=4)+ as.factor(sex) 
                       + offset(log(value1)) + as.factor(risk1) + as.factor(SITE) + as.factor(art_status) + as.factor(mig_ende),
                       family = poisson(link = "log")))

a_age <-D1(imp.la, imp.la.age)
p_age <- round(a_age$result[4],digits=3)

imp.la.cd <- with(imp,
                  glm(histo_num ~   ns(age_b,df =4) + ns(log_rna,df =4) + ns(year,df=4)+ as.factor(sex)  + offset(log(value1)) + as.factor(risk1) + as.factor(SITE) + as.factor(art_status) + as.factor(mig_ende), 
                      family = poisson(link = "log")))

a_cd  <- D1(imp.la, imp.la.cd)
p_cd <- round(a_cd$result[4],digits=3)

imp.la.rna <- with(imp,
                   glm(histo_num ~ ns(age_b,df =4) + ns(sq_cd4,df=2) + ns(year,df=4)+ as.factor(sex)  + offset(log(value1)) + as.factor(risk1) + as.factor(SITE) + as.factor(art_status) + as.factor(mig_ende), 
                       family = poisson(link = "log")))

a_rna <- D1(imp.la, imp.la.rna)
p_rna <- round(a_rna$result[4],digits=3)

imp.la.yr <- with(imp,
                  glm(histo_num ~    ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + as.factor(sex)  + offset(log(value1)) + as.factor(risk1) + as.factor(SITE) + as.factor(art_status) + as.factor(mig_ende),
                      family = poisson(link = "log")))

a_yr <- D1(imp.la, imp.la.yr)
p_yr <- round(a_yr$result[4],digits=3)


imp.la.gen <- with(imp,
                   glm(histo_num ~   ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(year,df=4)  + offset(log(value1)) + as.factor(risk1) + as.factor(SITE) + as.factor(art_status) + as.factor(mig_ende),
                       family = poisson(link = "log")))
a_gen <- D1(imp.la, imp.la.gen)
p_gen <- round(a_gen$result[4],digits=3)

imp.la.risk <- with(imp,
                    glm(histo_num ~ ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(year,df=4)+ as.factor(sex)  + offset(log(value1))  + as.factor(SITE) + as.factor(art_status) + as.factor(mig_ende),
                        family = poisson(link = "log")))

a_risk <- D1(imp.la, imp.la.risk)
p_risk <- round(a_risk$result[4],digits=4)

imp.la.site <- with(imp,
                    glm(histo_num ~ ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(year,df=4)+ as.factor(sex)  + offset(log(value1)) + as.factor(risk1)  + as.factor(art_status) + as.factor(mig_ende), 
                        family = poisson(link = "log")))

a_site <- D1(imp.la,imp.la.site)
p_site <- round(a_site$result[4],digits=3)


imp.la.art <- with(imp,
                   glm(histo_num ~   ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(year,df=4)+ as.factor(sex)  + offset(log(value1)) + as.factor(risk1) + as.factor(SITE)  + as.factor(mig_ende), 
                       family = poisson(link = "log")))


a_art <- D1(imp.la, imp.la.art)
p_art <- round(a_art$result[4],digits=3)


imp.la.mig <- with(imp,
                   glm(histo_num ~   ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(year,df=4)+ as.factor(sex)  + offset(log(value1)) + as.factor(risk1) + as.factor(SITE) + as.factor(art_status), 
                       family = poisson(link = "log")))


a_mig <- D1(imp.la, imp.la.mig)
p_mig <- round(a_mig$result[4],digits=3)

## Variance covariance
t2 <- pool(imp.la)
m <- t2$m

## Applying rubin's rule
vcov_list <- list()
for (i in 1:m) {
  vcov_list[[i]] <- vcovHC(imp.la$analyses[[i]], cluster = imp$data$patient_id)
} 

vw <- Reduce('+',vcov_list)/m
#vcov_array <- array(unlist(vcov_list), dim = c(dim(vcov_list[[1]]),length(vcov_list)))
## Computing within impu2taion variance 
#var_vcov <- apply(vcov_array,c(1,2),var)
#vw <- var_vcov
## Between impu2taion variance covariance
qbar <- getqbar(t2)
qhats <- sapply(imp.la$analyses, coef)
vb <- (1 / (m-1)) * (qhats - qbar) %*% t(qhats - qbar)
vt <- vw + (1 + 1 / (m)) * vb 
vt <- vt[-c(1), -c(1)]


impu2.la <- function(model,var) {
  s <-  summary(pool(model))
  s <-  s$estimate[-1]
  sand <- sqrt(diag(var))
  r = matrix(0,nrow = length(s), ncol = 4)
  for (i in 1:(length(s))){
    r[i,] <- c(exp(s[i]), exp(s[i] - 1.96*sand[i]),exp(s[i] + 1.96*sand[i]),round((1 - pnorm(abs(s[i]/sand[i]))) * 2,digits=5))
  }
  
  r <- as.data.frame(r)
  colnames(r) <- c( "RR","Lower CI", "Upper CI","p-value") 
  rownam <- rownames(s)
  rownames(r) <- rownam
  r$ind <- seq_len(nrow(r))
  r <- rbind(r,c(NA,NA,NA,p_gen,ind=14.1),c(NA,NA,NA,p_risk,ind=15.1),c(NA,NA,NA,p_art,ind=18.1),c(NA,NA,NA,p_mig,ind=19.1),c(NA,NA,NA,p_site,ind=20.1))
  r <- r %>% arrange(ind) %>% dplyr::select(-ind)
  return(r)
}


impu2.la <- impu2.la(imp.la,vt)
cd4s <- seq(min(ld$sq_cd4,na.rm = T),max(ld$sq_cd4,na.rm=T), by = 0.01)
sp_cd4 <- as.data.frame(ns(cd4s, df=4))
sp_cd4 <- cbind(cd4s, sp_cd4) 
sp_cd4 <- as.data.frame(sp_cd4)
rownames(sp_cd4) <- sp_cd4$cd4
sp_cd4 <- sp_cd4[, -1]
colnames(sp_cd4) <- c("X1", "X2", "X3", "X4")

ages <- round(min(ld$age_b)):round(max(ld$age_b))
sp <- as.data.frame(ns(ages, df=4))
sp <- cbind(ages, sp) 
sp <- as.data.frame(sp)
rownames(sp) <- sp$age
sp <- sp[, -1]
colnames(sp) <- c("X1", "X2", "X3", "X4")

yr <- round(min(ld$year,na.rm = TRUE)):round(max(ld$year,na.rm = TRUE))
sp_yr <- as.data.frame(ns(yr, df=4))
sp_yr <- cbind(yr, sp_yr) 
sp_yr <- as.data.frame(sp_yr)
rownames(sp_yr) <- sp_yr$yr
sp_yr <- sp_yr[, -1]
colnames(sp_yr) <- c("X1", "X2", "X3", "X4")

rna <- round(min(ld$log_rna,na.rm = TRUE)):round(max(ld$log_rna,na.rm = TRUE))
sp_rna <- as.data.frame(ns(rna, df=4))
sp_rna <- cbind(rna,sp_rna) 
sp_rna <- as.data.frame(sp_rna)
rownames(sp_rna) <- sp_rna$rna
sp_rna <- sp_rna[, -1]
colnames(sp_rna) <- c("X1", "X2", "X3", "X4")

summ_spline  <- summary(pool(imp.la))
v_spline <- vt
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- summ_spline[-1,c(1,2)]
summ_spline <- as.data.frame(cbind(summ_spline, se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Variables","Coefs", "SE" )

rr_age_20 <- rr(20,summ_spline)
rr_age_25 <- rr(25,summ_spline)
rr_age_30 <- rr(30,summ_spline)
rr_age_40 <- rr(40,summ_spline)
rr_age_45 <- rr(45,summ_spline)
rr_age_50 <- rr(50,summ_spline)
rr_age_55 <- rr(55,summ_spline)
rr_age_60 <- rr(60,summ_spline)

SE_20 <- spline_se_fun(20)
SE_25 <- spline_se_fun(25)
SE_30 <- spline_se_fun(30)
SE_40 <- spline_se_fun(40)
SE_45 <- spline_se_fun(45)
SE_50 <- spline_se_fun(50)
SE_55 <- spline_se_fun(55)
SE_60 <- spline_se_fun(60)

CI_20 <- CI(rr_age_20,SE_20)
CI_25 <- CI(rr_age_25,SE_25)
CI_30 <- CI(rr_age_30,SE_30)
CI_40 <- CI(rr_age_40,SE_40)
CI_45 <- CI(rr_age_45,SE_45)
CI_50 <- CI(rr_age_50,SE_50)
CI_55 <- CI(rr_age_55,SE_55)
CI_60 <- CI(rr_age_60,SE_60)

abs_20 <- abs_val(20,summ_spline)
abs_25 <- abs_val(25,summ_spline)
abs_30 <- abs_val(30,summ_spline)
abs_40 <- abs_val(40,summ_spline)
abs_45 <- abs_val(45,summ_spline)
abs_50 <- abs_val(50,summ_spline)
abs_55 <- abs_val(55,summ_spline)
abs_60 <- abs_val(60,summ_spline)

rr_sq_5 <- rr_sq_sh(5,summ_spline)
rr_sq_10 <- rr_sq_sh(10,summ_spline)
rr_sq_18.70 <- rr_sq_sh(18.70,summ_spline)
rr_sq_22.40 <- rr_sq_sh(22.40,summ_spline)


SE_sq5 <- spline_se_sh(5)
SE_sq10 <- spline_se_sh(10)
SE_sq18.70 <- spline_se_sh(18.70)
SE_sq22.40 <- spline_se_sh(22.40)

CI_sq5 <- CI(rr_sq_5,SE_5)
CI_sq10 <- CI(rr_sq_10,SE_10)
CI_sq18.70 <- CI(rr_sq_18.70,SE_18.70)
CI_sq22.40 <- CI(rr_sq_22.40,SE_22.40)

abs_sq5 <- abs_sh(5,summ_spline)
abs_sq10 <- abs_sh(10,summ_spline)
abs_sq18.70 <- abs_sh(18.70,summ_spline)
abs_sq22.40 <- abs_sh(22.40,summ_spline)


rr_yr_2005 <- rr_yr2(2005,summ_spline)
rr_yr_2000 <- rr_yr2(2000,summ_spline)
rr_yr_2015 <- rr_yr2(2015,summ_spline)
rr_yr_2020 <- rr_yr2(2020,summ_spline)

SE_2005 <- spline_se_yr2(2005)
SE_2000 <- spline_se_yr2(2000)
SE_2015 <- spline_se_yr2(2015)
SE_2020 <- spline_se_yr2(2020)


CI_2005 <- CI(rr_yr_2005,SE_2005)
CI_2000 <- CI(rr_yr_2000,SE_2000)
CI_2015 <- CI(rr_yr_2015,SE_2015)
CI_2020 <- CI(rr_yr_2020,SE_2020)



abs_2005 <- abs_val_yrs2(2005,summ_spline)
abs_2000 <- abs_val_yrs2(2000,summ_spline)
abs_2015 <- abs_val_yrs2(2015,summ_spline)
abs_2020 <- abs_val_yrs2(2020,summ_spline)

rr_rna_1 <- rr_rna(1,summ_spline)
rr_rna_3 <- rr_rna(3,summ_spline)
rr_rna_7 <- rr_rna(7,summ_spline)


SE_1 <- spline_se_rna(1)
SE_3 <- spline_se_rna(3)
SE_7 <- spline_se_rna(7)



CI_1 <- CI(rr_rna_1,SE_1)
CI_3 <- CI(rr_rna_3,SE_3)
CI_7 <- CI(rr_rna_7,SE_7)




abs_1 <- abs_val_rna(1,summ_spline)
abs_3 <- abs_val_rna(3,summ_spline)
abs_7 <- abs_val_rna(7,summ_spline)

age_spl_multi <- data.frame("RR" = c(rr_age_20,rr_age_25,rr_age_30,rr_age_40,rr_age_45,rr_age_50,rr_age_55,rr_age_60),
                            "Lower CI"= c(CI_20[1],CI_25[1],CI_30[1],CI_40[1],CI_45[1],CI_50[1],CI_55[1],CI_60[1]),
                            "Upper CI"= c(CI_20[2],CI_25[2],CI_30[2],CI_40[2],CI_45[2],CI_50[2],CI_55[2],CI_60[2]),
                            "p-value" = c(2*(1-pnorm(abs(abs_20/SE_20))), 2*(1-pnorm(abs(abs_25/SE_25))), 2*(1-pnorm(abs(abs_30/SE_30))),
                                          2*(1-pnorm(abs(abs_40/SE_40))), 2*(1-pnorm(abs(abs_45/SE_45))), 2*(1-pnorm(abs(abs_50/SE_50))),
                                          2*(1-pnorm(abs(abs_55/SE_55))), 2*(1-pnorm(abs(abs_60/SE_60)))))
colnames(age_spl_multi) <- c("RR","Lower CI","Upper CI","p-value")


cd4_spline_multi <- data.frame("RR" = c(rr_sq_5,rr_sq_10,rr_sq_18.70,rr_sq_22.40),
                               "Lower CI"= c(CI_sq5[1],CI_sq10[1],CI_sq18.70[1],CI_sq22.40[1]),
                               "Upper CI"= c(CI_sq5[2],CI_sq10[2],CI_sq18.70[2],CI_sq22.40[2]),
                               "p-value" = c(2*(1-pnorm(abs(abs_sq5/SE_sq5))), 2*(1-pnorm(abs(abs_sq10/SE_sq10))), 2*(1-pnorm(abs(abs_sq18.70/SE_sq18.70))),
                                             2*(1-pnorm(abs(abs_sq22.40/SE_sq22.40)))))
colnames(cd4_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")


yr_spline_multi <- data.frame("RR" = c(rr_yr_2000,rr_yr_2005,rr_yr_2015,rr_yr_2020),
                              "Lower CI"= c(CI_2000[1],CI_2005[1],CI_2015[1],CI_2020[1]),
                              "Upper CI"= c(CI_2000[2],CI_2005[2],CI_2015[2],CI_2020[2]),
                              "p-value" = c(2*(1-pnorm(abs(abs_2000/SE_2005))), 2*(1-pnorm(abs(abs_2005/SE_2000))),
                                            2*(1-pnorm(abs(abs_2015/SE_2015))), 2*(1-pnorm(abs(abs_2020/SE_2020)))))
colnames(yr_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")

rna_spline_multi <- data.frame("RR" = c(rr_rna_1,rr_rna_3,rr_rna_7),
                               "Lower CI"= c(CI_1[1],CI_3[1],CI_7[1]),
                               "Upper CI"= c(CI_1[2],CI_3[2],CI_7[2]),
                               "p-value" = c(2*(1-pnorm(abs(abs_1/SE_1))), 2*(1-pnorm(abs(abs_3/SE_3))),
                                             2*(1-pnorm(abs(abs_7/SE_7)))))
colnames(rna_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")

impu2.la <- impu2.la[-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14),]
tab5_extra <- rbind(c(NA,NA,NA,p_age),age_spl_multi,c(NA,NA,NA,p_cd),cd4_spline_multi,c(NA,NA,NA,p_rna),rna_spline_multi,c(NA,NA,NA,p_yr),yr_spline_multi,impu2.la)
rownames(tab5_extra) <- c("Age (ref=35)","20","25","30","40","45","50","55","60",
                          "CD4 (ref=200)","25 ","100","350","500","log HIV RNA(ref=5)"," 1",  "3","7"
                          ,"Year (ref=2010)", "2000","2005","2015","2020","Sex","Male ", "Risk factor (ref =MSM)","Hetero ","IDU ", "Other/Unknown ",
                          "ART initiation status (ref=post)","pre-ART initiation", "Migration status (ref=No)","Yes",
                          "Country (ref=Brazil)", "Chile ", "Honduras ","Mexico ", "Peru ")

save(tab5_extra, file="tab5_extra.Rdata")
## imputation model Tn ##

imp_t <- tn_t %>% select(CFAR_PID,age_b,sq_cd4,log_rna, baseline_y,birthsex,risk1,art_status,mig_ende,follow_up,histo_num)
imp_t <- mice(imp_t,m=5,printFlag = FALSE, idvars = "patient_id")
set.seed(7)

#tn_imp <- aregImpute(~ age_b + sq_cd4+ log_rna + baseline_y+ as.factor(birthsex)  + as.factor(risk1)  + as.factor(art_status) + as.factor(mig_ende) + follow_up, data = tn_t, n.impute=20,x=TRUE, nk=0)


imp.tn <- with(imp_t,
               glm(histo_num ~  age_b + sq_cd4 + log_rna + baseline_y + as.factor(birthsex)  + offset(log(follow_up)) + as.factor(risk1)  + as.factor(art_status) + as.factor(mig_ende),
                   family = poisson(link = "log")))

imp.tn.age <-  with(imp_t,
                    glm(histo_num ~  sq_cd4 + log_rna + baseline_y + as.factor(birthsex)  + offset(log(follow_up)) + as.factor(risk1)  + as.factor(art_status) + as.factor(mig_ende), 
                        family = poisson(link = "log")))        

a_tng <- D1(imp.tn,imp.tn.age)
p_tng <- round(a_tng$result[4],digits=3)

imp.tn.cd <- with(imp_t,
                  glm(histo_num ~  age_b + log_rna + baseline_y + as.factor(birthsex)  + offset(log(follow_up)) + as.factor(risk1)  + as.factor(art_status) + as.factor(mig_ende), 
                      family = poisson(link = "log")))

a_tncd <- D1(imp.tn,imp.tn.cd)
p_tncd <- round(a_tncd$result[4],digits=3)

imp.tn.rna <- with(imp_t,
                   glm(histo_num ~ age_b + sq_cd4 + baseline_y + as.factor(birthsex)  + offset(log(follow_up)) + as.factor(risk1)  + as.factor(art_status) + as.factor(mig_ende), 
                       family = poisson(link = "log")))

a_tnrna <- D1(imp.tn,imp.tn.rna)
p_tnrna <- round(a_tnrna$result[4],digits=3)

imp.tn.yr <- with(imp_t,
                  glm(histo_num ~ age_b + sq_cd4+ log_rna + as.factor(birthsex)  + offset(log(follow_up)) + as.factor(risk1)  + as.factor(art_status) + as.factor(mig_ende), 
                      family = poisson(link = "log")))

a_tnyr <- D1(imp.tn,imp.tn.yr)
p_tnyr <- round(a_tnyr$result[4],digits=3)


imp.tn.gen <- with(imp_t,
                   glm(histo_num ~ age_b + sq_cd4+ log_rna + baseline_y   + offset(log(follow_up)) + as.factor(risk1)  + as.factor(art_status) + as.factor(mig_ende), 
                       family = poisson(link = "log")))

a_tngen <- D1(imp.tn,imp.tn.gen)
p_tngen <- round(a_tngen$result[4],digits=3)

imp.tn.risk <- with(imp_t,
                    glm(histo_num ~ age_b + sq_cd4+ log_rna + baseline_y + as.factor(birthsex)  + offset(log(follow_up))   + as.factor(art_status) + as.factor(mig_ende), 
                        family = poisson(link = "log")))

a_tnrisk <- D1(imp.tn,imp.tn.risk)
p_tnrisk <- round(a_tnrisk$result[4],digits=3)


imp.tn.art <- with(imp_t,
                   glm(histo_num ~  age_b + sq_cd4+ log_rna + baseline_y + as.factor(birthsex)  + offset(log(follow_up)) + as.factor(risk1)   + as.factor(mig_ende), 
                       family = poisson(link = "log")))


a_tnart <- D1(imp.tn,imp.tn.art)
p_tnart <- round(a_tnart$result[4],digits=3)


imp.tn.mig <- with(imp_t,
                   glm(histo_num ~  age_b + sq_cd4+ log_rna + baseline_y + as.factor(birthsex)  + offset(log(follow_up)) + as.factor(risk1)  + as.factor(art_status), 
                       family = poisson(link = "log")))


a_tnmig <- D1(imp.tn,imp.tn.mig)
p_tnmig <- round(a_tnmig$result[4],digits=3)

## Variance covariance
t2 <- pool(imp.tn)
m <- t2$m

## Applying rubin's rule
vcov_list <- list()
for (i in 1:m) {
  vcov_list[[i]] <- vcovHC(imp.tn$analyses[[i]], cluster = imp$data$patient_id)
} 

vw <- Reduce('+',vcov_list)/m
#vcov_array <- array(unlist(vcov_list), dim = c(dim(vcov_list[[1]]),length(vcov_list)))
## Computing within imputaion variance 
#var_vcov <- apply(vcov_array,c(1,2),var)
#vw <- var_vcov
## Between imputaion variance covariance
qbar <- getqbar(t2)
qhats <- sapply(imp.tn$analyses, coef)
vb <- (1 / (m-1)) * (qhats - qbar) %*% t(qhats - qbar)
vt <- vw + (1 + 1 / (m)) * vb 
vt <- vt[-c(1), -c(1)]

impu.tn <- function(model,var1) {
  s <-  summary(pool(model))
  s <-  s$estimate[-1]
  sand <- sqrt(diag(var1))
  r = matrix(0,nrow = 24, ncol = 4)
  r[1,] <- c(NA,NA,NA,p_tng)
  r[2,] <- c(exp(s[1]*(-10)), exp((s[1] + 1.96*sand[1])*(-10)),exp((s[1] - 1.96*sand[1])*(-10)),round((1 - pnorm(abs(s[1]/sand[1]))) * 2,digits=5))
  r[3,] <- c(exp(s[1]*10), exp((s[1] - 1.96*sand[1])*10),exp((s[1] + 1.96*sand[1])*10),round((1 - pnorm(abs(s[1]/sand[1]))) * 2,digits=5))
  r[4,] <- c(exp(s[1]*20), exp((s[1] - 1.96*sand[1])*20),exp((s[1] + 1.96*sand[1])*20),round((1 - pnorm(abs(s[1]/sand[1]))) * 2,digits=5))
  
  r[5,] <- c(NA,NA,NA,p_tncd)
  r[6,] <- c(exp(s[2]*(-9.14)), exp((s[2] + 1.96*sand[2])*(-9.14)),exp((s[2] - 1.96*sand[2])*(-9.14)),round((1 - pnorm(abs(s[2]/sand[2]))) * 2,digits=5))
  r[7,] <- c(exp(s[2]*(-4.14)), exp((s[2] + 1.96*sand[2])*(-4.14)),exp((s[2] - 1.96*sand[2])*(-4.14)),round((1 - pnorm(abs(s[2]/sand[2]))) * 2,digits=5))
  r[8,] <- c(exp(s[2]*(4.56)), exp((s[2] - 1.96*sand[2])*(4.56)),exp((s[2] + 1.96*sand[2])*(4.56)),round((1 - pnorm(abs(s[2]/sand[2]))) * 2,digits=5))
  r[9,] <- c(exp(s[2]*(8.218)), exp((s[2] - 1.96*sand[2])*(8.218)),exp((s[2] + 1.96*sand[2])*(8.218)),round((1 - pnorm(abs(s[2]/sand[2]))) * 2,digits=5))
  
  r[10,] <- c(NA,NA,NA,p_tnrna)
  r[11,] <- c(exp(s[3]*(-4)), exp((s[3] + 1.96*sand[3])*(-4)),exp((s[3] - 1.96*sand[3])*(-4)),round((1 - pnorm(abs(s[3]/sand[3]))) * 2,digits=5))
  r[12,] <- c(exp(s[3]*(-2)), exp((s[3] + 1.96*sand[3])*(-2)),exp((s[3] - 1.96*sand[3])*(-2)),round((1 - pnorm(abs(s[3]/sand[3]))) * 2,digits=5))
  r[13,] <- c(exp(s[3]*2), exp((s[3] - 1.96*sand[3])*2),exp((s[3] + 1.96*sand[3])*2),round((1 - pnorm(abs(s[3]/sand[3]))) * 2,digits=5))
  r[14,] <- c(NA,NA,NA,p_tnyr)
  r[15,] <- c(exp(s[4]*(-10)), exp((s[4] + 1.96*sand[4])*(-10)),exp((s[4] - 1.96*sand[4])*(-10)),round((1 - pnorm(abs(s[4]/sand[4]))) * 2,digits=5))
  r[16,] <-  c(exp(s[4]*(-5)), exp((s[4] + 1.96*sand[4])*(-5)),exp((s[4] - 1.96*sand[4])*(-5)),round((1 - pnorm(abs(s[4]/sand[4]))) * 2,digits=5))
  r[17,] <- c(exp(s[4]*5), exp((s[4] - 1.96*sand[4])*5),exp((s[4] + 1.96*sand[4])*5),round((1 - pnorm(abs(s[4]/sand[4]))) * 2,digits=5))
  r[18,] <- c(exp(s[4]*10), exp((s[4] - 1.96*sand[4])*10),exp((s[4] + 1.96*sand[4])*10),round((1 - pnorm(abs(s[4]/sand[4]))) * 2,digits=5))
  
  for (i in 5:10){
    r[i+14,] <- c(exp(s[i]), exp(s[i] - 1.96*sand[i]),exp(s[i] + 1.96*sand[i]),round((1 - pnorm(abs(s[i]/sand[i]))) * 2,digits=5))
  }
  
  r <- as.data.frame(r)
  colnames(r) <- c( "RR","Lower CI", "Upper CI","p-value") 
  rownam <- rownames(s)
  rownames(r) <- rownam
  r$ind <- seq_len(nrow(r))
  r <- rbind(r,c(NA,NA,NA,p_tngen,ind=18.1),c(NA,NA,NA,p_tnrisk,ind=19.1),c(NA,NA,NA,p_tnart,ind=22.1),c(NA,NA,NA,p_tnmig,ind=23.1))
  r <- r %>% arrange(ind) %>% dplyr::select(-ind)
  return(r)
}


tab6 <- impu.tn(imp.tn,vt)

rownames(tab6) <- c("Age (ref=35)","25","45","55",
                    "CD4 (ref=200)","25  ","  100 ","350","500  ","log HIV RNA(ref=  5 )", "1", "3" ,"7",
                    "Year (ref=2010)", "2000","2005","2015","2020","Sex","Male ", "Risk factor (ref =MSM)","Hetero ","IDU ", "Other/Unknown ",
                    "ART initiation status (ref=post)","pre-ART initiation", "Migration status (ref=No)","Yes")

save(tab6, file="tab6.Rdata")
incid_new <- incid_tn_long %>% select(CFAR_PID,value,year,follow_time_yrs)
tn_d <- full_join(tn_t,incid_new, by = "CFAR_PID") %>% select(CFAR_PID,value,year,follow_time_yrs,art_sd,last_d,baseline,age_b,sq_cd4,log_rna,baseline_y,birthsex,follow_up,risk1,art_status,mig_ende,
                                                              histo,histo_d,histo_num)
#tn_d <- tn_d %>% distinct("CFAR_PID", .keep_all = T)
tn_d$art_y <- year(tn_d$art_sd)
#tn_d <- tn_d %>% filter(!is.na(art_y))
tn_d$year_end_d <- as.Date(paste(tn_d$year,12,31, sep="-"))
tn_d$year_beg_d <- as.Date(paste(tn_d$year,01,01, sep="-"))
tn_d <- tn_d %>% group_by(CFAR_PID) %>% mutate(d1 = pmin(baseline,art_sd, na.rm=T))
tn_d <- tn_d %>% group_by(CFAR_PID) %>% mutate(d2 = pmax(baseline,art_sd, na.rm=T))
tn_d$t1 <- difftime(tn_d$d2,tn_d$d1,units = "days")/365.25
tn_d$t2 <- difftime(tn_d$year_end_d,tn_d$d2,units = "days")/365.25
tn_d$indi <- with(tn_d,ifelse(year == art_y & !is.na(art_y),1,0))

tn_d <- tn_d %>% group_by(CFAR_PID) %>% dplyr::mutate(indi_row = cumsum(indi),
                                                      value1 = ifelse(indi ==1 & indi_row ==1 ,t1,
                                                                      ifelse(indi == 1 & indi_row == 2,t2,NA)))

tn_d$value1 <- with(tn_d,ifelse(is.na(value1),value,value1))

tn_d <- tn_d %>% distinct(CFAR_PID,value1,year, .keep_all = T)

tn_d <- tn_d %>% group_by(CFAR_PID) %>% arrange(year,value1)


## redefining art_atatus
tn_d$indi <- with(tn_d,ifelse(art_status == "post-ART initiation",1,0))
tn_d <- tn_d %>% group_by(CFAR_PID) %>% dplyr::mutate(indi_row = cumsum(indi))
tn_d$art_status <- with(tn_d,ifelse(indi_row >= 1,"post-ART initiation","pre-ART initiation"))
tn_d$value1 <- with(tn_d,ifelse(value1 == 0,0.002737851,value1))
# ms1 <- glm(histo_num ~  ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(year,df=4)+ as.factor(birthsex)  + offset(log(value1))
#            + as.factor(risk1) + as.factor(art_status) + as.factor(mig_ende),family = poisson(link = "log"), data = tn_d)
imp_t <- tn_d %>% select(CFAR_PID,age_b,sq_cd4,log_rna, year,birthsex,risk1,art_status,mig_ende,follow_up,histo_num)
imp_t <- mice(imp_t,m=5,printFlag = FALSE, idvars = "patient_id")
set.seed(7)

#tn_imp <- aregImpute(~ age_b + sq_cd4+ log_rna + baseline_y+ as.factor(birthsex)  + as.factor(risk1)  + as.factor(art_status) + as.factor(mig_ende) + follow_up, data = tn_d, n.impute=20,x=TRUE, nk=0)


imp.tn <- with(imp_t,
               glm(histo_num ~  age_b + sq_cd4 + log_rna + year + as.factor(birthsex)  + offset(log(follow_up)) + as.factor(risk1)  + as.factor(art_status) + as.factor(mig_ende),
                   family = poisson(link = "log")))

imp.tn.age <-  with(imp_t,
                    glm(histo_num ~  sq_cd4 + log_rna + year + as.factor(birthsex)  + offset(log(follow_up)) + as.factor(risk1)  + as.factor(art_status) + as.factor(mig_ende), 
                        family = poisson(link = "log")))        

a_tng <- D1(imp.tn,imp.tn.age)
p_tng <- round(a_tng$result[4],digits=3)

imp.tn.cd <- with(imp_t,
                  glm(histo_num ~  age_b + log_rna + year + as.factor(birthsex)  + offset(log(follow_up)) + as.factor(risk1)  + as.factor(art_status) + as.factor(mig_ende), 
                      family = poisson(link = "log")))

a_tncd <- D1(imp.tn,imp.tn.cd)
p_tncd <- round(a_tncd$result[4],digits=3)

imp.tn.rna <- with(imp_t,
                   glm(histo_num ~ age_b + sq_cd4 + year + as.factor(birthsex)  + offset(log(follow_up)) + as.factor(risk1)  + as.factor(art_status) + as.factor(mig_ende), 
                       family = poisson(link = "log")))

a_tnrna <- D1(imp.tn,imp.tn.rna)
p_tnrna <- round(a_tnrna$result[4],digits=3)

imp.tn.yr <- with(imp_t,
                  glm(histo_num ~ age_b + sq_cd4+ log_rna + as.factor(birthsex)  + offset(log(follow_up)) + as.factor(risk1)  + as.factor(art_status) + as.factor(mig_ende), 
                      family = poisson(link = "log")))

a_tnyr <- D1(imp.tn,imp.tn.yr)
p_tnyr <- round(a_tnyr$result[4],digits=3)


imp.tn.gen <- with(imp_t,
                   glm(histo_num ~ age_b + sq_cd4+ log_rna + year   + offset(log(follow_up)) + as.factor(risk1)  + as.factor(art_status) + as.factor(mig_ende), 
                       family = poisson(link = "log")))

a_tngen <- D1(imp.tn,imp.tn.gen)
p_tngen <- round(a_tngen$result[4],digits=3)

imp.tn.risk <- with(imp_t,
                    glm(histo_num ~ age_b + sq_cd4+ log_rna + year + as.factor(birthsex)  + offset(log(follow_up))   + as.factor(art_status) + as.factor(mig_ende), 
                        family = poisson(link = "log")))

a_tnrisk <- D1(imp.tn,imp.tn.risk)
p_tnrisk <- round(a_tnrisk$result[4],digits=3)


imp.tn.art <- with(imp_t,
                   glm(histo_num ~  age_b + sq_cd4+ log_rna + year + as.factor(birthsex)  + offset(log(follow_up)) + as.factor(risk1)   + as.factor(mig_ende), 
                       family = poisson(link = "log")))


a_tnart <- D1(imp.tn,imp.tn.art)
p_tnart <- round(a_tnart$result[4],digits=3)


imp.tn.mig <- with(imp_t,
                   glm(histo_num ~  age_b + sq_cd4+ log_rna + year + as.factor(birthsex)  + offset(log(follow_up)) + as.factor(risk1)  + as.factor(art_status), 
                       family = poisson(link = "log")))


a_tnmig <- D1(imp.tn,imp.tn.mig)
p_tnmig <- round(a_tnmig$result[4],digits=3)

## Variance covariance
t2 <- pool(imp.tn)
m <- t2$m

## Applying rubin's rule
vcov_list <- list()
for (i in 1:m) {
  vcov_list[[i]] <- vcovHC(imp.tn$analyses[[i]], cluster = imp$data$patient_id)
} 

vw <- Reduce('+',vcov_list)/m
#vcov_array <- array(unlist(vcov_list), dim = c(dim(vcov_list[[1]]),length(vcov_list)))
## Computing within imputaion variance 
#var_vcov <- apply(vcov_array,c(1,2),var)
#vw <- var_vcov
## Between imputaion variance covariance
qbar <- getqbar(t2)
qhats <- sapply(imp.tn$analyses, coef)
vb <- (1 / (m-1)) * (qhats - qbar) %*% t(qhats - qbar)
vt <- vw + (1 + 1 / (m)) * vb 
vt <- vt[-c(1), -c(1)]

impu.tn <- function(model,var1) {
  s <-  summary(pool(model))
  s <-  s$estimate[-1]
  sand <- sqrt(diag(var1))
  r = matrix(0,nrow = 24, ncol = 4)
  r[1,] <- c(NA,NA,NA,p_tng)
  r[2,] <- c(exp(s[1]*(-10)), exp((s[1] + 1.96*sand[1])*(-10)),exp((s[1] - 1.96*sand[1])*(-10)),round((1 - pnorm(abs(s[1]/sand[1]))) * 2,digits=5))
  r[3,] <- c(exp(s[1]*10), exp((s[1] - 1.96*sand[1])*10),exp((s[1] + 1.96*sand[1])*10),round((1 - pnorm(abs(s[1]/sand[1]))) * 2,digits=5))
  r[4,] <- c(exp(s[1]*20), exp((s[1] - 1.96*sand[1])*20),exp((s[1] + 1.96*sand[1])*20),round((1 - pnorm(abs(s[1]/sand[1]))) * 2,digits=5))
  
  r[5,] <- c(NA,NA,NA,p_tncd)
  r[6,] <- c(exp(s[2]*(-9.14)), exp((s[2] + 1.96*sand[2])*(-9.14)),exp((s[2] - 1.96*sand[2])*(-9.14)),round((1 - pnorm(abs(s[2]/sand[2]))) * 2,digits=5))
  r[7,] <- c(exp(s[2]*(-4.14)), exp((s[2] + 1.96*sand[2])*(-4.14)),exp((s[2] - 1.96*sand[2])*(-4.14)),round((1 - pnorm(abs(s[2]/sand[2]))) * 2,digits=5))
  r[8,] <- c(exp(s[2]*(4.56)), exp((s[2] - 1.96*sand[2])*(4.56)),exp((s[2] + 1.96*sand[2])*(4.56)),round((1 - pnorm(abs(s[2]/sand[2]))) * 2,digits=5))
  r[9,] <- c(exp(s[2]*(8.218)), exp((s[2] - 1.96*sand[2])*(8.218)),exp((s[2] + 1.96*sand[2])*(8.218)),round((1 - pnorm(abs(s[2]/sand[2]))) * 2,digits=5))
  
  r[10,] <- c(NA,NA,NA,p_tnrna)
  r[11,] <- c(exp(s[3]*(-4)), exp((s[3] + 1.96*sand[3])*(-4)),exp((s[3] - 1.96*sand[3])*(-4)),round((1 - pnorm(abs(s[3]/sand[3]))) * 2,digits=5))
  r[12,] <- c(exp(s[3]*(-2)), exp((s[3] + 1.96*sand[3])*(-2)),exp((s[3] - 1.96*sand[3])*(-2)),round((1 - pnorm(abs(s[3]/sand[3]))) * 2,digits=5))
  r[13,] <- c(exp(s[3]*2), exp((s[3] - 1.96*sand[3])*2),exp((s[3] + 1.96*sand[3])*2),round((1 - pnorm(abs(s[3]/sand[3]))) * 2,digits=5))
  r[14,] <- c(NA,NA,NA,p_tnyr)
  r[15,] <- c(exp(s[4]*(-10)), exp((s[4] + 1.96*sand[4])*(-10)),exp((s[4] - 1.96*sand[4])*(-10)),round((1 - pnorm(abs(s[4]/sand[4]))) * 2,digits=5))
  r[16,] <-  c(exp(s[4]*(-5)), exp((s[4] + 1.96*sand[4])*(-5)),exp((s[4] - 1.96*sand[4])*(-5)),round((1 - pnorm(abs(s[4]/sand[4]))) * 2,digits=5))
  r[17,] <- c(exp(s[4]*5), exp((s[4] - 1.96*sand[4])*5),exp((s[4] + 1.96*sand[4])*5),round((1 - pnorm(abs(s[4]/sand[4]))) * 2,digits=5))
  r[18,] <- c(exp(s[4]*10), exp((s[4] - 1.96*sand[4])*10),exp((s[4] + 1.96*sand[4])*10),round((1 - pnorm(abs(s[4]/sand[4]))) * 2,digits=5))
  
  for (i in 5:10){
    r[i+14,] <- c(exp(s[i]), exp(s[i] - 1.96*sand[i]),exp(s[i] + 1.96*sand[i]),round((1 - pnorm(abs(s[i]/sand[i]))) * 2,digits=5))
  }
  
  r <- as.data.frame(r)
  colnames(r) <- c( "RR","Lower CI", "Upper CI","p-value") 
  rownam <- rownames(s)
  rownames(r) <- rownam
  r$ind <- seq_len(nrow(r))
  r <- rbind(r,c(NA,NA,NA,p_tngen,ind=18.1),c(NA,NA,NA,p_tnrisk,ind=19.1),c(NA,NA,NA,p_tnart,ind=22.1),c(NA,NA,NA,p_tnmig,ind=23.1))
  r <- r %>% arrange(ind) %>% dplyr::select(-ind)
  return(r)
}


tab6_extra <- impu.tn(imp.tn,vt)

rownames(tab6_extra) <- c("Age (ref=35)","25","45","55",
                          "CD4 (ref=200)","25  ","  100 ","350","500  ","log HIV RNA(ref=  5 )", "1", "3" ,"7",
                          "Year (ref=2010)", "2000","2005","2015","2020","Sex","Male ", "Risk factor (ref =MSM)","Hetero ","IDU ", "Other/Unknown ",
                          "ART initiation status (ref=post)","pre-ART initiation", "Migration status (ref=No)","Yes")

save(tab6_extra, file="tab6_extra.Rdata")

# Relevel
la_t_pre$sex <- relevel(as.factor(la_t_pre$sex), ref = "Female")
la_t_pre$SITE <- relevel(as.factor(la_t_pre$SITE), ref = "Brazil")
la_t_pre$risk1 <- relevel(la_t_pre$risk1, ref = "MSM")
la_t_pre$mig_ende <- relevel(as.factor(la_t_pre$mig_ende), ref = "No")

## Non spline univariate models
uni0 <- glm(histo_num ~ offset(log(follow_up)),family = poisson(link = "log"), data = la_t_pre)
uni.gen <- glm(histo_num ~as.factor(sex) + offset(log(follow_up)),family = poisson(link = "log"), data = la_t_pre)
uni.risk <- glm(histo_num ~as.factor(risk1) + offset(log(follow_up)),family = poisson(link = "log"), data = la_t_pre)
uni.site <- glm(histo_num ~as.factor(SITE) + offset(log(follow_up)),family = poisson(link = "log"), data = la_t_pre)
uni.mig <- glm(histo_num ~as.factor(mig_ende) + offset(log(follow_up)),family = poisson(link = "log"), data = la_t_pre)

## Splines ##
ages <- round(min(la_t_pre$age_b)):round(max(la_t_pre$age_b))
sp <- as.data.frame(ns(ages, df=4))
sp <- cbind(ages, sp) 
sp <- as.data.frame(sp)
rownames(sp) <- sp$age
sp <- sp[, -1]
colnames(sp) <- c("X1", "X2", "X3", "X4")

uni.age <- glm(histo_num ~ ns(age_b,df =4) + offset(log(follow_up)),family = poisson(link =log), data = la_t_pre)
summ_spline  <- summary(uni.age)
v_spline <- sandwich(uni.age)
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Coefs", "SE", "Variables")

rr_age_20 <- rr(20,summ_spline)
rr_age_25 <- rr(25,summ_spline)
rr_age_30 <- rr(30,summ_spline)
rr_age_40 <- rr(40,summ_spline)
rr_age_45 <- rr(45,summ_spline)
rr_age_50 <- rr(50,summ_spline)
rr_age_55 <- rr(55,summ_spline)
rr_age_60 <- rr(60,summ_spline)


SE_20 <- spline_se_fun(20)
SE_25 <- spline_se_fun(25)
SE_30 <- spline_se_fun(30)
SE_40 <- spline_se_fun(40)
SE_45 <- spline_se_fun(45)
SE_50 <- spline_se_fun(50)
SE_55 <- spline_se_fun(55)
SE_60 <- spline_se_fun(60)

CI <- function(rr,SE) 
{ exp(log(rr) + c(-1,1)*1.96*SE)}


CI_20 <- CI(rr_age_20,SE_20)
CI_25 <- CI(rr_age_25,SE_25)
CI_30 <- CI(rr_age_30,SE_30)
CI_40 <- CI(rr_age_40,SE_40)
CI_45 <- CI(rr_age_45,SE_45)
CI_50 <- CI(rr_age_50,SE_50)
CI_55 <- CI(rr_age_55,SE_55)
CI_60 <- CI(rr_age_60,SE_60)

abs_20 <- abs_val(20,summ_spline)
abs_25 <- abs_val(25,summ_spline)
abs_30 <- abs_val(30,summ_spline)
abs_40 <- abs_val(40,summ_spline)
abs_45 <- abs_val(45,summ_spline)
abs_50 <- abs_val(50,summ_spline)
abs_55 <- abs_val(55,summ_spline)
abs_60 <- abs_val(60,summ_spline)

age_spline <- data.frame("RR" = c(rr_age_20,rr_age_25,rr_age_30,rr_age_40,rr_age_45,rr_age_50,rr_age_55,rr_age_60),
                         "Lower CI"= c(CI_20[1],CI_25[1],CI_30[1],CI_40[1],CI_45[1],CI_50[1],CI_55[1],CI_60[1]),
                         "Upper CI"= c(CI_20[2],CI_25[2],CI_30[2],CI_40[2],CI_45[2],CI_50[2],CI_55[2],CI_60[2]),
                         "p-value" = c(2*(1-pnorm(abs(abs_20/SE_20))), 2*(1-pnorm(abs(abs_25/SE_25))), 2*(1-pnorm(abs(abs_30/SE_30))),
                                       2*(1-pnorm(abs(abs_40/SE_40))), 2*(1-pnorm(abs(abs_45/SE_45))), 2*(1-pnorm(abs(abs_50/SE_50))),
                                       2*(1-pnorm(abs(abs_55/SE_55))), 2*(1-pnorm(abs(abs_60/SE_60)))))
colnames(age_spline) <- c("RR","Lower CI","Upper CI","p-value")


cd4s <- seq(min(la_t_pre$sq_cd4,na.rm = T),max(la_t_pre$sq_cd4,na.rm=T), by = 0.01)
sp_cd4 <- as.data.frame(ns(cd4s, df=4))
sp_cd4 <- cbind(cd4s, sp_cd4) 
sp_cd4 <- as.data.frame(sp_cd4)
rownames(sp_cd4) <- sp_cd4$cd4
sp_cd4 <- sp_cd4[, -1]
colnames(sp_cd4) <- c("X1", "X2", "X3", "X4")

uni.cd <- glm(histo_num ~ ns(sq_cd4,df =4) + offset(log(follow_up)),family = poisson(link =log), data = la_t_pre)
summ_spline  <- summary(uni.cd)
v_spline <- sandwich(uni.cd)
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Coefs", "SE", "Variables")

rr_sq_5 <- rr_sq(5,summ_spline)
rr_sq_10 <- rr_sq(10,summ_spline)
rr_sq_18.70 <- rr_sq(18.70,summ_spline)
rr_sq_22.40 <- rr_sq(22.40,summ_spline)


SE_5 <- spline_se_fun_sq_cd(5)
SE_10 <- spline_se_fun_sq_cd(10)
SE_18.70 <- spline_se_fun_sq_cd(18.70)
SE_22.40 <- spline_se_fun_sq_cd(22.40)



CI_5 <- CI(rr_sq_5,SE_5)
CI_10 <- CI(rr_sq_10,SE_10)
CI_18.70 <- CI(rr_sq_18.70,SE_18.70)
CI_22.40 <- CI(rr_sq_22.40,SE_22.40)



abs_5 <- abs_sq_cd(5,summ_spline)
abs_10 <- abs_sq_cd(10,summ_spline)
abs_18.70 <- abs_sq_cd(18.70,summ_spline)
abs_22.40 <- abs_sq_cd(22.40,summ_spline)



cd4_spline <- data.frame("RR" = c(rr_sq_5,rr_sq_10,rr_sq_18.70,rr_sq_22.40),
                         "Lower CI"= c(CI_5[1],CI_10[1],CI_18.70[1],CI_20[1]),
                         "Upper CI"= c(CI_5[2],CI_10[2],CI_18.70[2],CI_20[2]),
                         "p-value" = c(2*(1-pnorm(abs(abs_5/SE_5))), 2*(1-pnorm(abs(abs_10/SE_10))), 2*(1-pnorm(abs(abs_18.70/SE_18.70))),
                                       2*(1-pnorm(abs(abs_25/SE_22.40)))))
colnames(cd4_spline) <- c("RR","Lower CI","Upper CI","p-value")

yr <- round(min(la_t_pre$baseline_y,na.rm = TRUE)):round(max(la_t_pre$baseline_y,na.rm = TRUE))
sp_yr <- as.data.frame(ns(yr, df=4))
sp_yr <- cbind(yr, sp_yr) 
sp_yr <- as.data.frame(sp_yr)
rownames(sp_yr) <- sp_yr$yr
sp_yr <- sp_yr[, -1]
colnames(sp_yr) <- c("X1", "X2", "X3", "X4")

uni.yrs <- glm(histo_num ~ ns(baseline_y,df =4) + offset(log(follow_up)),family = poisson(link =log), data = la_t_pre)
summ_spline  <- summary(uni.yrs)
v_spline <- sandwich(uni.yrs)
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Coefs", "SE", "Variables")


rr_yr_2005 <- rr_yr(2005,summ_spline)
rr_yr_2000 <- rr_yr(2000,summ_spline)
rr_yr_2015 <- rr_yr(2015,summ_spline)
rr_yr_2020 <- rr_yr(2020,summ_spline)

SE_2005 <- spline_se_yr(2005)
SE_2000 <- spline_se_yr(2000)
SE_2015 <- spline_se_yr(2015)
SE_2020 <- spline_se_yr(2020)


CI_2005 <- CI(rr_yr_2005,SE_2005)
CI_2000 <- CI(rr_yr_2000,SE_2000)
CI_2015 <- CI(rr_yr_2015,SE_2015)
CI_2020 <- CI(rr_yr_2020,SE_2020)



abs_2005 <- abs_val_yrs(2005,summ_spline)
abs_2000 <- abs_val_yrs(2000,summ_spline)
abs_2015 <- abs_val_yrs(2015,summ_spline)
abs_2020 <- abs_val_yrs(2020,summ_spline)



yr_spline <- data.frame("RR" = c(rr_yr_2005,rr_yr_2000,rr_yr_2015,rr_yr_2020),
                        "Lower CI"= c(CI_2005[1],CI_2000[1],CI_2015[1],CI_2020[1]),
                        "Upper CI"= c(CI_2005[2],CI_2000[2],CI_2015[2],CI_2020[2]),
                        "p-value" = c(2*(1-pnorm(abs(abs_2005/SE_2005))), 2*(1-pnorm(abs(abs_2000/SE_2000))),
                                      2*(1-pnorm(abs(abs_2015/SE_2015))), 2*(1-pnorm(abs(abs_2020/SE_2020)))))
colnames(yr_spline) <- c("RR","Lower CI","Upper CI","p-value")


rna <- round(min(la_t_pre$log_rna,na.rm = TRUE)):round(max(la_t_pre$log_rna,na.rm = TRUE))
sp_rna <- as.data.frame(ns(rna, df=4))
sp_rna <- cbind(rna,sp_rna) 
sp_rna <- as.data.frame(sp_rna)
rownames(sp_rna) <- sp_rna$rna
sp_rna <- sp_rna[, -1]
colnames(sp_rna) <- c("X1", "X2", "X3", "X4")

uni.rna <- glm(histo_num ~ ns(log_rna,df =4) + offset(log(follow_up)),family = poisson(link =log), data = la_t_pre)
summ_spline  <- summary(uni.rna)
v_spline <- sandwich(uni.rna)
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Coefs", "SE", "Variables")


rr_rna_1 <- rr_rna(1,summ_spline)
rr_rna_3 <- rr_rna(3,summ_spline)
rr_rna_7 <- rr_rna(7,summ_spline)


SE_1 <- spline_se_rna(1)
SE_3 <- spline_se_rna(3)
SE_7 <- spline_se_rna(7)



CI_1 <- CI(rr_rna_1,SE_1)
CI_3 <- CI(rr_rna_3,SE_3)
CI_7 <- CI(rr_rna_7,SE_7)




abs_1 <- abs_val_rna(1,summ_spline)
abs_3 <- abs_val_rna(3,summ_spline)
abs_7 <- abs_val_rna(7,summ_spline)




rna_spline <- data.frame("RR" = c(rr_rna_1,rr_rna_3,rr_rna_7),
                         "Lower CI"= c(CI_1[1],CI_3[1],CI_7[1]),
                         "Upper CI"= c(CI_1[2],CI_3[2],CI_7[2]),
                         "p-value" = c(2*(1-pnorm(abs(abs_1/SE_1))), 2*(1-pnorm(abs(abs_3/SE_3))),
                                       2*(1-pnorm(abs(abs_7/SE_7)))))
colnames(rna_spline) <- c("RR","Lower CI","Upper CI","p-value")



uni.tab1 <- rbind(age_spline ,cd4_spline,rna_spline,yr_spline,uni.fit(uni.gen),uni.fit3(uni.risk),uni.fit3(uni.site),uni.fit(uni.mig))

ms1 <- glm(histo_num ~  ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE)  + as.factor(mig_ende),family = poisson(link = "log"), data = la_t_pre)
#ms2 <- glm(histo_num ~   ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE)  + as.factor(mig_ende),family = poisson(link = "log"), data = la_t_pre)
summ_spline  <- summary(ms1)
v_spline <- sandwich(ms1)
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Coefs", "SE", "Variables")

rr_age_20 <- rr(20,summ_spline)
rr_age_25 <- rr(25,summ_spline)
rr_age_30 <- rr(30,summ_spline)
rr_age_40 <- rr(40,summ_spline)
rr_age_45 <- rr(45,summ_spline)
rr_age_50 <- rr(50,summ_spline)
rr_age_55 <- rr(55,summ_spline)
rr_age_60 <- rr(60,summ_spline)

SE_20 <- spline_se_fun(20)
SE_25 <- spline_se_fun(25)
SE_30 <- spline_se_fun(30)
SE_40 <- spline_se_fun(40)
SE_45 <- spline_se_fun(45)
SE_50 <- spline_se_fun(50)
SE_55 <- spline_se_fun(55)
SE_60 <- spline_se_fun(60)

CI_20 <- CI(rr_age_20,SE_20)
CI_25 <- CI(rr_age_25,SE_25)
CI_30 <- CI(rr_age_30,SE_30)
CI_40 <- CI(rr_age_40,SE_40)
CI_45 <- CI(rr_age_45,SE_45)
CI_50 <- CI(rr_age_50,SE_50)
CI_55 <- CI(rr_age_55,SE_55)
CI_60 <- CI(rr_age_60,SE_60)

abs_20 <- abs_val(20,summ_spline)
abs_25 <- abs_val(25,summ_spline)
abs_30 <- abs_val(30,summ_spline)
abs_40 <- abs_val(40,summ_spline)
abs_45 <- abs_val(45,summ_spline)
abs_50 <- abs_val(50,summ_spline)
abs_55 <- abs_val(55,summ_spline)
abs_60 <- abs_val(60,summ_spline)

rr_sq_5 <- rr_sq_sh(5,summ_spline)
rr_sq_10 <- rr_sq_sh(10,summ_spline)
rr_sq_18.70 <- rr_sq_sh(18.70,summ_spline)
rr_sq_22.40 <- rr_sq_sh(22.40,summ_spline)


SE_sq5 <- spline_se_sh(5)
SE_sq10 <- spline_se_sh(10)
SE_sq18.70 <- spline_se_sh(18.70)
SE_sq22.40 <- spline_se_sh(22.40)


CI_sq5 <- CI(rr_sq_5,SE_sq5)
CI_sq10 <- CI(rr_sq_10,SE_sq10)
CI_sq18.70 <- CI(rr_sq_18.70,SE_sq18.70)
CI_sq22.40 <- CI(rr_sq_22.40,SE_sq22.40)


abs_sq5 <- abs_sh(5,summ_spline)
abs_sq10 <- abs_sh(10,summ_spline)
abs_sq18.70 <- abs_sh(18.70,summ_spline)
abs_sq22.40 <- abs_sh(22.40,summ_spline)



rr_yr_2005 <- rr_yr(2005,summ_spline)
rr_yr_2000 <- rr_yr(2000,summ_spline)
rr_yr_2015 <- rr_yr(2015,summ_spline)
rr_yr_2020 <- rr_yr(2020,summ_spline)

SE_2005 <- spline_se_yr(2005)
SE_2000 <- spline_se_yr(2000)
SE_2015 <- spline_se_yr(2015)
SE_2020 <- spline_se_yr(2020)


CI_2005 <- CI(rr_yr_2005,SE_2005)
CI_2000 <- CI(rr_yr_2000,SE_2000)
CI_2015 <- CI(rr_yr_2015,SE_2015)
CI_2020 <- CI(rr_yr_2020,SE_2020)



abs_2005 <- abs_val_yrs(2005,summ_spline)
abs_2000 <- abs_val_yrs(2000,summ_spline)
abs_2015 <- abs_val_yrs(2015,summ_spline)
abs_2020 <- abs_val_yrs(2020,summ_spline)

rr_rna_1 <- rr_rna(1,summ_spline)
rr_rna_3 <- rr_rna(3,summ_spline)
rr_rna_7 <- rr_rna(7,summ_spline)


SE_1 <- spline_se_rna(1)
SE_3 <- spline_se_rna(3)
SE_7 <- spline_se_rna(7)



CI_1 <- CI(rr_rna_1,SE_1)
CI_3 <- CI(rr_rna_3,SE_3)
CI_7 <- CI(rr_rna_7,SE_7)




abs_1 <- abs_val_rna(1,summ_spline)
abs_3 <- abs_val_rna(3,summ_spline)
abs_7 <- abs_val_rna(7,summ_spline)

age_spl_multi <- data.frame("RR" = c(rr_age_20,rr_age_25,rr_age_30,rr_age_40,rr_age_45,rr_age_50,rr_age_55,rr_age_60),
                            "Lower CI"= c(CI_20[1],CI_25[1],CI_30[1],CI_40[1],CI_45[1],CI_50[1],CI_55[1],CI_60[1]),
                            "Upper CI"= c(CI_20[2],CI_25[2],CI_30[2],CI_40[2],CI_45[2],CI_50[2],CI_55[2],CI_60[2]),
                            "p-value" = c(2*(1-pnorm(abs(abs_20/SE_20))), 2*(1-pnorm(abs(abs_25/SE_25))), 2*(1-pnorm(abs(abs_30/SE_30))),
                                          2*(1-pnorm(abs(abs_40/SE_40))), 2*(1-pnorm(abs(abs_45/SE_45))), 2*(1-pnorm(abs(abs_50/SE_50))),
                                          2*(1-pnorm(abs(abs_55/SE_55))), 2*(1-pnorm(abs(abs_60/SE_60)))))
colnames(age_spl_multi) <- c("RR","Lower CI","Upper CI","p-value")


cd4_spline_multi <- data.frame("RR" = c(rr_sq_5,rr_sq_10,rr_sq_18.70,rr_sq_22.40),
                               "Lower CI"= c(CI_sq5[1],CI_sq10[1],CI_sq18.70[1],CI_sq22.40[1]),
                               "Upper CI"= c(CI_sq5[2],CI_sq10[2],CI_sq18.70[2],CI_sq22.40[2]),
                               "p-value" = c(2*(1-pnorm(abs(abs_sq5/SE_sq5))), 2*(1-pnorm(abs(abs_sq10/SE_sq10))), 2*(1-pnorm(abs(abs_sq18.70/SE_sq18.70))),
                                             2*(1-pnorm(abs(abs_sq22.40/SE_sq22.40)))))

colnames(cd4_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")


yr_spline_multi <- data.frame("RR" = c(rr_yr_2005,rr_yr_2000,rr_yr_2015,rr_yr_2020),
                              "Lower CI"= c(CI_2005[1],CI_2000[1],CI_2015[1],CI_2020[1]),
                              "Upper CI"= c(CI_2005[2],CI_2000[2],CI_2015[2],CI_2020[2]),
                              "p-value" = c(2*(1-pnorm(abs(abs_2005/SE_2005))), 2*(1-pnorm(abs(abs_2000/SE_2000))),
                                            2*(1-pnorm(abs(abs_2015/SE_2015))), 2*(1-pnorm(abs(abs_2020/SE_2020)))))
colnames(yr_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")

rna_spline_multi <- data.frame("RR" = c(rr_rna_1,rr_rna_3,rr_rna_7),
                               "Lower CI"= c(CI_1[1],CI_3[1],CI_7[1]),
                               "Upper CI"= c(CI_1[2],CI_3[2],CI_7[2]),
                               "p-value" = c(2*(1-pnorm(abs(abs_1/SE_1))), 2*(1-pnorm(abs(abs_3/SE_3))),
                                             2*(1-pnorm(abs(abs_7/SE_7)))))
colnames(rna_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")

multi.tab1 <- uni.fit2(ms1)
multi.tab1 <- multi.tab1[-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14),]
multi.tab1 <- rbind(age_spl_multi,cd4_spline_multi,rna_spline_multi,yr_spline_multi,multi.tab1)


tab7 <- cbind(uni.tab1,multi.tab1)
rownames(tab7) <- c("Age 20","Age 25","Age 30","Age 40","Age 45","Age 50","Age 55","Age 60",
                    "CD4 25","CD4 100","CD4 350","CD4 500",
                    "RNA 1","RNA 3","RNA 7", "Year 2000","Year 2005"," Year 2015","Year 2020",
                    "Male ", "Hetero ","IDU ", "Other/Unknown ", "Chile ", "Honduras ","Mexico ", 
                    "Peru ", "Migration status")

save(tab7, file="tab7.Rdata")

### Regression Models ###



## Imputation ##

imp_la_pre <- la_t_pre %>% select(patient_id,age_b,sq_cd4,log_rna, baseline_y,sex,risk1,art_status,mig_ende,follow_up,histo_num,sex,SITE)
imp_la_pre <- mice(imp_la_pre,m=5,printFlag = FALSE, idvars = "patient_id")
set.seed(7)
#la_imp <- aregImpute(~ age_b + sq_cd4+ log_rna + baseline_y+ as.factor(sex)  + as.factor(risk1)   + as.factor(mig_ende) + follow_up, data = la_t_pre, n.impute=20,x=TRUE, nk=0)
#imp.la <- fit.mult.impute(histo_num ~  ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE)  + as.factor(mig_ende), data=la_t_pre,
#family = poisson(link = "log")))

imp_pre <- with(imp_la_pre,
                glm(histo_num ~  ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex) 
                    + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE) + as.factor(mig_ende),
                    family = poisson(link = "log")))

imp.la.age <- with(imp_la_pre,
                   glm(histo_num ~   ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE)  + as.factor(mig_ende), 
                       family = poisson(link = "log")))

a_age <- D1(imp_pre,imp.la.age)
p_age <- round(a_age$result[4],digits=3)

imp.la.cd <- with(imp_la_pre,
                  glm(histo_num ~  ns(age_b,df =4) + ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE)  + as.factor(mig_ende), 
                      family = poisson(link = "log")))

a_cd <- D1(imp_pre,imp.la.cd)
p_cd <- round(a_cd$result[4],digits=3)

imp.la.rna <- with(imp_la_pre,
                   glm(histo_num ~   ns(age_b,df =4) + ns(sq_cd4,df=2) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE)  + as.factor(mig_ende), 
                       family = poisson(link = "log")))

a_rna <- D1(imp_pre,imp.la.rna)
p_rna <- round(a_rna$result[4],digits=3)

imp.la.yr <- with(imp_la_pre,
                  glm(histo_num ~  ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE)  + as.factor(mig_ende), 
                      family = poisson(link = "log")))

a_yr <- D1(imp_pre,imp.la.yr)
p_yr <- round(a_yr$result[4],digits=3)


with(imp_la_pre,
     glm(histo_num ~   ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE)  + as.factor(mig_ende), 
         family = poisson(link = "log")))
a_gen <- D1(imp_pre,imp.la.gen)
p_gen <- round(a_gen$result[4],digits=3)

imp.la.risk <- with(imp_la_pre,
                    glm(histo_num ~ ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up))  + as.factor(SITE)  + as.factor(mig_ende), 
                        family = poisson(link = "log")))

a_risk <- D1(imp_pre,imp.la.risk)
p_risk <- round(a_risk$result[4],digits=3)

imp.la.site <- with(imp_la_pre,
                    glm(histo_num ~   ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1)   + as.factor(mig_ende), 
                        family = poisson(link = "log")))

a_risk <- D1(imp_pre,imp.la.site)
p_risk <- round(a_risk$result[4],digits=3)


a_site <- D1(imp_pre,imp.la.site)
p_site <- round(a_site$result[4],digits=3)

imp.la.mig <- with(imp_la_pre,
                   glm(histo_num ~  ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE),
                       family = poisson(link = "log")))


a_mig <- D1(imp_pre,imp.la.mig)
p_mig <- round(a_mig$result[4],digits=3)


t2 <- pool(imp_pre)
m <- t2$m

## Applying rubin's rule
vcov_list <- list()
for (i in 1:m) {
  vcov_list[[i]] <- vcovHC(imp_pre$analyses[[i]])
} 

vw <- Reduce('+',vcov_list)/m
#vcov_array <- array(unlist(vcov_list), dim = c(dim(vcov_list[[1]]),length(vcov_list)))
## Computing within imputaion variance 
#var_vcov <- apply(vcov_array,c(1,2),var)
#vw <- var_vcov
## Between imputaion variance covariance
qbar <- getqbar(t2)
qhats <- sapply(imp_pre$analyses, coef)
vb <- (1 / (m-1)) * (qhats - qbar) %*% t(qhats - qbar)
vt <- vw + (1 + 1 / (m)) * vb 
vt <- vt[-c(1), -c(1)]

impu.la_pre <- function(model,var) {
  s <-  summary(pool(model))
  s <-  s$estimate[-1]
  sand <- sqrt(diag(var))
  r = matrix(0,nrow = length(s), ncol = 4)
  for (i in 1:(length(s))){
    r[i,] <- c(exp(s[i]), exp(s[i] - 1.96*sand[i]),exp(s[i] + 1.96*sand[i]),round((1 - pnorm(abs(s[i]/sand[i]))) * 2,digits=5))
  }
  
  r <- as.data.frame(r)
  colnames(r) <- c( "RR","Lower CI", "Upper CI","p-value") 
  rownam <- rownames(s)
  rownames(r) <- rownam
  r$ind <- seq_len(nrow(r))
  r <- rbind(r,c(NA,NA,NA,p_gen,ind=14.1),c(NA,NA,NA,p_risk,ind=15.1),c(NA,NA,NA,p_site,ind=18.1),c(NA,NA,NA,p_mig,ind=22.1))
  r <- r %>% arrange(ind) %>% dplyr::select(-ind)
  return(r)
}


impu.la <- impu.la_pre(imp_pre,vt)
cd4s <- seq(min(la_t_pre$sq_cd4,na.rm = T),max(la_t_pre$sq_cd4,na.rm=T), by = 0.01)
sp_cd4 <- as.data.frame(ns(cd4s, df=4))
sp_cd4 <- cbind(cd4s, sp_cd4) 
sp_cd4 <- as.data.frame(sp_cd4)
rownames(sp_cd4) <- sp_cd4$cd4
sp_cd4 <- sp_cd4[, -1]
colnames(sp_cd4) <- c("X1", "X2", "X3", "X4")

ages <- round(min(la_t_pre$age_b)):round(max(la_t_pre$age_b))
sp <- as.data.frame(ns(ages, df=4))
sp <- cbind(ages, sp) 
sp <- as.data.frame(sp)
rownames(sp) <- sp$age
sp <- sp[, -1]
colnames(sp) <- c("X1", "X2", "X3", "X4")

yr <- round(min(la_t_pre$baseline_y,na.rm = TRUE)):round(max(la_t_pre$baseline_y,na.rm = TRUE))
sp_yr <- as.data.frame(ns(yr, df=4))
sp_yr <- cbind(yr, sp_yr) 
sp_yr <- as.data.frame(sp_yr)
rownames(sp_yr) <- sp_yr$yr
sp_yr <- sp_yr[, -1]
colnames(sp_yr) <- c("X1", "X2", "X3", "X4")

rna <- round(min(la_t_pre$log_rna,na.rm = TRUE)):round(max(la_t_pre$log_rna,na.rm = TRUE))
sp_rna <- as.data.frame(ns(rna, df=4))
sp_rna <- cbind(rna,sp_rna) 
sp_rna <- as.data.frame(sp_rna)
rownames(sp_rna) <- sp_rna$rna
sp_rna <- sp_rna[, -1]
colnames(sp_rna) <- c("X1", "X2", "X3", "X4")

summ_spline  <- summary(pool(imp_pre))
v_spline <- vt
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- summ_spline[-1,c(1,2)]
summ_spline <- as.data.frame(cbind(summ_spline, se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Variables","Coefs", "SE" )

rr_age_20 <- rr(20,summ_spline)
rr_age_25 <- rr(25,summ_spline)
rr_age_30 <- rr(30,summ_spline)
rr_age_40 <- rr(40,summ_spline)
rr_age_45 <- rr(45,summ_spline)
rr_age_50 <- rr(50,summ_spline)
rr_age_55 <- rr(55,summ_spline)
rr_age_60 <- rr(60,summ_spline)

SE_20 <- spline_se_fun(20)
SE_25 <- spline_se_fun(25)
SE_30 <- spline_se_fun(30)
SE_40 <- spline_se_fun(40)
SE_45 <- spline_se_fun(45)
SE_50 <- spline_se_fun(50)
SE_55 <- spline_se_fun(55)
SE_60 <- spline_se_fun(60)

CI_20 <- CI(rr_age_20,SE_20)
CI_25 <- CI(rr_age_25,SE_25)
CI_30 <- CI(rr_age_30,SE_30)
CI_40 <- CI(rr_age_40,SE_40)
CI_45 <- CI(rr_age_45,SE_45)
CI_50 <- CI(rr_age_50,SE_50)
CI_55 <- CI(rr_age_55,SE_55)
CI_60 <- CI(rr_age_60,SE_60)

abs_20 <- abs_val(20,summ_spline)
abs_25 <- abs_val(25,summ_spline)
abs_30 <- abs_val(30,summ_spline)
abs_40 <- abs_val(40,summ_spline)
abs_45 <- abs_val(45,summ_spline)
abs_50 <- abs_val(50,summ_spline)
abs_55 <- abs_val(55,summ_spline)
abs_60 <- abs_val(60,summ_spline)

rr_sq_5 <- rr_sq_sh(5,summ_spline)
rr_sq_10 <- rr_sq_sh(10,summ_spline)
rr_sq_18.70 <- rr_sq_sh(18.70,summ_spline)
rr_sq_22.40 <- rr_sq_sh(22.40,summ_spline)


SE_sq5 <- spline_se_sh(5)
SE_sq10 <- spline_se_sh(10)
SE_sq18.70 <- spline_se_sh(18.70)
SE_sq22.40 <- spline_se_sh(22.40)

CI_sq5 <- CI(rr_sq_5,SE_5)
CI_sq10 <- CI(rr_sq_10,SE_10)
CI_sq18.70 <- CI(rr_sq_18.70,SE_18.70)
CI_sq22.40 <- CI(rr_sq_22.40,SE_22.40)

abs_sq5 <- abs_sh(5,summ_spline)
abs_sq10 <- abs_sh(10,summ_spline)
abs_sq18.70 <- abs_sh(18.70,summ_spline)
abs_sq22.40 <- abs_sh(22.40,summ_spline)


rr_yr_2005 <- rr_yr(2005,summ_spline)
rr_yr_2000 <- rr_yr(2000,summ_spline)
rr_yr_2015 <- rr_yr(2015,summ_spline)
rr_yr_2020 <- rr_yr(2020,summ_spline)

SE_2005 <- spline_se_yr(2005)
SE_2000 <- spline_se_yr(2000)
SE_2015 <- spline_se_yr(2015)
SE_2020 <- spline_se_yr(2020)


CI_2005 <- CI(rr_yr_2005,SE_2005)
CI_2000 <- CI(rr_yr_2000,SE_2000)
CI_2015 <- CI(rr_yr_2015,SE_2015)
CI_2020 <- CI(rr_yr_2020,SE_2020)



abs_2005 <- abs_val_yrs(2005,summ_spline)
abs_2000 <- abs_val_yrs(2000,summ_spline)
abs_2015 <- abs_val_yrs(2015,summ_spline)
abs_2020 <- abs_val_yrs(2020,summ_spline)

rr_rna_1 <- rr_rna(1,summ_spline)
rr_rna_3 <- rr_rna(3,summ_spline)
rr_rna_7 <- rr_rna(7,summ_spline)


SE_1 <- spline_se_rna(1)
SE_3 <- spline_se_rna(3)
SE_7 <- spline_se_rna(7)



CI_1 <- CI(rr_rna_1,SE_1)
CI_3 <- CI(rr_rna_3,SE_3)
CI_7 <- CI(rr_rna_7,SE_7)




abs_1 <- abs_val_rna(1,summ_spline)
abs_3 <- abs_val_rna(3,summ_spline)
abs_7 <- abs_val_rna(7,summ_spline)

age_spl_multi <- data.frame("RR" = c(rr_age_20,rr_age_25,rr_age_30,rr_age_40,rr_age_45,rr_age_50,rr_age_55,rr_age_60),
                            "Lower CI"= c(CI_20[1],CI_25[1],CI_30[1],CI_40[1],CI_45[1],CI_50[1],CI_55[1],CI_60[1]),
                            "Upper CI"= c(CI_20[2],CI_25[2],CI_30[2],CI_40[2],CI_45[2],CI_50[2],CI_55[2],CI_60[2]),
                            "p-value" = c(2*(1-pnorm(abs(abs_20/SE_20))), 2*(1-pnorm(abs(abs_25/SE_25))), 2*(1-pnorm(abs(abs_30/SE_30))),
                                          2*(1-pnorm(abs(abs_40/SE_40))), 2*(1-pnorm(abs(abs_45/SE_45))), 2*(1-pnorm(abs(abs_50/SE_50))),
                                          2*(1-pnorm(abs(abs_55/SE_55))), 2*(1-pnorm(abs(abs_60/SE_60)))))
colnames(age_spl_multi) <- c("RR","Lower CI","Upper CI","p-value")


cd4_spline_multi <- data.frame("RR" = c(rr_sq_5,rr_sq_10,rr_sq_18.70,rr_sq_22.40),
                               "Lower CI"= c(CI_sq5[1],CI_sq10[1],CI_sq18.70[1],CI_sq22.40[1]),
                               "Upper CI"= c(CI_sq5[2],CI_sq10[2],CI_sq18.70[2],CI_sq22.40[2]),
                               "p-value" = c(2*(1-pnorm(abs(abs_sq5/SE_sq5))), 2*(1-pnorm(abs(abs_sq10/SE_sq10))), 2*(1-pnorm(abs(abs_sq18.70/SE_sq18.70))),
                                             2*(1-pnorm(abs(abs_sq22.40/SE_sq22.40)))))
colnames(cd4_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")


yr_spline_multi <- data.frame("RR" = c(rr_yr_2005,rr_yr_2000,rr_yr_2015,rr_yr_2020),
                              "Lower CI"= c(CI_2005[1],CI_2000[1],CI_2015[1],CI_2020[1]),
                              "Upper CI"= c(CI_2005[2],CI_2000[2],CI_2015[2],CI_2020[2]),
                              "p-value" = c(2*(1-pnorm(abs(abs_2005/SE_2005))), 2*(1-pnorm(abs(abs_2000/SE_2000))),
                                            2*(1-pnorm(abs(abs_2015/SE_2015))), 2*(1-pnorm(abs(abs_2020/SE_2020)))))
colnames(yr_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")

rna_spline_multi <- data.frame("RR" = c(rr_rna_1,rr_rna_3,rr_rna_7),
                               "Lower CI"= c(CI_1[1],CI_3[1],CI_7[1]),
                               "Upper CI"= c(CI_1[2],CI_3[2],CI_7[2]),
                               "p-value" = c(2*(1-pnorm(abs(abs_1/SE_1))), 2*(1-pnorm(abs(abs_3/SE_3))),
                                             2*(1-pnorm(abs(abs_7/SE_7)))))
colnames(rna_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")

impu.la <- impu.la[-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14),]
tab9 <- rbind(c(NA,NA,NA,p_age),age_spl_multi,c(NA,NA,NA,p_cd),cd4_spline_multi,c(NA,NA,NA,p_rna),rna_spline_multi,c(NA,NA,NA,p_yr),yr_spline_multi,impu.la)
rownames(tab9) <- c("Age (ref=35)","20","25","30","40","45","50","55","60",
                    "CD4 (ref=200)","25 ","100","350","500",
                    "log HIV RNA(ref=5)","1","3 " ,"7","Year (ref=2010)", "2000","2005","2015","2020","Sex",
                    "Male ", "Risk factor (ref =MSM)","Hetero ","IDU ", "Other/Unknown ", "Country (ref=Brazil)", "Chile ", "Honduras ","Mexico ", "Peru ",
                    "Migration status (ref=No)","Yes")

save(tab9, file="tab9.Rdata")


# Relevel
la_t_post$sex <- relevel(as.factor(la_t_post$sex), ref = "Female")
la_t_post$SITE <- relevel(as.factor(la_t_post$SITE), ref = "Brazil")
la_t_post$risk1 <- relevel(la_t_post$risk1, ref = "MSM")
la_t_post$art_status <- relevel(as.factor(la_t_post$art_status), ref = "post-ART initiation")
la_t_post$mig_ende <- relevel(as.factor(la_t_post$mig_ende), ref = "No")

## Non spline univariate models
uni0 <- glm(histo_num ~ offset(log(follow_up)),family = poisson(link = "log"), data = la_t_post)
uni.gen <- glm(histo_num ~as.factor(sex) + offset(log(follow_up)),family = poisson(link = "log"), data = la_t_post)
uni.risk <- glm(histo_num ~as.factor(risk1) + offset(log(follow_up)),family = poisson(link = "log"), data = la_t_post)
uni.site <- glm(histo_num ~as.factor(SITE) + offset(log(follow_up)),family = poisson(link = "log"), data = la_t_post)
uni.mig <- glm(histo_num ~as.factor(mig_ende) + offset(log(follow_up)),family = poisson(link = "log"), data = la_t_post)

## Splines ##
ages <- round(min(la_t_post$age_b)):round(max(la_t_post$age_b))
sp <- as.data.frame(ns(ages, df=4))
sp <- cbind(ages, sp) 
sp <- as.data.frame(sp)
rownames(sp) <- sp$age
sp <- sp[, -1]
colnames(sp) <- c("X1", "X2", "X3", "X4")

uni.age <- glm(histo_num ~ ns(age_b,df =4) + offset(log(follow_up)),family = poisson(link =log), data = la_t_post)
summ_spline  <- summary(uni.age)
v_spline <- sandwich(uni.age)
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Coefs", "SE", "Variables")

rr_age_20 <- rr(20,summ_spline)
rr_age_25 <- rr(25,summ_spline)
rr_age_30 <- rr(30,summ_spline)
rr_age_40 <- rr(40,summ_spline)
rr_age_45 <- rr(45,summ_spline)
rr_age_50 <- rr(50,summ_spline)
rr_age_55 <- rr(55,summ_spline)
rr_age_60 <- rr(60,summ_spline)


SE_20 <- spline_se_fun(20)
SE_25 <- spline_se_fun(25)
SE_30 <- spline_se_fun(30)
SE_40 <- spline_se_fun(40)
SE_45 <- spline_se_fun(45)
SE_50 <- spline_se_fun(50)
SE_55 <- spline_se_fun(55)
SE_60 <- spline_se_fun(60)

CI <- function(rr,SE) 
{ exp(log(rr) + c(-1,1)*1.96*SE)}


CI_20 <- CI(rr_age_20,SE_20)
CI_25 <- CI(rr_age_25,SE_25)
CI_30 <- CI(rr_age_30,SE_30)
CI_40 <- CI(rr_age_40,SE_40)
CI_45 <- CI(rr_age_45,SE_45)
CI_50 <- CI(rr_age_50,SE_50)
CI_55 <- CI(rr_age_55,SE_55)
CI_60 <- CI(rr_age_60,SE_60)

abs_20 <- abs_val(20,summ_spline)
abs_25 <- abs_val(25,summ_spline)
abs_30 <- abs_val(30,summ_spline)
abs_40 <- abs_val(40,summ_spline)
abs_45 <- abs_val(45,summ_spline)
abs_50 <- abs_val(50,summ_spline)
abs_55 <- abs_val(55,summ_spline)
abs_60 <- abs_val(60,summ_spline)

age_spline <- data.frame("RR" = c(rr_age_20,rr_age_25,rr_age_30,rr_age_40,rr_age_45,rr_age_50,rr_age_55,rr_age_60),
                         "Lower CI"= c(CI_20[1],CI_25[1],CI_30[1],CI_40[1],CI_45[1],CI_50[1],CI_55[1],CI_60[1]),
                         "Upper CI"= c(CI_20[2],CI_25[2],CI_30[2],CI_40[2],CI_45[2],CI_50[2],CI_55[2],CI_60[2]),
                         "p-value" = c(2*(1-pnorm(abs(abs_20/SE_20))), 2*(1-pnorm(abs(abs_25/SE_25))), 2*(1-pnorm(abs(abs_30/SE_30))),
                                       2*(1-pnorm(abs(abs_40/SE_40))), 2*(1-pnorm(abs(abs_45/SE_45))), 2*(1-pnorm(abs(abs_50/SE_50))),
                                       2*(1-pnorm(abs(abs_55/SE_55))), 2*(1-pnorm(abs(abs_60/SE_60)))))
colnames(age_spline) <- c("RR","Lower CI","Upper CI","p-value")


cd4s <- seq(min(la_t_post$sq_cd4,na.rm = T),max(la_t_post$sq_cd4,na.rm=T), by = 0.01)
sp_cd4 <- as.data.frame(ns(cd4s, df=4))
sp_cd4 <- cbind(cd4s, sp_cd4) 
sp_cd4 <- as.data.frame(sp_cd4)
rownames(sp_cd4) <- sp_cd4$cd4
sp_cd4 <- sp_cd4[, -1]
colnames(sp_cd4) <- c("X1", "X2", "X3", "X4")

uni.cd <- glm(histo_num ~ ns(sq_cd4,df =4) + offset(log(follow_up)),family = poisson(link =log), data = la_t_post)
summ_spline  <- summary(uni.cd)
v_spline <-sandwich(uni.cd)
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Coefs", "SE", "Variables")

rr_sq_5 <- rr_sq(5,summ_spline)
rr_sq_10 <- rr_sq(10,summ_spline)
rr_sq_18.70 <- rr_sq(18.70,summ_spline)
rr_sq_22.40 <- rr_sq(22.40,summ_spline)


SE_5 <- spline_se_fun_sq_cd(5)
SE_10 <- spline_se_fun_sq_cd(10)
SE_18.70 <- spline_se_fun_sq_cd(18.70)
SE_22.40 <- spline_se_fun_sq_cd(22.40)



CI_5 <- CI(rr_sq_5,SE_5)
CI_10 <- CI(rr_sq_10,SE_10)
CI_18.70 <- CI(rr_sq_18.70,SE_18.70)
CI_22.40 <- CI(rr_sq_22.40,SE_22.40)



abs_5 <- abs_sq_cd(5,summ_spline)
abs_10 <- abs_sq_cd(10,summ_spline)
abs_18.70 <- abs_sq_cd(18.70,summ_spline)
abs_22.40 <- abs_sq_cd(22.40,summ_spline)



cd4_spline <- data.frame("RR" = c(rr_sq_5,rr_sq_10,rr_sq_18.70,rr_sq_22.40),
                         "Lower CI"= c(CI_5[1],CI_10[1],CI_18.70[1],CI_20[1]),
                         "Upper CI"= c(CI_5[2],CI_10[2],CI_18.70[2],CI_20[2]),
                         "p-value" = c(2*(1-pnorm(abs(abs_5/SE_5))), 2*(1-pnorm(abs(abs_10/SE_10))), 2*(1-pnorm(abs(abs_18.70/SE_18.70))),
                                       2*(1-pnorm(abs(abs_25/SE_22.40)))))
colnames(cd4_spline) <- c("RR","Lower CI","Upper CI","p-value")

yr <- round(min(la_t_post$baseline_y,na.rm = TRUE)):round(max(la_t_post$baseline_y,na.rm = TRUE))
sp_yr <- as.data.frame(ns(yr, df=4))
sp_yr <- cbind(yr, sp_yr) 
sp_yr <- as.data.frame(sp_yr)
rownames(sp_yr) <- sp_yr$yr
sp_yr <- sp_yr[, -1]
colnames(sp_yr) <- c("X1", "X2", "X3", "X4")

uni.yrs <- glm(histo_num ~ ns(baseline_y,df =4) + offset(log(follow_up)),family = poisson(link =log), data = la_t_post)
summ_spline  <- summary(uni.yrs)
v_spline <- sandwich(uni.yrs)
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Coefs", "SE", "Variables")


rr_yr_2005 <- rr_yr(2005,summ_spline)
rr_yr_2000 <- rr_yr(2000,summ_spline)
rr_yr_2015 <- rr_yr(2015,summ_spline)
rr_yr_2020 <- rr_yr(2020,summ_spline)

SE_2005 <- spline_se_yr(2005)
SE_2000 <- spline_se_yr(2000)
SE_2015 <- spline_se_yr(2015)
SE_2020 <- spline_se_yr(2020)


CI_2005 <- CI(rr_yr_2005,SE_2005)
CI_2000 <- CI(rr_yr_2000,SE_2000)
CI_2015 <- CI(rr_yr_2015,SE_2015)
CI_2020 <- CI(rr_yr_2020,SE_2020)



abs_2005 <- abs_val_yrs(2005,summ_spline)
abs_2000 <- abs_val_yrs(2000,summ_spline)
abs_2015 <- abs_val_yrs(2015,summ_spline)
abs_2020 <- abs_val_yrs(2020,summ_spline)



yr_spline <- data.frame("RR" = c(rr_yr_2005,rr_yr_2000,rr_yr_2015,rr_yr_2020),
                        "Lower CI"= c(CI_2005[1],CI_2000[1],CI_2015[1],CI_2020[1]),
                        "Upper CI"= c(CI_2005[2],CI_2000[2],CI_2015[2],CI_2020[2]),
                        "p-value" = c(2*(1-pnorm(abs(abs_2005/SE_2005))), 2*(1-pnorm(abs(abs_2000/SE_2000))),
                                      2*(1-pnorm(abs(abs_2015/SE_2015))), 2*(1-pnorm(abs(abs_2020/SE_2020)))))
colnames(yr_spline) <- c("RR","Lower CI","Upper CI","p-value")


rna <- round(min(la_t_post$log_rna,na.rm = TRUE)):round(max(la_t_post$log_rna,na.rm = TRUE))
sp_rna <- as.data.frame(ns(rna, df=4))
sp_rna <- cbind(rna,sp_rna) 
sp_rna <- as.data.frame(sp_rna)
rownames(sp_rna) <- sp_rna$rna
sp_rna <- sp_rna[, -1]
colnames(sp_rna) <- c("X1", "X2", "X3", "X4")

uni.rna <- glm(histo_num ~ ns(log_rna,df =4) + offset(log(follow_up)),family = poisson(link =log), data = la_t_post)
summ_spline  <- summary(uni.rna)
v_spline <- sandwich(uni.rna)
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Coefs", "SE", "Variables")


rr_rna_1 <- rr_rna(1,summ_spline)
rr_rna_3 <- rr_rna(3,summ_spline)
rr_rna_7 <- rr_rna(7,summ_spline)


SE_1 <- spline_se_rna(1)
SE_3 <- spline_se_rna(3)
SE_7 <- spline_se_rna(7)



CI_1 <- CI(rr_rna_1,SE_1)
CI_3 <- CI(rr_rna_3,SE_3)
CI_7 <- CI(rr_rna_7,SE_7)




abs_1 <- abs_val_rna(1,summ_spline)
abs_3 <- abs_val_rna(3,summ_spline)
abs_7 <- abs_val_rna(7,summ_spline)




rna_spline <- data.frame("RR" = c(rr_rna_1,rr_rna_3,rr_rna_7),
                         "Lower CI"= c(CI_1[1],CI_3[1],CI_7[1]),
                         "Upper CI"= c(CI_1[2],CI_3[2],CI_7[2]),
                         "p-value" = c(2*(1-pnorm(abs(abs_1/SE_1))), 2*(1-pnorm(abs(abs_3/SE_3))),
                                       2*(1-pnorm(abs(abs_7/SE_7)))))
colnames(rna_spline) <- c("RR","Lower CI","Upper CI","p-value")



uni.tab1 <- rbind(age_spline ,cd4_spline,rna_spline,yr_spline,uni.fit(uni.gen),uni.fit3(uni.risk),uni.fit3(uni.site),uni.fit(uni.mig))

ms1 <- glm(histo_num ~  ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE)  + as.factor(mig_ende),family = poisson(link = "log"), data = la_t_post)
#ms2 <- glm(histo_num ~   ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE)  + as.factor(mig_ende),family = poisson(link = "log"), data = la_t_post)
summ_spline  <- summary(ms1)
v_spline <- sandwich(ms1)
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- as.data.frame(cbind(summ_spline$coefficients[,1], se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Coefs", "SE", "Variables")

rr_age_20 <- rr(20,summ_spline)
rr_age_25 <- rr(25,summ_spline)
rr_age_30 <- rr(30,summ_spline)
rr_age_40 <- rr(40,summ_spline)
rr_age_45 <- rr(45,summ_spline)
rr_age_50 <- rr(50,summ_spline)
rr_age_55 <- rr(55,summ_spline)
rr_age_60 <- rr(60,summ_spline)

SE_20 <- spline_se_fun(20)
SE_25 <- spline_se_fun(25)
SE_30 <- spline_se_fun(30)
SE_40 <- spline_se_fun(40)
SE_45 <- spline_se_fun(45)
SE_50 <- spline_se_fun(50)
SE_55 <- spline_se_fun(55)
SE_60 <- spline_se_fun(60)

CI_20 <- CI(rr_age_20,SE_20)
CI_25 <- CI(rr_age_25,SE_25)
CI_30 <- CI(rr_age_30,SE_30)
CI_40 <- CI(rr_age_40,SE_40)
CI_45 <- CI(rr_age_45,SE_45)
CI_50 <- CI(rr_age_50,SE_50)
CI_55 <- CI(rr_age_55,SE_55)
CI_60 <- CI(rr_age_60,SE_60)

abs_20 <- abs_val(20,summ_spline)
abs_25 <- abs_val(25,summ_spline)
abs_30 <- abs_val(30,summ_spline)
abs_40 <- abs_val(40,summ_spline)
abs_45 <- abs_val(45,summ_spline)
abs_50 <- abs_val(50,summ_spline)
abs_55 <- abs_val(55,summ_spline)
abs_60 <- abs_val(60,summ_spline)

rr_sq_5 <- rr_sq_sh(5,summ_spline)
rr_sq_10 <- rr_sq_sh(10,summ_spline)
rr_sq_18.70 <- rr_sq_sh(18.70,summ_spline)
rr_sq_22.40 <- rr_sq_sh(22.40,summ_spline)


SE_sq5 <- spline_se_sh(5)
SE_sq10 <- spline_se_sh(10)
SE_sq18.70 <- spline_se_sh(18.70)
SE_sq22.40 <- spline_se_sh(22.40)


CI_sq5 <- CI(rr_sq_5,SE_sq5)
CI_sq10 <- CI(rr_sq_10,SE_sq10)
CI_sq18.70 <- CI(rr_sq_18.70,SE_sq18.70)
CI_sq22.40 <- CI(rr_sq_22.40,SE_sq22.40)


abs_sq5 <- abs_sh(5,summ_spline)
abs_sq10 <- abs_sh(10,summ_spline)
abs_sq18.70 <- abs_sh(18.70,summ_spline)
abs_sq22.40 <- abs_sh(22.40,summ_spline)



rr_yr_2005 <- rr_yr(2005,summ_spline)
rr_yr_2000 <- rr_yr(2000,summ_spline)
rr_yr_2015 <- rr_yr(2015,summ_spline)
rr_yr_2020 <- rr_yr(2020,summ_spline)

SE_2005 <- spline_se_yr(2005)
SE_2000 <- spline_se_yr(2000)
SE_2015 <- spline_se_yr(2015)
SE_2020 <- spline_se_yr(2020)


CI_2005 <- CI(rr_yr_2005,SE_2005)
CI_2000 <- CI(rr_yr_2000,SE_2000)
CI_2015 <- CI(rr_yr_2015,SE_2015)
CI_2020 <- CI(rr_yr_2020,SE_2020)



abs_2005 <- abs_val_yrs(2005,summ_spline)
abs_2000 <- abs_val_yrs(2000,summ_spline)
abs_2015 <- abs_val_yrs(2015,summ_spline)
abs_2020 <- abs_val_yrs(2020,summ_spline)

rr_rna_1 <- rr_rna(1,summ_spline)
rr_rna_3 <- rr_rna(3,summ_spline)
rr_rna_7 <- rr_rna(7,summ_spline)


SE_1 <- spline_se_rna(1)
SE_3 <- spline_se_rna(3)
SE_7 <- spline_se_rna(7)



CI_1 <- CI(rr_rna_1,SE_1)
CI_3 <- CI(rr_rna_3,SE_3)
CI_7 <- CI(rr_rna_7,SE_7)




abs_1 <- abs_val_rna(1,summ_spline)
abs_3 <- abs_val_rna(3,summ_spline)
abs_7 <- abs_val_rna(7,summ_spline)

age_spl_multi <- data.frame("RR" = c(rr_age_20,rr_age_25,rr_age_30,rr_age_40,rr_age_45,rr_age_50,rr_age_55,rr_age_60),
                            "Lower CI"= c(CI_20[1],CI_25[1],CI_30[1],CI_40[1],CI_45[1],CI_50[1],CI_55[1],CI_60[1]),
                            "Upper CI"= c(CI_20[2],CI_25[2],CI_30[2],CI_40[2],CI_45[2],CI_50[2],CI_55[2],CI_60[2]),
                            "p-value" = c(2*(1-pnorm(abs(abs_20/SE_20))), 2*(1-pnorm(abs(abs_25/SE_25))), 2*(1-pnorm(abs(abs_30/SE_30))),
                                          2*(1-pnorm(abs(abs_40/SE_40))), 2*(1-pnorm(abs(abs_45/SE_45))), 2*(1-pnorm(abs(abs_50/SE_50))),
                                          2*(1-pnorm(abs(abs_55/SE_55))), 2*(1-pnorm(abs(abs_60/SE_60)))))
colnames(age_spl_multi) <- c("RR","Lower CI","Upper CI","p-value")


cd4_spline_multi <- data.frame("RR" = c(rr_sq_5,rr_sq_10,rr_sq_18.70,rr_sq_22.40),
                               "Lower CI"= c(CI_sq5[1],CI_sq10[1],CI_sq18.70[1],CI_sq22.40[1]),
                               "Upper CI"= c(CI_sq5[2],CI_sq10[2],CI_sq18.70[2],CI_sq22.40[2]),
                               "p-value" = c(2*(1-pnorm(abs(abs_sq5/SE_sq5))), 2*(1-pnorm(abs(abs_sq10/SE_sq10))), 2*(1-pnorm(abs(abs_sq18.70/SE_sq18.70))),
                                             2*(1-pnorm(abs(abs_sq22.40/SE_sq22.40)))))

colnames(cd4_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")


yr_spline_multi <- data.frame("RR" = c(rr_yr_2005,rr_yr_2000,rr_yr_2015,rr_yr_2020),
                              "Lower CI"= c(CI_2005[1],CI_2000[1],CI_2015[1],CI_2020[1]),
                              "Upper CI"= c(CI_2005[2],CI_2000[2],CI_2015[2],CI_2020[2]),
                              "p-value" = c(2*(1-pnorm(abs(abs_2005/SE_2005))), 2*(1-pnorm(abs(abs_2000/SE_2000))),
                                            2*(1-pnorm(abs(abs_2015/SE_2015))), 2*(1-pnorm(abs(abs_2020/SE_2020)))))
colnames(yr_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")

rna_spline_multi <- data.frame("RR" = c(rr_rna_1,rr_rna_3,rr_rna_7),
                               "Lower CI"= c(CI_1[1],CI_3[1],CI_7[1]),
                               "Upper CI"= c(CI_1[2],CI_3[2],CI_7[2]),
                               "p-value" = c(2*(1-pnorm(abs(abs_1/SE_1))), 2*(1-pnorm(abs(abs_3/SE_3))),
                                             2*(1-pnorm(abs(abs_7/SE_7)))))
colnames(rna_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")

multi.tab1 <- uni.fit2(ms1)
multi.tab1 <- multi.tab1[-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14),]
multi.tab1 <- rbind(age_spl_multi,cd4_spline_multi,rna_spline_multi,yr_spline_multi,multi.tab1)


tab11 <- cbind(uni.tab1,multi.tab1)
rownames(tab11) <- c("Age 20","Age 25","Age 30","Age 40","Age 45","Age 50","Age 55","Age 60",
                     "CD4 25","CD4 100","CD4 350","CD4 500",
                     "RNA 1","RNA 3","RNA 7", "Year 2000","Year 2005"," Year 2015","Year 2020",
                     "Male ", "Hetero ","IDU ", "Other/Unknown ", "Chile ", "Honduras ","Mexico ", 
                     "Peru ","Migration status")

save(tab11, file="tab11.Rdata")

### Regression Models ###

## removing transgender cate. because of error 


# Relevel
tn_t_post$birthsex <- relevel(as.factor(tn_t_post$birthsex), ref = "Female")
tn_t_post$risk1 <- relevel(tn_t_post$risk1, ref = "MSM")
tn_t_post$mig_ende <- relevel(as.factor(tn_t_post$mig_ende), ref = "No")

## Non spline univariate models
uni.gen <- glm(histo_num ~as.factor(birthsex) + offset(log(follow_up)),family = poisson(link = "log"), data = tn_t_post)
uni.risk <- glm(histo_num ~as.factor(risk1) + offset(log(follow_up)),family = poisson(link = "log"), data = tn_t_post)
uni.mig <- glm(histo_num ~as.factor(mig_ende) + offset(log(follow_up)),family = poisson(link = "log"), data = tn_t_post)
uni.age <- glm(histo_num ~ age_b + offset(log(follow_up)),family = poisson(link =log), data = tn_t_post)
uni.cd <- glm(histo_num ~ sq_cd4 + offset(log(follow_up)),family = poisson(link =log), data = tn_t_post)
uni.yrs <- glm(histo_num ~ baseline_y + offset(log(follow_up)),family = poisson(link =log), data = tn_t_post)
uni.rna <- glm(histo_num ~ log_rna + offset(log(follow_up)),family = poisson(link =log), data = tn_t_post)
uni.tab1 <- rbind(uni.fit.conti(uni.age,10),uni.fit.conti(uni.cd,50),uni.fit.conti(uni.yrs,10),uni.fit.conti(uni.rna,4),uni.fit(uni.gen),uni.fit3(uni.risk),uni.fit(uni.mig))
ms1 <- glm(histo_num ~  age_b + sq_cd4 + log_rna + baseline_y + as.factor(birthsex)  + offset(log(follow_up)) + as.factor(risk1)   + as.factor(mig_ende),family = poisson(link = "log"), data = tn_t_post)
multi.tab1 <- uni.fit2(ms1)
multi.tab1 <- multi.tab1[-c(1,2,3,4),]
multi.tab1 <- rbind(multi.fit.conti2(ms1),multi.tab1)


tab12 <- cbind(uni.tab1,multi.tab1)
rownames(tab12) <- c("Age 35 (ref=25) ",
                     "CD4 100 (ref =50) ",
                     "RNA 5 (ref=1) ","Year 2015 (ref=2005)",
                     "Male ", "Hetero ","IDU ", "Other/Unknown ","Migration status")

save(tab12, file="tab12.Rdata")

## Imputation ##

imp_post <- la_t_post %>% select(patient_id,age_b,sq_cd4,log_rna, baseline_y,sex,risk1,art_status,mig_ende,follow_up,histo_num,sex,SITE)
imp_post <- mice(imp_post,m=5,printFlag = FALSE, idvars = "patient_id")

set.seed(7)
#la_imp <- aregImpute(~ age_b + sq_cd4+ log_rna + baseline_y+ as.factor(sex)  + as.factor(risk1)   + as.factor(mig_ende) + follow_up, data = la_t_post, n.impute=20,x=TRUE, nk=0)

imp_la_post <- with(imp_post,
                    glm(histo_num ~  ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE)  + as.factor(mig_ende), 
                        family = poisson(link = "log")))


imp.la.age <- with(imp_post,
                   glm(histo_num ~  ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE)  + as.factor(mig_ende),
                       family = poisson(link = "log")))

a_age <- D1(imp_la_post,imp.la.age)
p_age <- round(a_age$result[4],digits=3)

imp.la.cd <- with(imp_post,
                  glm(histo_num ~ ns(age_b,df =4) + ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE)  + as.factor(mig_ende), 
                      family = poisson(link = "log")))

a_cd <- D1(imp_la_post,imp.la.cd)
p_cd <- round(a_cd$result[4],digits=3)

imp.la.rna <- with(imp_post,
                   glm(histo_num ~ ns(age_b,df =4) + ns(sq_cd4,df=2) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE)  + as.factor(mig_ende), 
                       family = poisson(link = "log")))

a_rna <- D1(imp_la_post,imp.la.rna)
p_rna <- round(a_rna$result[4],digits=3)

imp.la.yr <- with(imp_post,
                  glm(histo_num ~  ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE)  + as.factor(mig_ende), 
                      family = poisson(link = "log")))

a_yr <- D1(imp_la_post,imp.la.yr)
p_yr <- round(a_yr$result[4],digits=3)


imp.la.gen <- with(imp_post,
                   glm(histo_num ~ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE)  + as.factor(mig_ende), 
                       family = poisson(link = "log")))

a_gen <- D1(imp_la_post,imp.la.gen)
p_gen <- round(a_gen$result[4],digits=3)

imp.la.risk <- with(imp_post,
                    glm(histo_num ~ ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up))  + as.factor(SITE)  + as.factor(mig_ende), 
                        family = poisson(link = "log")))

a_risk <- D1(imp_la_post,imp.la.risk)
p_risk <- round(a_risk$result[4],digits=3)

imp.la.site <- with(imp_post,
                    glm(histo_num ~  ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1)   + as.factor(mig_ende), 
                        family = poisson(link = "log")))

a_site <- D1(imp_la_post,imp.la.site)
p_site <- round(a_site$result[4],digits=3)

imp.la.mig <- with(imp_post,
                   glm(histo_num ~ ns(age_b,df =4) + ns(sq_cd4,df=2)+ ns(log_rna,df =4) + ns(baseline_y,df=4)+ as.factor(sex)  + offset(log(follow_up)) + as.factor(risk1) + as.factor(SITE) , 
                       family = poisson(link = "log")))


a_mig <- D1(imp_la_post,imp.la.mig)
p_mig <- round(a_mig$result[4],digits=3)

## Variance covariance
t2 <- pool(imp_la_post)
m <- t2$m

## Applying rubin's rule
vcov_list <- list()
for (i in 1:m) {
  vcov_list[[i]] <- vcovHC(imp_la_post$analyses[[i]])
} 

vw <- Reduce('+',vcov_list)/m
#vcov_array <- array(unlist(vcov_list), dim = c(dim(vcov_list[[1]]),length(vcov_list)))
## Computing within imputaion variance 
#var_vcov <- apply(vcov_array,c(1,2),var)
#vw <- var_vcov
## Between imputaion variance covariance
qbar <- getqbar(t2)
qhats <- sapply(imp_la_post$analyses, coef)
vb <- (1 / (m-1)) * (qhats - qbar) %*% t(qhats - qbar)
vt <- vw + (1 + 1 / (m)) * vb 
vt <- vt[-c(1), -c(1)]

impu.la.post <- function(model,var) {
  s <-  summary(pool(model))
  s <-  s$estimate[-1]
  sand <- sqrt(diag(var))
  r = matrix(0,nrow = length(s), ncol = 4)
  for (i in 1:(length(s))){
    r[i,] <- c(exp(s[i]), exp(s[i] - 1.96*sand[i]),exp(s[i] + 1.96*sand[i]),round((1 - pnorm(abs(s[i]/sand[i]))) * 2,digits=5))
  }
  
  r <- as.data.frame(r)
  colnames(r) <- c( "RR","Lower CI", "Upper CI","p-value") 
  rownam <- rownames(s)
  rownames(r) <- rownam
  r$ind <- seq_len(nrow(r))
  r <- rbind(r,c(NA,NA,NA,p_gen,ind=14.1),c(NA,NA,NA,p_risk,ind=15.1),c(NA,NA,NA,p_site,ind=18.1),c(NA,NA,NA,p_mig,ind=22.1))
  r <- r %>% arrange(ind) %>% dplyr::select(-ind)
  return(r)
}


impu.la <- impu.la.post(imp_la_post,vt)
cd4s <- seq(min(la_t_post$sq_cd4,na.rm = T),max(la_t_post$sq_cd4,na.rm=T), by = 0.01)
sp_cd4 <- as.data.frame(ns(cd4s, df=4))
sp_cd4 <- cbind(cd4s, sp_cd4) 
sp_cd4 <- as.data.frame(sp_cd4)
rownames(sp_cd4) <- sp_cd4$cd4
sp_cd4 <- sp_cd4[, -1]
colnames(sp_cd4) <- c("X1", "X2", "X3", "X4")

ages <- round(min(la_t_post$age_b)):round(max(la_t_post$age_b))
sp <- as.data.frame(ns(ages, df=4))
sp <- cbind(ages, sp) 
sp <- as.data.frame(sp)
rownames(sp) <- sp$age
sp <- sp[, -1]
colnames(sp) <- c("X1", "X2", "X3", "X4")

yr <- round(min(la_t_post$baseline_y,na.rm = TRUE)):round(max(la_t_post$baseline_y,na.rm = TRUE))
sp_yr <- as.data.frame(ns(yr, df=4))
sp_yr <- cbind(yr, sp_yr) 
sp_yr <- as.data.frame(sp_yr)
rownames(sp_yr) <- sp_yr$yr
sp_yr <- sp_yr[, -1]
colnames(sp_yr) <- c("X1", "X2", "X3", "X4")

rna <- round(min(la_t_post$log_rna,na.rm = TRUE)):round(max(la_t_post$log_rna,na.rm = TRUE))
sp_rna <- as.data.frame(ns(rna, df=4))
sp_rna <- cbind(rna,sp_rna) 
sp_rna <- as.data.frame(sp_rna)
rownames(sp_rna) <- sp_rna$rna
sp_rna <- sp_rna[, -1]
colnames(sp_rna) <- c("X1", "X2", "X3", "X4")

summ_spline  <- summary(pool(imp_la_post))
v_spline <- vt
se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
summ_spline <- summ_spline[-1,c(1,2)]
summ_spline <- as.data.frame(cbind(summ_spline, se))
summ_spline$vars <- rownames(summ_spline)
colnames(summ_spline) <- c("Variables","Coefs", "SE" )

rr_age_20 <- rr(20,summ_spline)
rr_age_25 <- rr(25,summ_spline)
rr_age_30 <- rr(30,summ_spline)
rr_age_40 <- rr(40,summ_spline)
rr_age_45 <- rr(45,summ_spline)
rr_age_50 <- rr(50,summ_spline)
rr_age_55 <- rr(55,summ_spline)
rr_age_60 <- rr(60,summ_spline)

SE_20 <- spline_se_fun(20)
SE_25 <- spline_se_fun(25)
SE_30 <- spline_se_fun(30)
SE_40 <- spline_se_fun(40)
SE_45 <- spline_se_fun(45)
SE_50 <- spline_se_fun(50)
SE_55 <- spline_se_fun(55)
SE_60 <- spline_se_fun(60)

CI_20 <- CI(rr_age_20,SE_20)
CI_25 <- CI(rr_age_25,SE_25)
CI_30 <- CI(rr_age_30,SE_30)
CI_40 <- CI(rr_age_40,SE_40)
CI_45 <- CI(rr_age_45,SE_45)
CI_50 <- CI(rr_age_50,SE_50)
CI_55 <- CI(rr_age_55,SE_55)
CI_60 <- CI(rr_age_60,SE_60)

abs_20 <- abs_val(20,summ_spline)
abs_25 <- abs_val(25,summ_spline)
abs_30 <- abs_val(30,summ_spline)
abs_40 <- abs_val(40,summ_spline)
abs_45 <- abs_val(45,summ_spline)
abs_50 <- abs_val(50,summ_spline)
abs_55 <- abs_val(55,summ_spline)
abs_60 <- abs_val(60,summ_spline)

rr_sq_5 <- rr_sq_sh(5,summ_spline)
rr_sq_10 <- rr_sq_sh(10,summ_spline)
rr_sq_18.70 <- rr_sq_sh(18.70,summ_spline)
rr_sq_22.40 <- rr_sq_sh(22.40,summ_spline)


SE_sq5 <- spline_se_sh(5)
SE_sq10 <- spline_se_sh(10)
SE_sq18.70 <- spline_se_sh(18.70)
SE_sq22.40 <- spline_se_sh(22.40)

CI_sq5 <- CI(rr_sq_5,SE_5)
CI_sq10 <- CI(rr_sq_10,SE_10)
CI_sq18.70 <- CI(rr_sq_18.70,SE_18.70)
CI_sq22.40 <- CI(rr_sq_22.40,SE_22.40)

abs_sq5 <- abs_sh(5,summ_spline)
abs_sq10 <- abs_sh(10,summ_spline)
abs_sq18.70 <- abs_sh(18.70,summ_spline)
abs_sq22.40 <- abs_sh(22.40,summ_spline)


rr_yr_2005 <- rr_yr(2005,summ_spline)
rr_yr_2000 <- rr_yr(2000,summ_spline)
rr_yr_2015 <- rr_yr(2015,summ_spline)
rr_yr_2020 <- rr_yr(2020,summ_spline)

SE_2005 <- spline_se_yr(2005)
SE_2000 <- spline_se_yr(2000)
SE_2015 <- spline_se_yr(2015)
SE_2020 <- spline_se_yr(2020)


CI_2005 <- CI(rr_yr_2005,SE_2005)
CI_2000 <- CI(rr_yr_2000,SE_2000)
CI_2015 <- CI(rr_yr_2015,SE_2015)
CI_2020 <- CI(rr_yr_2020,SE_2020)



abs_2005 <- abs_val_yrs(2005,summ_spline)
abs_2000 <- abs_val_yrs(2000,summ_spline)
abs_2015 <- abs_val_yrs(2015,summ_spline)
abs_2020 <- abs_val_yrs(2020,summ_spline)

rr_rna_1 <- rr_rna(1,summ_spline)
rr_rna_3 <- rr_rna(3,summ_spline)
rr_rna_7 <- rr_rna(7,summ_spline)


SE_1 <- spline_se_rna(1)
SE_3 <- spline_se_rna(3)
SE_7 <- spline_se_rna(7)



CI_1 <- CI(rr_rna_1,SE_1)
CI_3 <- CI(rr_rna_3,SE_3)
CI_7 <- CI(rr_rna_7,SE_7)




abs_1 <- abs_val_rna(1,summ_spline)
abs_3 <- abs_val_rna(3,summ_spline)
abs_7 <- abs_val_rna(7,summ_spline)

age_spl_multi <- data.frame("RR" = c(rr_age_20,rr_age_25,rr_age_30,rr_age_40,rr_age_45,rr_age_50,rr_age_55,rr_age_60),
                            "Lower CI"= c(CI_20[1],CI_25[1],CI_30[1],CI_40[1],CI_45[1],CI_50[1],CI_55[1],CI_60[1]),
                            "Upper CI"= c(CI_20[2],CI_25[2],CI_30[2],CI_40[2],CI_45[2],CI_50[2],CI_55[2],CI_60[2]),
                            "p-value" = c(2*(1-pnorm(abs(abs_20/SE_20))), 2*(1-pnorm(abs(abs_25/SE_25))), 2*(1-pnorm(abs(abs_30/SE_30))),
                                          2*(1-pnorm(abs(abs_40/SE_40))), 2*(1-pnorm(abs(abs_45/SE_45))), 2*(1-pnorm(abs(abs_50/SE_50))),
                                          2*(1-pnorm(abs(abs_55/SE_55))), 2*(1-pnorm(abs(abs_60/SE_60)))))
colnames(age_spl_multi) <- c("RR","Lower CI","Upper CI","p-value")


cd4_spline_multi <- data.frame("RR" = c(rr_sq_5,rr_sq_10,rr_sq_18.70,rr_sq_22.40),
                               "Lower CI"= c(CI_sq5[1],CI_sq10[1],CI_sq18.70[1],CI_sq22.40[1]),
                               "Upper CI"= c(CI_sq5[2],CI_sq10[2],CI_sq18.70[2],CI_sq22.40[2]),
                               "p-value" = c(2*(1-pnorm(abs(abs_sq5/SE_sq5))), 2*(1-pnorm(abs(abs_sq10/SE_sq10))), 2*(1-pnorm(abs(abs_sq18.70/SE_sq18.70))),
                                             2*(1-pnorm(abs(abs_sq22.40/SE_sq22.40)))))
colnames(cd4_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")


yr_spline_multi <- data.frame("RR" = c(rr_yr_2005,rr_yr_2000,rr_yr_2015,rr_yr_2020),
                              "Lower CI"= c(CI_2005[1],CI_2000[1],CI_2015[1],CI_2020[1]),
                              "Upper CI"= c(CI_2005[2],CI_2000[2],CI_2015[2],CI_2020[2]),
                              "p-value" = c(2*(1-pnorm(abs(abs_2005/SE_2005))), 2*(1-pnorm(abs(abs_2000/SE_2000))),
                                            2*(1-pnorm(abs(abs_2015/SE_2015))), 2*(1-pnorm(abs(abs_2020/SE_2020)))))
colnames(yr_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")

rna_spline_multi <- data.frame("RR" = c(rr_rna_1,rr_rna_3,rr_rna_7),
                               "Lower CI"= c(CI_1[1],CI_3[1],CI_7[1]),
                               "Upper CI"= c(CI_1[2],CI_3[2],CI_7[2]),
                               "p-value" = c(2*(1-pnorm(abs(abs_1/SE_1))), 2*(1-pnorm(abs(abs_3/SE_3))),
                                             2*(1-pnorm(abs(abs_7/SE_7)))))
colnames(rna_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")

impu.la <- impu.la[-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14),]
tab13 <- rbind(c(NA,NA,NA,p_age),age_spl_multi,c(NA,NA,NA,p_cd),cd4_spline_multi,c(NA,NA,NA,p_rna),rna_spline_multi,c(NA,NA,NA,p_yr),yr_spline_multi,impu.la)
rownames(tab13) <- c("Age (ref=35)","20","25","30","40","45","50","55","60",
                     "CD4 (ref=200)","25 ","100","350","500",
                     "log HIV RNA(ref=5)","1","3 "  ,"7","Year (ref=2010)", "2000","2005","2015","2020","Sex",
                     "Male ", "Risk factor (ref =MSM)","Hetero ","IDU ", "Other/Unknown ", "Country (ref=Brazil)", "Chile ", "Honduras ","Mexico ", "Peru "
                     ,"Migration status (ref=No)","Yes")

save(tab13, file="tab13.Rdata")

## imputation model Tn ##
imp_tn_post <- tn_t_post %>% select(CFAR_PID,age_b,sq_cd4,log_rna, baseline_y,risk1,art_status,mig_ende,follow_up,histo_num,birthsex)
imp_tn_post <- mice(imp_tn_post,m=5,printFlag = FALSE, idvars = "CFAR_PID")
set.seed(7)
#tn_imp <- aregImpute(~ age_b + sq_cd4+ log_rna + baseline_y+ as.factor(birthsex)  + as.factor(risk1)   + as.factor(mig_ende) + follow_up, data = tn_t_post, n.impute=20,x=TRUE, nk=0)
imp.tn.post <- with(imp_tn_post,
                    glm(histo_num ~ age_b + sq_cd4+ log_rna + baseline_y+ as.factor(birthsex)  + offset(log(follow_up)) + as.factor(risk1)   + as.factor(mig_ende), 
                        family = poisson(link = "log")))


imp.tn.age <-  with(imp_tn_post,
                    glm(histo_num ~sq_cd4+ log_rna + baseline_y+ as.factor(birthsex)  + offset(log(follow_up)) + as.factor(risk1)   + as.factor(mig_ende), 
                        family = poisson(link = "log")))

a_tng <- D1(imp.tn.post,imp.tn.age)
p_tng <- round(a_tng$result[4],digits=3)

imp.tn.cd <- with(imp_tn_post,
                  glm(histo_num ~ age_b + log_rna + baseline_y+ as.factor(birthsex)  + offset(log(follow_up)) + as.factor(risk1)   + as.factor(mig_ende), 
                      family = poisson(link = "log")))

a_tncd <- D1(imp.tn.post,imp.tn.cd)
p_tncd <- round(a_tncd$result[4],digits=3)

imp.tn.rna <- with(imp_tn_post,
                   glm(histo_num ~ age_b + sq_cd4 + baseline_y+ as.factor(birthsex)  + offset(log(follow_up)) + as.factor(risk1)   + as.factor(mig_ende), 
                       family = poisson(link = "log")))
a_tnrna <- D1(imp.tn.post,imp.tn.rna)
p_tnrna <- round(a_tnrna$result[4],digits=3)

imp.tn.yr <- with(imp_tn_post,
                  glm(histo_num ~ age_b + sq_cd4+ log_rna + as.factor(birthsex)  + offset(log(follow_up)) + as.factor(risk1)   + as.factor(mig_ende),
                      family = poisson(link = "log")))
a_tnyr <- D1(imp.tn.post,imp.tn.yr)
p_tnyr <- round(a_tnyr$result[4],digits=3)


imp.tn.gen <- with(imp_tn_post,
                   glm(histo_num ~  age_b + sq_cd4+ log_rna + baseline_y  + offset(log(follow_up)) + as.factor(risk1)   + as.factor(mig_ende), 
                       family = poisson(link = "log")))

a_tngen <- D1(imp.tn.post,imp.tn.gen)
p_tngen <- round(a_tngen$result[4],digits=3)

imp.tn.risk <- with(imp_tn_post,
                    glm(histo_num ~  age_b + sq_cd4+ log_rna + baseline_y+ as.factor(birthsex)  + offset(log(follow_up))    + as.factor(mig_ende), 
                        family = poisson(link = "log")))

a_tnrisk <- D1(imp.tn.post,imp.tn.risk)
p_tnrisk <- round(a_tnrisk$result[4],digits=3)


imp.tn.mig <- with(imp_tn_post,
                   glm(histo_num ~ age_b + sq_cd4+ log_rna + baseline_y+ as.factor(birthsex)  + offset(log(follow_up)) + as.factor(risk1)  ,
                       family = poisson(link = "log")))


a_tnmig <- D1(imp.tn.post,imp.tn.mig)
p_tnmig <- round(a_tnmig$result[4],digits=3)

## Variance covariance
t2 <- pool(imp.tn.post)
m <- t2$m

## Applying rubin's rule
vcov_list <- list()
for (i in 1:m) {
  vcov_list[[i]] <- vcovHC(imp.tn.post$analyses[[i]])
} 

vw <- Reduce('+',vcov_list)/m
#vcov_array <- array(unlist(vcov_list), dim = c(dim(vcov_list[[1]]),length(vcov_list)))
## Computing within imputaion variance 
#var_vcov <- apply(vcov_array,c(1,2),var)
#vw <- var_vcov
## Between imputaion variance covariance
qbar <- getqbar(t2)
qhats <- sapply(imp.tn.post$analyses, coef)
vb <- (1 / (m-1)) * (qhats - qbar) %*% t(qhats - qbar)
vt <- vw + (1 + 1 / (m)) * vb 
vt <- vt[-c(1), -c(1)]

impu.tn_post <- function(model,var) {
  s <-  summary(pool(model))
  s <-  s$estimate[-1]
  sand <- sqrt(diag(var))
  r = matrix(0,nrow = length(s), ncol = 4)
  r[1,] <- c(exp(s[1]*10), exp((s[1] - 1.96*sand[2])*10),exp((s[1] + 1.96*sand[2])*10),round(2*(1-pnorm(abs(s[1]/sand[2]))),digits =5))
  r[2,] <- c(exp(s[2]*50), exp((s[2] - 1.96*sand[3])*50),exp((s[2] + 1.96*sand[3])*50),round(2*(1-pnorm(abs(s[2]/sand[3]))),digits =5))
  r[3,] <- c(exp(s[2]*4), exp((s[2] - 1.96*sand[3])*4),exp((s[2] + 1.96*sand[3])*4),round(2*(1-pnorm(abs(s[2]/sand[3]))),digits =5))
  
  for (i in 4:(length(s))){
    r[i,] <- c(exp(s[i]), exp(s[i] - 1.96*sand[i]),exp(s[i] + 1.96*sand[i]),round((1 - pnorm(abs(s[i]/sand[i]))) * 2,digits=5))
  }
  
  r <- as.data.frame(r)
  colnames(r) <- c( "RR","Lower CI", "Upper CI","p-value") 
  rownam <- rownames(s)
  rownames(r) <- rownam
  r$ind <- seq_len(nrow(r))
  r <- rbind(r,c(NA,NA,NA,p_tngen,ind=4.1),c(NA,NA,NA,p_tnrisk,ind=5.1),c(NA,NA,NA,p_tnmig,ind=8.1))
  r <- r %>% arrange(ind) %>% dplyr::select(-ind)
  return(r)
}


# cd4s <- seq(min(tn_t_post$sq_cd4,na.rm = T),max(tn_t_post$sq_cd4,na.rm=T), by = 0.01)
# sp_cd4 <- as.data.frame(ns(cd4s, df=4))
# sp_cd4 <- cbind(cd4s, sp_cd4) 
# sp_cd4 <- as.data.frame(sp_cd4)
# rownames(sp_cd4) <- sp_cd4$cd4
# sp_cd4 <- sp_cd4[, -1]
# colnames(sp_cd4) <- c("X1", "X2", "X3", "X4")
# 
# ages <- round(min(tn_t_post$age_b)):round(max(tn_t_post$age_b))
# sp <- as.data.frame(ns(ages, df=4))
# sp <- cbind(ages, sp) 
# sp <- as.data.frame(sp)
# rownames(sp) <- sp$age
# sp <- sp[, -1]
# colnames(sp) <- c("X1", "X2", "X3", "X4")
# 
# yr <- round(min(tn_t_post$baseline_y,na.rm = TRUE)):round(max(tn_t_post$baseline_y,na.rm = TRUE))
# sp_yr <- as.data.frame(ns(yr, df=4))
# sp_yr <- cbind(yr, sp_yr) 
# sp_yr <- as.data.frame(sp_yr)
# rownames(sp_yr) <- sp_yr$yr
# sp_yr <- sp_yr[, -1]
# colnames(sp_yr) <- c("X1", "X2", "X3", "X4")
# 
# rna <- round(min(tn_t_post$log_rna,na.rm = TRUE)):round(max(tn_t_post$log_rna,na.rm = TRUE))
# sp_rna <- as.data.frame(ns(rna, df=4))
# sp_rna <- cbind(rna,sp_rna) 
# sp_rna <- as.data.frame(sp_rna)
# rownames(sp_rna) <- sp_rna$rna
# sp_rna <- sp_rna[, -1]
# colnames(sp_rna) <- c("X1", "X2", "X3", "X4")
# 
# summ_spline  <- summary(pool(imp.tn.post))
# v_spline <- vt
# se <- sqrt(v_spline[row(v_spline)==col(v_spline)])
# summ_spline <- summ_spline[-1,c(1,2)]
# summ_spline <- as.data.frame(cbind(summ_spline, se))
# summ_spline$vars <- rownames(summ_spline)
# colnames(summ_spline) <- c("Variables","Coefs", "SE" )
# 
# rr_age_20 <- rr(20,summ_spline)
# rr_age_25 <- rr(25,summ_spline)
# rr_age_30 <- rr(30,summ_spline)
# rr_age_40 <- rr(40,summ_spline)
# rr_age_45 <- rr(45,summ_spline)
# rr_age_50 <- rr(50,summ_spline)
# rr_age_55 <- rr(55,summ_spline)
# rr_age_60 <- rr(60,summ_spline)
# 
# SE_20 <- spline_se_fun(20)
# SE_25 <- spline_se_fun(25)
# SE_30 <- spline_se_fun(30)
# SE_40 <- spline_se_fun(40)
# SE_45 <- spline_se_fun(45)
# SE_50 <- spline_se_fun(50)
# SE_55 <- spline_se_fun(55)
# SE_60 <- spline_se_fun(60)
# 
# CI_20 <- CI(rr_age_20,SE_20)
# CI_25 <- CI(rr_age_25,SE_25)
# CI_30 <- CI(rr_age_30,SE_30)
# CI_40 <- CI(rr_age_40,SE_40)
# CI_45 <- CI(rr_age_45,SE_45)
# CI_50 <- CI(rr_age_50,SE_50)
# CI_55 <- CI(rr_age_55,SE_55)
# CI_60 <- CI(rr_age_60,SE_60)
# 
# abs_20 <- abs_val(20,summ_spline)
# abs_25 <- abs_val(25,summ_spline)
# abs_30 <- abs_val(30,summ_spline)
# abs_40 <- abs_val(40,summ_spline)
# abs_45 <- abs_val(45,summ_spline)
# abs_50 <- abs_val(50,summ_spline)
# abs_55 <- abs_val(55,summ_spline)
# abs_60 <- abs_val(60,summ_spline)
# 
# rr_sq_5 <- rr_sq_sh(5,summ_spline)
# rr_sq_10 <- rr_sq_sh(10,summ_spline)
# rr_sq_18.70 <- rr_sq_sh(18.70,summ_spline)
# rr_sq_22.40 <- rr_sq_sh(22.40,summ_spline)
# 
# 
# SE_sq5 <- spline_se_sh(5)
# SE_sq10 <- spline_se_sh(10)
# SE_sq18.70 <- spline_se_sh(18.70)
# SE_sq22.40 <- spline_se_sh(22.40)
# 
# CI_sq5 <- CI(rr_sq_5,SE_5)
# CI_sq10 <- CI(rr_sq_10,SE_10)
# CI_sq18.70 <- CI(rr_sq_18.70,SE_18.70)
# CI_sq22.40 <- CI(rr_sq_22.40,SE_22.40)
# 
# abs_sq5 <- abs_sh(5,summ_spline)
# abs_sq10 <- abs_sh(10,summ_spline)
# abs_sq18.70 <- abs_sh(18.70,summ_spline)
# abs_sq22.40 <- abs_sh(22.40,summ_spline)
# 
# 
# rr_yr_2005 <- rr_yr(2005,summ_spline)
# rr_yr_2000 <- rr_yr(2000,summ_spline)
# rr_yr_2015 <- rr_yr(2015,summ_spline)
# rr_yr_2020 <- rr_yr(2020,summ_spline)
# 
# SE_2005 <- spline_se_yr(2005)
# SE_2000 <- spline_se_yr(2000)
# SE_2015 <- spline_se_yr(2015)
# SE_2020 <- spline_se_yr(2020)
# 
# 
# CI_2005 <- CI(rr_yr_2005,SE_2005)
# CI_2000 <- CI(rr_yr_2000,SE_2000)
# CI_2015 <- CI(rr_yr_2015,SE_2015)
# CI_2020 <- CI(rr_yr_2020,SE_2020)
# 
# 
# 
# abs_2005 <- abs_val_yrs(2005,summ_spline)
# abs_2000 <- abs_val_yrs(2000,summ_spline)
# abs_2015 <- abs_val_yrs(2015,summ_spline)
# abs_2020 <- abs_val_yrs(2020,summ_spline)
# 
# rr_rna_1 <- rr_rna(1,summ_spline)
# rr_rna_3 <- rr_rna(3,summ_spline)
# rr_rna_7 <- rr_rna(7,summ_spline)
# 
# 
# SE_1 <- spline_se_rna(1)
# SE_3 <- spline_se_rna(3)
# SE_7 <- spline_se_rna(7)
# 
# 
# 
# CI_1 <- CI(rr_rna_1,SE_1)
# CI_3 <- CI(rr_rna_3,SE_3)
# CI_7 <- CI(rr_rna_7,SE_7)
# 
# 
# 
# 
# abs_1 <- abs_val_rna(1,summ_spline)
# abs_3 <- abs_val_rna(3,summ_spline)
# abs_7 <- abs_val_rna(7,summ_spline)
# 
# age_spl_multi <- data.frame("RR" = c(rr_age_20,rr_age_25,rr_age_30,rr_age_40,rr_age_45,rr_age_50,rr_age_55,rr_age_60),
#                             "Lower CI"= c(CI_20[1],CI_25[1],CI_30[1],CI_40[1],CI_45[1],CI_50[1],CI_55[1],CI_60[1]),
#                             "Upper CI"= c(CI_20[2],CI_25[2],CI_30[2],CI_40[2],CI_45[2],CI_50[2],CI_55[2],CI_60[2]),
#                             "p-value" = c(2*(1-pnorm(abs(abs_20/SE_20))), 2*(1-pnorm(abs(abs_25/SE_25))), 2*(1-pnorm(abs(abs_30/SE_30))),
#                                           2*(1-pnorm(abs(abs_40/SE_40))), 2*(1-pnorm(abs(abs_45/SE_45))), 2*(1-pnorm(abs(abs_50/SE_50))),
#                                           2*(1-pnorm(abs(abs_55/SE_55))), 2*(1-pnorm(abs(abs_60/SE_60)))))
# colnames(age_spl_multi) <- c("RR","Lower CI","Upper CI","p-value")
# 
# 
# cd4_spline_multi <- data.frame("RR" = c(rr_sq_5,rr_sq_10,rr_sq_18.70,rr_sq_22.40),
#                                "Lower CI"= c(CI_sq5[1],CI_sq10[1],CI_sq18.70[1],CI_sq22.40[1]),
#                                "Upper CI"= c(CI_sq5[2],CI_sq10[2],CI_sq18.70[2],CI_sq22.40[2]),
#                                "p-value" = c(2*(1-pnorm(abs(abs_sq5/SE_sq5))), 2*(1-pnorm(abs(abs_sq10/SE_sq10))), 2*(1-pnorm(abs(abs_sq18.70/SE_sq18.70))),
#                                              2*(1-pnorm(abs(abs_sq22.40/SE_sq22.40)))))
# colnames(cd4_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")
# 
# 
# yr_spline_multi <- data.frame("RR" = c(rr_yr_2005,rr_yr_2000,rr_yr_2015,rr_yr_2020),
#                               "Lower CI"= c(CI_2005[1],CI_2000[1],CI_2015[1],CI_2020[1]),
#                               "Upper CI"= c(CI_2005[2],CI_2000[2],CI_2015[2],CI_2020[2]),
#                               "p-value" = c(2*(1-pnorm(abs(abs_2005/SE_2005))), 2*(1-pnorm(abs(abs_2000/SE_2000))),
#                                             2*(1-pnorm(abs(abs_2015/SE_2015))), 2*(1-pnorm(abs(abs_2020/SE_2020)))))
# colnames(yr_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")
# 
# rna_spline_multi <- data.frame("RR" = c(rr_rna_1,rr_rna_3,rr_rna_7),
#                                "Lower CI"= c(CI_1[1],CI_3[1],CI_7[1]),
#                                "Upper CI"= c(CI_1[2],CI_3[2],CI_7[2]),
#                                "p-value" = c(2*(1-pnorm(abs(abs_1/SE_1))), 2*(1-pnorm(abs(abs_3/SE_3))),
#                                              2*(1-pnorm(abs(abs_7/SE_7)))))
# colnames(rna_spline_multi) <- c("RR","Lower CI","Upper CI","p-value")
tab14 <- impu.tn_post(imp.tn.post,vt)

#impu.tn <- impu.tn[-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14),]
#tab14 <- rbind(c(NA,NA,NA,p_tng),age_spl_multi,c(NA,NA,NA,p_tncd),cd4_spline_multi,c(NA,NA,NA,p_tnrna),rna_spline_multi,c(NA,NA,NA,p_tnyr),yr_spline_multi,impu.tn)

rownames(tab14) <- c("Age 35(ref=25)",
                     "CD4 100 (ref=50)",
                     "log HIV RNA1 (ref=5)","Year 2005 (ref=2015) ","Sex",
                     "Male ", "Risk factor (ref =MSM)","Hetero ","IDU ", "Other/Unknown ", 
                     "Migration status (ref=No)","Yes")

save(tab14, file="tab14.Rdata")

### ------------------ ###
###  Survival Analysis ###
### ------------------ ###

tn$site <- "Tennessee"
tn$patient_id <- tn$CFAR_PID
tn$histo_y <- tn$HISTO_DX_YEAR 
tn$age_at_diag <- tn$age_at_histo
la$age_at_diag <- with(la,as.numeric(abs((difftime(as.Date(birth_d),as.Date(histo_d),units ="days")/365.25))))
tn$death_d <- with(tn,ifelse(!is.na(age_at_death), enrol_d + (age_at_death - age_at_first_visit)*365.25,NA))
tn_c <- tn %>% filter(!is.na(histo_d)) %>% dplyr::select(patient_id,baseline,baseline_y,histo_d,sex,site,histo_y,art_sd,ade_type,rna_diag_log,cd4_diag_sq,last_d,age_at_diag,death_d,mig_ende)
la_c <- la %>% filter(!is.na(histo_d)) %>% dplyr::select(patient_id,baseline,baseline_y,histo_d,sex,site,histo_y,art_sd,ade_type,rna_diag_log,cd4_diag_sq,last_d,age_at_diag,death_d,mig_ende)

#tn_c <- tn_c %>% group_by(patient_id) %>% dplyr::select(-CFAR_PID)
#c <- rbind(la_c, tn_c)

la_c$time <- with(la_c,ifelse(!is.na(death_d),difftime(as.Date(death_d),as.Date(histo_d),units ="days")/365.25,
                              difftime(as.Date(last_d),as.Date(histo_d),units ="days")/365.25))

la_c$death_y <- with(la_c,ifelse(!is.na(death_d), 1,0))
la_c$death <- as.factor(with(la_c,ifelse(death_y ==1,"Yes","No")))
la_c$art_status <- with(la_c,ifelse(!is.na(art_sd) & art_sd < histo_d, "post ART at histo", "pre ART at histo"))
dd <- datadist(la_c)
options(datadist= 'dd')


dd$limits$age_at_diag[2] <- 35
dd$limits$sex[2] <- "Female"
dd$limits$histo_y[2] <- 2010
dd$limits$ade_type[2] <- "No"
dd$limits$cd4_diag_sq[2] <- 50
dd$limits$rna_diag_log[2] <- 5 
dd$limits$art_status[2] <- "pre ART at histo"
dd$limits$mig_ende[2] <- "No"

surv.obj <- with(la_c, Surv(time, death_y==1))

my_cph_uni <- sapply(c("sex",  "art_status", "mig_ende"), function(x)
  as.formula(paste("surv.obj~ ", x, "+ strat(site)")))

uni_mods <- lapply(my_cph_uni, function(x) cph(x, data=la_c, x=TRUE, y=TRUE, surv=TRUE))

make.uni.table.catvar <- function(catvar_mod){
  library(Hmisc)
  library(rms)
  library(plyr)
  library(dplyr)
  clearvars <- c('Effect', 'Lower.0.95', 'Upper.0.95', 'CI', 'S.E.')
  
  uni_model <- as.data.frame(summary(catvar_mod, vnames='labels'), row.names = make.unique(rownames(summary(catvar_mod))))
  uni_model_aov <- as.data.frame(anova(catvar_mod))
  
  uni_model_aov$term <- rownames(uni_model_aov)
  uni_model$term <- rownames(uni_model)
  uni_model$term[1] <- rownames(uni_model_aov)[1]
  
  uni_model$Variable <- rownames(uni_model)
  uni_model_table <- join(uni_model, uni_model_aov)
  uni_model_table$Variable <- ifelse(uni_model_table$Type==2,  gsub("(.+?)(\\:.*)", "\\1", lag(uni_model_table$Variable)), sub("^[^:]*:", "", uni_model_table$Variable))
  uni_model_table$S.E. <- ifelse(uni_model_table$Type==2, lag(uni_model_table$S.E.), uni_model_table$S.E.)
  
  names(uni_model_table) <- gsub(" ", ".", names(uni_model_table))
  names(uni_model_table) <- gsub("-", ".", names(uni_model_table))
  uni_model_table$P <- ifelse(uni_model_table$P <.0001, "<.0001", format(round(uni_model_table$P, digits=4), nsmall=3))
  uni_model_table <- rapply(object = uni_model_table, f = round, classes = "numeric", how = "replace", digits=3)
  
  uni_model_table$CI <- with(uni_model_table, paste0("(", formatC(Lower.0.95, digits=3, format= "f"), ", ", formatC(Upper.0.95, digits=3, format= "f"), ")"), sep="")
  
  uni_model_table[which(uni_model_table$Type==1), clearvars] <- ""
  
  uni_model_table <- subset(uni_model_table, select= c("Variable", "Effect", "CI",  "P"))
  
  names(uni_model_table) <- c("Variable", "HR", "95% CI", "p")
  uni_model_table <- uni_model_table[-c(seq(3, 10, by=2)), ]
  uni_model_table
}

make.uni.table.age <- function(contvar_mod, ref, compare1, compare2){
  library(Hmisc)
  library(rms)
  library(plyr)
  library(dplyr)
  clearvars <- c('Effect', 'Lower.0.95', 'Upper.0.95', 'CI', 'S.E.')
  
  
  uni_cont2 <- as.data.frame(summary(contvar_mod, age_at_diag= c(ref, compare2), vnames="labels"), row.names = make.unique(rownames(summary(contvar_mod))))
  uni_cont1 <- as.data.frame(summary(contvar_mod, age_at_diag= c(ref, compare1), vnames="labels"), row.names = make.unique(rownames(summary(contvar_mod))))
  uni_cont <- rbind(uni_cont1, uni_cont2)
  uni_cont_aov <- as.data.frame(anova(contvar_mod))
  uni_cont_aov$term <- rownames(uni_cont_aov)
  
  rownames(uni_cont) <- ifelse(uni_cont$Type==2, uni_cont$High, rownames(uni_cont))
  uni_cont$term <- rownames(uni_cont)
  uni_cont$S.E. <- ifelse(uni_cont$Type==2, lag(uni_cont$S.E.), uni_cont$S.E.)
  uni_cont_table <- join(uni_cont, uni_cont_aov)
  names(uni_cont_table) <- gsub(" ", ".", names(uni_cont_table))
  names(uni_cont_table) <- gsub("-", ".", names(uni_cont_table))
  uni_cont_table$P <- ifelse(uni_cont_table$P <.0001, "<.0001", format(round(uni_cont_table$P, digits=4), nsmall=3))
  uni_cont_table <- rapply(object = uni_cont_table, f = round, classes = "numeric", how = "replace", digits=3)
  
  uni_cont_table$CI <- with(uni_cont_table, paste0("(", formatC(Lower.0.95, digits=3, format= "f"), ", ", formatC(Upper.0.95, digits=3, format= "f"), ")"), sep="")
  
  uni_cont_table[which(uni_cont_table$Type==1), clearvars] <- ""
  uni_cont_table <- subset(uni_cont_table, select= c("term", "Effect", "CI",  "P"))
  
  names(uni_cont_table) <-  c("Variable", "HR", "95% CI",  "p")
  uni_cont_table <- uni_cont_table[-3,]
  uni_cont_table$Variable <- c(paste("Age (ref: ", ref, ")", sep=''),  compare1, compare2)
  
  uni_cont_table
}

make.uni.table.yrs <- function(contvar_mod, ref, compare1, compare2){
  library(Hmisc)
  library(rms)
  library(plyr)
  library(dplyr)
  clearvars <- c('Effect', 'Lower.0.95', 'Upper.0.95', 'CI', 'S.E.')
  
  uni_cont1 <- as.data.frame(summary(contvar_mod, histo_y=c(ref, compare1), vnames="labels"), row.names = make.unique(rownames(summary(contvar_mod))))
  uni_cont2 <- as.data.frame(summary(contvar_mod, histo_y=c(ref, compare2), vnames="labels"), row.names = make.unique(rownames(summary(contvar_mod))))
  uni_cont <- rbind(uni_cont1, uni_cont2)
  uni_cont_aov <- as.data.frame(anova(contvar_mod))
  uni_cont_aov$term <- rownames(uni_cont_aov)
  
  rownames(uni_cont) <- ifelse(uni_cont$Type==2, uni_cont$High, rownames(uni_cont))
  uni_cont$term <- rownames(uni_cont)
  uni_cont$S.E. <- ifelse(uni_cont$Type==2, lag(uni_cont$S.E.), uni_cont$S.E.)
  uni_cont_table <- join(uni_cont, uni_cont_aov)
  names(uni_cont_table) <- gsub(" ", ".", names(uni_cont_table))
  names(uni_cont_table) <- gsub("-", ".", names(uni_cont_table))
  uni_cont_table$P <- ifelse(uni_cont_table$P <.0001, "<.0001", format(round(uni_cont_table$P, digits=4), nsmall=3))
  uni_cont_table <- rapply(object = uni_cont_table, f = round, classes = "numeric", how = "replace", digits=3)
  
  uni_cont_table$CI <- with(uni_cont_table, paste0("(", formatC(Lower.0.95, digits=3, format= "f"), ", ", formatC(Upper.0.95, digits=3, format= "f"), ")"), sep="")
  
  uni_cont_table[which(uni_cont_table$Type==1), clearvars] <- ""
  uni_cont_table <- subset(uni_cont_table, select= c("term", "Effect", "CI", "P"))
  # uni_cont_table$S.E. <- ifelse(uni_cont_table$Type==1, NA, uni_cont_table$S.E.)
  
  names(uni_cont_table) <- c("Variable", "HR", "95% CI","p")
  uni_cont_table <- uni_cont_table[-3,]
  uni_cont_table$Variable <- c(paste("Years since Histo diagnosis (ref: ", ref, ")", sep=''),  compare1, compare2)
  
  uni_cont_table
}


make.uni.table.art <- function(contvar_mod, ref, compare1, compare2){
  library(Hmisc)
  library(rms)
  library(plyr)
  library(dplyr)
  clearvars <- c('Effect', 'Lower.0.95', 'Upper.0.95', 'CI', 'S.E.')
  
  uni_cont1 <- as.data.frame(summary(contvar_mod, art_missed=c(ref, compare1), vnames="labels"), row.names = make.unique(rownames(summary(contvar_mod))))
  uni_cont2 <- as.data.frame(summary(contvar_mod, art_missed=c(ref, compare2), vnames="labels"), row.names = make.unique(rownames(summary(contvar_mod))))
  uni_cont <- rbind(uni_cont1, uni_cont2)
  uni_cont_aov <- as.data.frame(anova(contvar_mod))
  uni_cont_aov$term <- rownames(uni_cont_aov)
  
  rownames(uni_cont) <- ifelse(uni_cont$Type==2, uni_cont$High, rownames(uni_cont))
  uni_cont$term <- rownames(uni_cont)
  uni_cont$S.E. <- ifelse(uni_cont$Type==2, lag(uni_cont$S.E.), uni_cont$S.E.)
  uni_cont_table <- join(uni_cont, uni_cont_aov)
  names(uni_cont_table) <- gsub(" ", ".", names(uni_cont_table))
  names(uni_cont_table) <- gsub("-", ".", names(uni_cont_table))
  uni_cont_table$P <- ifelse(uni_cont_table$P <.0001, "<.0001", format(round(uni_cont_table$P, digits=4), nsmall=3))
  uni_cont_table <- rapply(object = uni_cont_table, f = round, classes = "numeric", how = "replace", digits=3)
  
  uni_cont_table$CI <- with(uni_cont_table, paste0("(", formatC(Lower.0.95, digits=3, format= "f"), ", ", formatC(Upper.0.95, digits=3, format= "f"), ")"), sep="")
  
  uni_cont_table[which(uni_cont_table$Type==1), clearvars] <- ""
  uni_cont_table <- subset(uni_cont_table, select= c("term", "Effect", "CI", "P"))
  # uni_cont_table$S.E. <- ifelse(uni_cont_table$Type==1, NA, uni_cont_table$S.E.)
  
  names(uni_cont_table) <- c("Variable", "HR", "95% CI","p")
  uni_cont_table <- uni_cont_table[-3,]
  uni_cont_table$Variable <- c(paste("ART doses missed in past 7 days (ref: ", ref, ")", sep=''),  compare1, compare2)
  
  uni_cont_table
}

make.uni.table.cd4 <- function(contvar_mod, ref, compare1, compare2){
  library(Hmisc)
  library(rms)
  library(plyr)
  library(dplyr)
  clearvars <- c('Effect', 'Lower.0.95', 'Upper.0.95', 'CI', 'S.E.')
  
  uni_cont1 <- as.data.frame(summary(contvar_mod, cd4_diag_sq=c(ref, compare1), vnames="labels"), row.names = make.unique(rownames(summary(contvar_mod))))
  uni_cont2 <- as.data.frame(summary(contvar_mod, cd4_diag_sq= c(ref, compare2), vnames="labels"), row.names = make.unique(rownames(summary(contvar_mod))))
  uni_cont <- rbind(uni_cont1, uni_cont2)
  uni_cont_aov <- as.data.frame(anova(contvar_mod))
  uni_cont_aov$term <- rownames(uni_cont_aov)
  
  rownames(uni_cont) <- ifelse(uni_cont$Type==2, uni_cont$High, rownames(uni_cont))
  uni_cont$term <- rownames(uni_cont)
  uni_cont$S.E. <- ifelse(uni_cont$Type==2, lag(uni_cont$S.E.), uni_cont$S.E.)
  uni_cont_table <- join(uni_cont, uni_cont_aov)
  names(uni_cont_table) <- gsub(" ", ".", names(uni_cont_table))
  names(uni_cont_table) <- gsub("-", ".", names(uni_cont_table))
  uni_cont_table$P <- ifelse(uni_cont_table$P <.0001, "<.0001", format(round(uni_cont_table$P, digits=4), nsmall=3))
  uni_cont_table <- rapply(object = uni_cont_table, f = round, classes = "numeric", how = "replace", digits=3)
  
  uni_cont_table$CI <- with(uni_cont_table, paste0("(", formatC(Lower.0.95, digits=3, format= "f"), ", ", formatC(Upper.0.95, digits=3, format= "f"), ")"), sep="")
  
  uni_cont_table[which(uni_cont_table$Type==1), clearvars] <- ""
  uni_cont_table <- subset(uni_cont_table, select= c("term", "Effect", "CI",  "P"))
  
  names(uni_cont_table) <-  c("Variable", "HR", "95% CI",  "p")
  uni_cont_table <- uni_cont_table[-3,]
  uni_cont_table$Variable <- c(paste("CD4 sqrt (ref: ", ref, ")", sep=''),  compare1, compare2)
  
  uni_cont_table
}

make.uni.table.rna <- function(contvar_mod, ref, compare1, compare2){
  library(Hmisc)
  library(rms)
  library(plyr)
  library(dplyr)
  clearvars <- c('Effect', 'Lower.0.95', 'Upper.0.95', 'CI', 'S.E.')
  
  uni_cont1 <- as.data.frame(summary(contvar_mod, rna_diag_log=c(ref, compare1), vnames="labels"), row.names = make.unique(rownames(summary(contvar_mod))))
  uni_cont2 <- as.data.frame(summary(contvar_mod, rna_diag_log= c(ref, compare2), vnames="labels"), row.names = make.unique(rownames(summary(contvar_mod))))
  uni_cont <- rbind(uni_cont1, uni_cont2)
  uni_cont_aov <- as.data.frame(anova(contvar_mod))
  uni_cont_aov$term <- rownames(uni_cont_aov)
  
  rownames(uni_cont) <- ifelse(uni_cont$Type==2, uni_cont$High, rownames(uni_cont))
  uni_cont$term <- rownames(uni_cont)
  uni_cont$S.E. <- ifelse(uni_cont$Type==2, lag(uni_cont$S.E.), uni_cont$S.E.)
  uni_cont_table <- join(uni_cont, uni_cont_aov)
  names(uni_cont_table) <- gsub(" ", ".", names(uni_cont_table))
  names(uni_cont_table) <- gsub("-", ".", names(uni_cont_table))
  uni_cont_table$P <- ifelse(uni_cont_table$P <.0001, "<.0001", format(round(uni_cont_table$P, digits=4), nsmall=3))
  uni_cont_table <- rapply(object = uni_cont_table, f = round, classes = "numeric", how = "replace", digits=3)
  
  uni_cont_table$CI <- with(uni_cont_table, paste0("(", formatC(Lower.0.95, digits=3, format= "f"), ", ", formatC(Upper.0.95, digits=3, format= "f"), ")"), sep="")
  
  uni_cont_table[which(uni_cont_table$Type==1), clearvars] <- ""
  uni_cont_table <- subset(uni_cont_table, select= c("term", "Effect", "CI",  "P"))
  
  names(uni_cont_table) <-  c("Variable", "HR", "95% CI",  "p")
  uni_cont_table <- uni_cont_table[-3,]
  uni_cont_table$Variable <- c(paste("RNA at diagnosis (ref: ", ref, ")", sep=''),  compare1, compare2)
  
  uni_cont_table
}

# gender
cph_uni_sex <- make.uni.table.catvar(uni_mods$sex)
cph_uni_sex$Variable <- c("Sex (ref: Female)", "Male")

# ART status
cph_uni_art<- make.uni.table.catvar(uni_mods$art_status)
cph_uni_art$Variable <- c("ART initiation status at histo(ref=pre-ART)",  "post ART at histo")

# Migration status
cph_uni_reg <- make.uni.table.catvar(uni_mods$mig_ende)
cph_uni_reg$Variable <- c("Migratied from Histo endemic country(ref=no)",  "Yes")

# age
uni_age <- cph(surv.obj~ age_at_diag + strat(site), data=la_c)
cph_uni_age <- make.uni.table.age(uni_age, 40, 20, 30)
cph_uni_age2 <- make.uni.table.age(uni_age, 40, 50, 60)
# Histo year
uni_yrs <- cph(surv.obj~ histo_y + strat(site), data=la_c)
cph_uni_yr <- make.uni.table.yrs(uni_yrs, 2010, 2000, 2006)
cph_uni_yr2 <- make.uni.table.yrs(uni_yrs, 2010, 2015, 2018)

# cd4
mymod_cd4 <- cph(surv.obj~ cd4_diag_sq + strat(site), data=la_c)
cph_uni_cd4 <- make.uni.table.cd4(mymod_cd4, 7, 3, 5)
cph_uni_cd42 <- make.uni.table.cd4(mymod_cd4, 7,9,14.14214)

# rna log
mymod_rna <- cph(surv.obj~ rna_diag_log + strat(site), data=la_c)
cph_uni_rna <- make.uni.table.rna(mymod_rna, 4, 2, 3)
cph_uni_rna2 <- make.uni.table.rna(mymod_rna, 4, 5, 6)

cph_uni_tab_la <- rbind(cph_uni_sex, cph_uni_age, cph_uni_age2[-1, ], cph_uni_reg,cph_uni_yr, cph_uni_yr2[-1, ], cph_uni_cd4, cph_uni_cd42[-1, ],cph_uni_rna, cph_uni_rna2[-1, ],cph_uni_art)
rownames(cph_uni_tab_la) <- NULL


## MODELS for TN ##

tn_c$time <- with(tn_c,ifelse(!is.na(death_d),difftime(as.Date(death_d),as.Date(histo_d),units ="days")/365.25,
                              difftime(as.Date(last_d),as.Date(histo_d),units ="days")/365.25))

tn_c$death_y <- with(tn_c,ifelse(!is.na(death_d), 1,0))
tn_c$death <- as.factor(with(tn_c,ifelse(death_y ==1,"Yes","No")))
tn_c$art_status <- with(tn_c,ifelse(!is.na(art_sd) & art_sd < histo_d, "post ART at histo", "pre ART at histo"))
dd <- datadist(tn_c)
options(datadist= 'dd')


dd$limits$age_at_diag[2] <- 35
dd$limits$sex[2] <- "Female"
dd$limits$histo_y[2] <- 2010
dd$limits$cd4_diag_sq[2] <- 50
dd$limits$rna_diag_log[2] <- 5 
dd$limits$art_status[2] <- "pre ART at histo"
dd$limits$mig_ende[2] <- "No"

surv.obj <- with(tn_c, Surv(time, death_y==1))


my_cph_uni <- sapply(c("sex",  "art_status", "mig_ende"), function(x)
  as.formula(paste("surv.obj~ ", x)))

uni_mods <- lapply(my_cph_uni, function(x) cph(x, data=tn_c, x=TRUE, y=TRUE, surv=TRUE))
# gender
cph_uni_sex <- make.uni.table.catvar(uni_mods$sex)
cph_uni_sex$Variable <- c("Sex (ref: Female)", "Male")

# ART status
cph_uni_art<- make.uni.table.catvar(uni_mods$art_status)
cph_uni_art$Variable <- c("ART initiation status at histo(ref=pre-ART)",  "post ART at histo")

# Migration status
cph_uni_reg <- make.uni.table.catvar(uni_mods$mig_ende)
cph_uni_reg$Variable <- c("Migratied from Histo endemic country(ref=no)",  "Yes")

# age
uni_age <- cph(surv.obj~ age_at_diag , data=tn_c)
cph_uni_age <- make.uni.table.age(uni_age, 40, 20, 30)
cph_uni_age2 <- make.uni.table.age(uni_age, 40, 50, 60)
# Histo year
uni_yrs <- cph(surv.obj~ histo_y , data=tn_c)
cph_uni_yr <- make.uni.table.yrs(uni_yrs, 2010, 2000, 2006)
cph_uni_yr2 <- make.uni.table.yrs(uni_yrs, 2010, 2015, 2018)

# cd4
mymod_cd4 <- cph(surv.obj~ cd4_diag_sq , data=tn_c)
cph_uni_cd4 <- make.uni.table.cd4(mymod_cd4, 7.071068, 3.162278, 5)
cph_uni_cd42 <- make.uni.table.cd4(mymod_cd4, 7.071068, 9, 14.14214)

# rna log
mymod_rna <- cph(surv.obj~ rna_diag_log , data=tn_c)
cph_uni_rna <- make.uni.table.rna(mymod_rna, 4, 2, 3)
cph_uni_rna2 <- make.uni.table.rna(mymod_rna, 4, 5, 6)

cph_uni_tab_tn <- rbind(cph_uni_sex, cph_uni_age, cph_uni_age2[-1, ], cph_uni_reg,cph_uni_yr, cph_uni_yr2[-1, ], cph_uni_cd4, cph_uni_cd42[-1, ],cph_uni_rna, cph_uni_rna2[-1, ],cph_uni_art)
rownames(cph_uni_tab_tn) <- NULL


## Overall
c <- rbind(la_c, tn_c)
dd <- datadist(c)
options(datadist= 'dd')


dd$limits$age_at_diag[2] <- 35
dd$limits$sex[2] <- "Female"
dd$limits$histo_y[2] <- 2010
dd$limits$cd4_diag_sq[2] <- 50
dd$limits$rna_diag_log[2] <- 5 
dd$limits$art_status[2] <- "pre ART at histo"
dd$limits$mig_ende[2] <- "No"
Variable <- c("Sex (ref: Female)", "Male", 
              "Migratied from Histo endemic country(ref=no)", "Yes",
              "Age at diagnosis (ref: 35)", "25",  "45", "55", "65",
              "Histo diagnosis year (ref: 2010)", "2000", "2005", "2015","2020",
              "CD4 at histo diagnosis (ref: 50)", "10" ," 25", "100", "200",
              "log HIV RNA at diagnosis (ref: 5)", "1","3","7","6","ART initiation status wrt histo ref:pre ART at histo","post ART at histo")
cph_uni_tab_la <- cph_uni_tab_la[,-1]
cph_uni_tab_tn <- cph_uni_tab_tn[,-1]
rownames(cph_uni_tab_la) <- Variable
rownames(cph_uni_tab_tn) <- Variable
save(cph_uni_tab_la, file="cph_uni_tab_la.Rdata")
save(cph_uni_tab_tn, file="cph_uni_tab_tn.Rdata")



## Multivariate analysis with TN as a site
c <- rbind(la_c, tn_c)

surv.obj <- with(c, Surv(time, death_y==1))
mymod <- cph(surv.obj~ sex + age_at_diag + mig_ende + histo_y + cd4_diag_sq + rna_diag_log + art_status + strat(site), data=c,x=TRUE, y=TRUE, surv=TRUE)

summ1 <- as.data.frame(summary(mymod, age_at_diag=c(35, 25), cd4_diag_sq=c(7.071068, 3.162278), histo_y=c(2010, 2000), rna_diag_log=c(5,1)))
summ2 <- as.data.frame(summary(mymod, age_at_diag=c(35, 45), cd4_diag_sq=c(7.071068, 5), histo_y=c(2010, 2005), rna_diag_log=c(5,3)))
summ3 <- as.data.frame(summary(mymod, age_at_diag=c(35, 55), cd4_diag_sq=c(7.071068, 10),histo_y=c(2010, 2015), rna_diag_log=c(5,7)))
summ4 <- as.data.frame(summary(mymod, age_at_diag=c(35, 65),cd4_diag_sq=c(7.071068, 14.14214),histo_y=c(2010, 2020), rna_diag_log=c(5,6)))



cph_summ <- rbind(summ1[9:10,], summ1[11:12, ], summ1[1:2,],summ2[2,], summ3[2,], summ4[2, ], 
                  summ1[3:4, ], summ2[4,], summ3[4, ], summ4[4, ],
                  summ1[5:6, ], summ2[6, ], summ3[6,], summ4[6,],   
                  summ1[7:8, ], summ2[8, ], summ3[8,], summ4[8,],summ1[13:14,])

cph_summ$Variable <- c("Sex (ref: Female)", "Male", 
                       "Migratied from Histo endemic country(ref=no)", "Yes",
                       "Age at diagnosis (ref: 35)", "25",  "45", "55", "65",
                       "Histo diagnosis year (ref: 2010)", "2000", "2005", "2015","2020",
                       "CD4 at histo diagnosis (ref: 50)", "10" ," 25", "100", "200",
                       "log HIV RNA at diagnosis (ref: 5)", "1","3","7","6","ART initiation status wrt histo ref:pre ART at histo","post ART at histo")

names(cph_summ)[names(cph_summ) %in% grep("Lower", names(cph_summ), value=TRUE)] <- "Lower.0.95"
names(cph_summ)[names(cph_summ) %in% grep("Upper", names(cph_summ), value=TRUE)] <- "Upper.0.95"

cph_summ$CI <- with(cph_summ, paste0("(", formatC(Lower.0.95, digits=3, format= "f"), ", ", formatC(Upper.0.95, digits=3, format= "f"), ")"), sep="")

cph_summ[which(cph_summ$Type==1), c('Effect', 'Lower.0.95', 'Upper.0.95', 'CI')] <- NA
cph_summ <- cph_summ[, c("Variable", "Effect", "CI")]
rownames(cph_summ) <- NULL

mymod_aov <- as.data.frame(anova(mymod))
mymod_aov <- mymod_aov[c(1,3,2,4:7), ]
mymod_aov$Variable <- c(grep("ref", cph_summ$Variable, value=TRUE))
mymod_aov <- mymod_aov[, c("Variable", "P")]
names(mymod_aov)[2] <- "p"

cph_summ <- plyr::join(cph_summ, mymod_aov, by="Variable")
cph_summ$p <- ifelse(cph_summ$p< 0.0001, "<0.0001", round(cph_summ$p, digits=3))

label(c$histo_y) <- "Histoplasmosis"
label(la_c$histo_y) <- "Histoplasmosis"
label(tn_c$histo_y) <- "Histoplasmosis"
label(tn_c$cd4_diag_sq)  <- "CD4 at daignosis (sq)"
label(la_c$cd4_diag_sq)  <- "CD4 at daignosis (sq)"
label(tn_c$rna_diag_log) <- "rna at daignosis (log)"
label(la_c$rna_diag_log) <- "rna at daignosis (log)"

save(cph_summ, file="cph_summ.Rdata")
save(c, file="c.Rdata")
save(la_c, file="la_c.Rdata")
save(tn_c, file="tn_c.Rdata")

## Proportional Hazards model with imputation overall

set.seed(5)

surv.obj <- with(c, Surv(time, death_y==1))
cph_multi_i <- aregImpute(~ (time*death_y) + sex + age_at_diag + mig_ende + histo_y + cd4_diag_sq + rna_diag_log ,group= c$site,data=c, n.impute = 20)

dd <- datadist(c)
options(datadist= 'dd')


dd$limits$age_at_diag[2] <- 35
dd$limits$sex[2] <- "Female"
dd$limits$histo_y[2] <- 2010
dd$limits$cd4_diag_sq[2] <- 50
dd$limits$rna_diag_log[2] <- 5 
dd$limits$art_status[2] <- "pre ART at histo"
dd$limits$mig_ende[2] <- "No"

cph_multi <- fit.mult.impute(surv.obj ~ + sex + age_at_diag + mig_ende + histo_y + cd4_diag_sq + rna_diag_log + art_status + strat(site), fitter=cph,
                             data = c,fit.reps=TRUE,xtrans=cph_multi_i,n.impute=20,fitargs=list(x=TRUE, y=TRUE))


summ1 <- as.data.frame(summary(cph_multi, age_at_diag=c(35, 25),cd4_diag_sq=c(7.071068, 3.162278),histo_y=c(2010, 2000), rna_diag_log=c(5,1)))
summ2 <- as.data.frame(summary(cph_multi, age_at_diag=c(35, 45),cd4_diag_sq=c(7.071068, 5),histo_y=c(2010, 2005), rna_diag_log=c(5,3)))
summ3 <- as.data.frame(summary(cph_multi, age_at_diag=c(35, 55),cd4_diag_sq=c(7.071068, 10),histo_y=c(2010, 2015), rna_diag_log=c(5,7)))
summ4 <- as.data.frame(summary(cph_multi, age_at_diag=c(35, 65),cd4_diag_sq=c(7.071068, 14.14214),histo_y=c(2010, 2020), rna_diag_log=c(5,6)))


cph_i <- rbind(summ1[9:10,], summ1[11:12, ], summ1[1:2,],summ2[2,], summ3[2,], summ4[2, ], 
               summ1[3:4, ], summ2[4,], summ3[4, ], summ4[4, ],
               summ1[5:6, ], summ2[6, ], summ3[6,], summ4[6,],   
               summ1[7:8, ], summ2[8, ], summ3[8,], summ4[8,],summ1[13:14,])

cph_i$Variable <- c("Sex (ref: Female)", "Male", 
                    "Migratied from Histo endemic country(ref=no)", "Yes",
                    "Age at diagnosis (ref: 35)", "25",  "45", "55", "65",
                    "Histo diagnosis year (ref: 2010)", "2000", "2005", "2015","2020",
                    "CD4 at histo diagnosis (ref: 50)", "10" ," 25", "100", "200",
                    "log HIV RNA at diagnosis (ref: 5)", "1","3","7","6","ART initiation status wrt histo ref:pre ART at histo","post ART at histo")

names(cph_i)[names(cph_i) %in% grep("Lower", names(cph_i), value=TRUE)] <- "Lower.0.95"
names(cph_i)[names(cph_i) %in% grep("Upper", names(cph_i), value=TRUE)] <- "Upper.0.95"

cph_i$CI <- with(cph_i, paste0("(", formatC(Lower.0.95, digits=3, format= "f"), ", ", formatC(Upper.0.95, digits=3, format= "f"), ")"), sep="")

cph_i[which(cph_i$Type==1), c('Effect', 'Lower.0.95', 'Upper.0.95', 'CI')] <- NA
cph_i <- cph_i[, c("Variable", "Effect", "CI")]
rownames(cph_i) <- NULL

cph_multi_i_aov <- as.data.frame(anova(cph_multi))
cph_multi_i_aov <- cph_multi_i_aov[c(1,3,2,4:7), ]
cph_multi_i_aov$Variable <- c(grep("ref", cph_i$Variable, value=TRUE))
cph_multi_i_aov <- cph_multi_i_aov[, c("Variable", "P")]
names(cph_multi_i_aov)[2] <- "p"

cph_i <- plyr::join(cph_i, cph_multi_i_aov, by="Variable")
cph_i$p <- ifelse(cph_i$p< 0.0001, "<0.0001", round(cph_i$p, digits=3))


save(cph_i, file="cph_i.Rdata")

## Proportional Hazards model with imputation with ccasanet

set.seed(5)

surv.obj <- with(la_c, Surv(time, death_y==1))
cph_multi_i <- aregImpute(~ (time*death_y) + sex + age_at_diag + mig_ende + histo_y + cd4_diag_sq + rna_diag_log + art_status,group= la_c$site,data=la_c, n.impute = 20)

dd <- datadist(la_c)
options(datadist= 'dd')


dd$limits$age_at_diag[2] <- 35
dd$limits$sex[2] <- "Female"
dd$limits$histo_y[2] <- 2010
dd$limits$cd4_diag_sq[2] <- 50
dd$limits$rna_diag_log[2] <- 5 
dd$limits$art_status[2] <- "pre ART at histo"
dd$limits$mig_ende[2] <- "No"

cph_multi <- fit.mult.impute(surv.obj ~ + sex + age_at_diag + mig_ende + histo_y + cd4_diag_sq + rna_diag_log  + art_status+ strat(site), fitter=cph,
                             data=la_c,fit.reps=TRUE,xtrans=cph_multi_i,n.impute=20,fitargs=list(x=TRUE, y=TRUE))


summ1 <- as.data.frame(summary(cph_multi, age_at_diag=c(35, 25),cd4_diag_sq=c(7.071068, 3.162278),histo_y=c(2010, 2000), rna_diag_log=c(5,1)))
summ2 <- as.data.frame(summary(cph_multi, age_at_diag=c(35, 45),cd4_diag_sq=c(7.071068, 5),histo_y=c(2010, 2005), rna_diag_log=c(5,3)))
summ3 <- as.data.frame(summary(cph_multi, age_at_diag=c(35, 55),cd4_diag_sq=c(7.071068, 9),histo_y=c(2010, 2015), rna_diag_log=c(5,7)))
summ4 <- as.data.frame(summary(cph_multi, age_at_diag=c(35, 65),cd4_diag_sq=c(7.071068, 14.14214),histo_y=c(2010, 2020), rna_diag_log=c(5,6)))


cph_icc <- rbind(summ1[9:10,], summ1[11:12, ], summ1[1:2,],summ2[2,], summ3[2,], summ4[2, ], 
                 summ1[3:4, ], summ2[4,], summ3[4, ], summ4[4, ],
                 summ1[5:6, ], summ2[6, ], summ3[6,], summ4[6,],   
                 summ1[7:8, ], summ2[8, ], summ3[8,], summ4[8,],summ1[13:14,])

cph_icc$Variable <- c("Sex (ref: Female)", "Male", 
                      "Migratied from Histo endemic country(ref=no)", "Yes",
                      "Age at diagnosis (ref: 35)", "25",  "45", "55", "65",
                      "Histo diagnosis year (ref: 2010)", "2000", "2005", "2015","2020",
                      "CD4 at histo diagnosis (ref: 50)", "10" ," 25", "100", "200",
                      "log HIV RNA at diagnosis (ref: 5)", "1","3","7","6","ART initiation status wrt histo ref:pre ART at histo","post ART at histo")

names(cph_icc)[names(cph_icc) %in% grep("Lower", names(cph_icc), value=TRUE)] <- "Lower.0.95"
names(cph_icc)[names(cph_icc) %in% grep("Upper", names(cph_icc), value=TRUE)] <- "Upper.0.95"

cph_icc$CI <- with(cph_icc, paste0("(", formatC(Lower.0.95, digits=3, format= "f"), ", ", formatC(Upper.0.95, digits=3, format= "f"), ")"), sep="")

cph_icc[which(cph_icc$Type==1), c('Effect', 'Lower.0.95', 'Upper.0.95', 'CI')] <- NA
cph_icc <- cph_icc[, c("Variable", "Effect", "CI")]
rownames(cph_icc) <- NULL

cph_multi_i_aov <- as.data.frame(anova(cph_multi))
cph_multi_i_aov <- cph_multi_i_aov[c(1,3,2,4:7), ]
cph_multi_i_aov$Variable <- c(grep("ref", cph_icc$Variable, value=TRUE))
cph_multi_i_aov <- cph_multi_i_aov[, c("Variable", "P")]
names(cph_multi_i_aov)[2] <- "p"

cph_icc <- plyr::join(cph_icc, cph_multi_i_aov, by="Variable")
cph_icc$p <- ifelse(cph_icc$p< 0.0001, "<0.0001", round(cph_icc$p, digits=3))


save(cph_icc, file="cph_icc.Rdata")