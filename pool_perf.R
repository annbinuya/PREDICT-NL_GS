#Pooling functions and other functions
#Author: Mary Ann Binuya
#Last updated: December 21, 2023

#1. Pool estimates and produce 95% confidence intervals:
  #est <- array of (performance) estimates
  #se <- array of standard errors
  pool_estimates <- function(est, se, logit_trans = FALSE, log_trans = FALSE, exp_trans = FALSE){
    RR_se <- function(est, se){
        m <- length(est)
        w_var <- mean(se^2) # within variance
        b_var <- var(est) # between variance
        t_var <- w_var + (1 + (1/m)) * b_var # total variance
        se_total <- sqrt(t_var) # total se
        r <- (1 + 1 / m) * (b_var / w_var)
        df <- (m - 1) * (1 + (1/r))^2 # degrees of freedom
        t <- qt(0.975, df) # inverse cdf of t dist (~N dist for large df)
        res <- c(se_total, t)
        return(res)
      }
      
      est <- unlist(est)
      se <- unlist(se)
      
      if(logit_trans){ #if logit transformation requested
        est_log <- log(est/(1-est)) 
        se_log <- se / (est * (1-est)) #for formula, see supplementary page 9 of Debray et al, 2017
        se_total <- RR_se(est_log, se_log) #RR pooling after logit transformation
        
        inv.est <- round(exp(mean(est_log))/(1+exp(mean(est_log))), 2) #back-transform
        inv.est.u <- round(exp(mean(est_log) + (se_total[2]*se_total[1])) /
          (1 + exp(mean(est_log) + (se_total[2]*se_total[1]))), 2)
        inv.est.l <- round(exp(mean(est_log) - (se_total[2]*se_total[1])) /
          (1 + exp(mean(est_log) - (se_total[2]*se_total[1]))), 2)
        res <- paste0(inv.est, " (", inv.est.l, ", ", inv.est.u, ")")
        
      } else if (exp_trans) { #if exp transformation requested for coefs (i.e., beta to hazard ratios)
        se_total <- RR_se(est, se) #RR pooling before exp transformation
        
        exp.est <- round(exp(mean(est)), 1)
        exp.est.u <- round(exp(mean(est) + (se_total[2]*se_total[1])), 1)
        exp.est.l <- round(exp(mean(est) - (se_total[2]*se_total[1])), 1)
      
        res <- paste0(exp.est, " (", exp.est.l, ", ", exp.est.u, ")")
        
      } else if (log_trans) { #if log transformation requested (e.g., for O/E ratios)
        est_log <- log(est) 
        se_log <- se #for O/E ratio, "se" here is se(obs risk)/obs risk (see supplementary page 9 of Debray et al, 2017); otherwise se(obs risk) / exp_risk
        se_total <- RR_se(est_log, se_log) # RR pooling after log transformation
        
        inv.est <- round(exp(mean(est_log)), 1) #back-transform
        inv.est.u <- round(exp(mean(est_log) + (se_total[2]*se_total[1])), 1)
        inv.est.l <- round(exp(mean(est_log) - (se_total[2]*se_total[1])), 1)
        res <- paste0(inv.est, " (", inv.est.l, ", ", inv.est.u, ")")
        
      } else {
        mean.est <- round(mean(est), 1)
        se_total <- RR_se(est, se)
        mean.est.u <- round(mean(est) + (se_total[2]*se_total[1]), 1)
        mean.est.l <- round(mean(est) - (se_total[2]*se_total[1]), 1)
        
        res <- paste0(mean.est, " (", mean.est.l, ", ", mean.est.u, ")")
      }
      return(res)
  }

#2. Create a calibration plot:
  #data <- data derived from pool_cal function
  #times <- calibration at t years
  #main <- plot title
  #AUC/C <- AUC/C-index
  #calslope <- calibration slope
  #OEratio <- OE ratio
  #size_legend <- size of texts
  #size_bintext <- size of texts next to histograms
  #line_bins <- position of histograms
  #triangles <- quantile means
  #g <- number of groups to define quantiles

  calplot_fn <- function(data, tmax, main = "", AUC, calslope, OEratio, limit = 0.5, size_lab = 1, size_legend = 0.45, size_bintext = 0.5, line_bins = 0, triangles = FALSE, g) {
    df <- data
    
    # Predicted risks
    df$x <- df$predicted_risk #predicted probabilities at tmax
    
    df$x.ll <- log(-log(1 - df$x)) #complementary log-log transformation of the predicted survival; improves linearity and lessens # of knots needed for rcs per Austin, et al, 2020
    
    model <- cph(Surv(time, event) ~ rcs(x.ll, 5), data = df, x = TRUE, y = TRUE, surv = TRUE)
    
    # Observed proportions
    xx <- seq(quantile(df$x, prob = 0.01), quantile(df$x, prob = 0.99), length = 100)
    xx.ll <- log(-log(1 - xx))
    xx.ll.df <- data.frame(x.ll = xx.ll)
    
    y <- 1 - survest(model, newdata = xx.ll.df, times = tmax)$surv
    y.lower <- 1 - survest(model, newdata = xx.ll.df, times = tmax)$upper
    y.upper <- 1 - survest(model, newdata = xx.ll.df, times = tmax)$lower
    
    # Plot parameters
    xlim <- c(0, limit + 0.01)
    ylim <- c(-0.04, limit + 0.01)
    xlab <- "Predicted probability"
    ylab <- "Observed proportion"
    
    # Plot
    par(las = 1, xaxs = "i", yaxs = "i")
    plot(0, 0, type = "n",
         xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = ylab,
         main = main, cex.main = 1,
         cex.lab = size_lab, cex.axis = size_lab, lwd = 2)
    polygon(c(xx, rev(xx), xx[1]),
            c(y.lower, rev(y.upper), y.lower[1]),
            border = NA, density = 50, angle = -20, col = "gray")
    abline(coef = c(0,1), lty = 1, col = "gray") #diagonal line
    lines(xx, y, lwd = 2, col = "black")
    
    # Triangles
    if (triangles) {
      q <- Hmisc::cut2(df$x, levels.mean = TRUE, g=g) #group predicted risks
      means <- as.double(levels(q))
      y1 <- 1-survest(model, newdata=df, times = tmax)$surv
      prop <- tapply(y1, q, mean) #mean Observed proportions
      points(means, prop, pch=17, cex=1, col="black") #triangles
    }
    
    # Histograms:
    #rect(0,-0.025,1,0,lwd=0,border=NA)
    lim <- c(min(df$x), max(df$x))
    bins <- seq(lim[1], lim[2], length=101) #bins
    
    f0	<- table(cut(df$x[df$event == 0], bins)) #non-events
    f1	<- table(cut(df$x[df$event == 1], bins)) #events
    
    bins0 <- (bins[-101])[f0 > 0] #x bins, non-events
    bins1 <- (bins[-101])[f1 > 0] #x bins, events
    
    pcty1 <- as.numeric(f1/sum(f1))*0.1 #y bin, fix height multiplier for now
    pctx1 <- rep(0, length(pcty1))
    for (i in 1:length(pcty1)) { #histograms
      if (pcty1[i]>0) 
        rect(bins1[i], line_bins,
             bins1[i+1], line_bins + pcty1[i],
             col="gray", border=NA)
    }
    
    pcty0 <- as.numeric(f0/sum(f0))*0.1 #y bin, fix height multiplier for now
    pctx0 <- rep(0, length(pcty0))
    for (i in 1:length(pcty0)) { #histograms
      if (pcty0[i]>0)
        rect(bins0[i], line_bins,
             bins0[i+1], line_bins - pcty0[i],
             col="gray", lwd = 0.5, border=NA)
    }
    
    abline(h = line_bins, lty = 3)
    text(x = limit + 0.01, y = (line_bins + 0.005),
         labels = "Events",
         cex = size_bintext, pos = 2, col = "darkgray")
    text(x = limit + 0.01, y = (line_bins - 0.005),
         labels = "Non-events",
         cex = size_bintext, pos = 2, col = "darkgray")
    
    #Texts
    text(x=0, y=(limit), labels=paste("Discrimination"), cex=size_legend, pos=4)
    text(x=0, y=(limit-0.01), labels=paste("...AUC: ", AUC, sep=""), cex=size_legend, pos=4)
    text(x=0, y=(limit-0.02), labels=paste("Calibration "), cex=size_legend, pos=4)
    text(x=0, y=(limit-0.03), labels=paste("...Slope: ", calslope, sep=""), cex=size_legend, pos=4)
    text(x=0, y=(limit-0.04), labels=paste("...O/E ratio: ", OEratio, sep=""), cex=size_legend, pos=4)
  }
  
  
#3. Plot 2 calibration plots
  calplot2df <- function(data1, data2, tmax, main, limit, size_lab = 1, size_legend = 0.5, triangles = FALSE, g, AUCdif){
    #data1 = risk_df file without grisk
    #data2 = risk_df file with grisk
    
    # Predicted risks for data1
    data1$x1 <- data1$predicted_risk
    data1$x1.ll <- log(-log(1 - data1$x1))
    
    # Predicted risks for data2
    data2$x2 <- data2$predicted_risk
    data2$x2.ll <- log(-log(1 - data2$x2))
    
    model1 <- cph(Surv(time, event) ~ rcs(x1.ll, 5), data = data1, x = TRUE, y = TRUE, surv = TRUE)
    model2 <- cph(Surv(time, event) ~ rcs(x2.ll, 5), data = data2, x = TRUE, y = TRUE, surv = TRUE)
    
    # Observed risks for data1
    xx1 <- seq(quantile(data1$x1, prob = 0.01), quantile(data1$x1, prob = 0.99), length = 100)
    xx1.ll <- log(-log(1 - xx1))
    xx1.ll.df <- data.frame(x1.ll = xx1.ll)
    y1 <- 1 - survest(model1, newdata = xx1.ll.df, times = tmax)$surv
    y1.lower <- 1 - survest(model1, newdata = xx1.ll.df, times = tmax)$upper
    y1.upper <- 1 - survest(model1, newdata = xx1.ll.df, times = tmax)$lower
    
    # Observed risks for data2
    xx2 <- seq(quantile(data2$x2, prob = 0.01), quantile(data2$x2, prob = 0.99), length = 100)
    xx2.ll <- log(-log(1 - xx2))
    xx2.ll.df <- data.frame(x2.ll = xx2.ll)
    y2 <- 1 - survest(model2, newdata = xx2.ll.df, times = tmax)$surv
    y2.lower <- 1 - survest(model2, newdata = xx2.ll.df, times = tmax)$upper
    y2.upper <- 1 - survest(model2, newdata = xx2.ll.df, times = tmax)$lower
    
    # Plot parameters
    xlim <- c(0, limit+0.01)
    ylim <- c(0, limit+0.01)
    xlab <- "Predicted probability"
    ylab <- "Observed proportion"
    
    # Plot
    par(las = 1, xaxs = "i", yaxs = "i")
    plot(0, 0, type = "n",
         xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = ylab,
         main = main, cex.main = 1,
         cex.lab = size_lab, cex.axis = size_lab, lwd = 2)
    polygon(c(xx1, rev(xx1), xx1[1]),
            c(y1.lower, rev(y1.upper), y1.lower[1]),
            border = NA, density = 50, angle = -20, col = "light gray")
    polygon(c(xx2, rev(xx2), xx2[1]),
            c(y2.lower, rev(y2.upper), y2.lower[1]),
            border = NA, density = 50, angle = 20, col = "light gray")
    abline(coef = c(0, 1), lty = 1, col = "gray") #diagonal line
    lines(xx1, y1, lwd = 2, col = "black")
    lines(xx2, y2, lty = 2, lwd = 2, col = "black")  # Add a second line (data2) in blue
    
    # Triangles
    if (triangles) {
      q1 <- Hmisc::cut2(data1$x1, levels.mean = TRUE, g = g)
      q2 <- Hmisc::cut2(data2$x2, levels.mean = TRUE, g = g)
      means1 <- as.double(levels(q1))
      means2 <- as.double(levels(q2))
      y1 <- 1-survest(model1, newdata=data1, times = tmax)$surv
      y2 <- 1-survest(model2, newdata=data2, times = tmax)$surv
      prop1 <- tapply(y1, q1, mean)
      prop2 <- tapply(y2, q2, mean)
      points(means1, prop1, pch = 17, cex = 1, col = "black")
      points(means2, prop2, pch = 17, cex = 1, col = "black")
    }
    
    # Legend
    legend("topleft", legend = c("Model without MammaPrint", "Model with MammaPrint"),
           col = c("black", "black"), lty = c(1, 2), lwd = 2, cex = size_legend, bty = "n")
    
    # Text
    text(x=limit-0.03, y=0.005, labels=paste("AUC difference: ", AUCdif, sep=""), cex=0.8, pos=4)
  }
  
  
#4. Calculate net benefit (NB):
  #In survival settings, we can manually calculate NB over a range of thresholds using this formula (Vickers et al, 2008):
  #NB = TP/n -  w x FP/n, where:
  #w= pt/(1-pt), where that pt = threshold probability
  #TP = [1 - (S(t) | x = 1)] * P(x = 1) * n,
  #FP = (S(t) | x = 1) * P(x = 1) * n,
  #(x = 1) = if a patient has a predicted risk from the model b	% pt
  #S(t) = the Kaplan-Meier survival probability at time t
  #n = sample size
  
  #data <- data derived from pool_cal function
  
  nb_fn_TP <- function(data, tmax, thresholdmax = 1) {
    
    thresholds <- seq(0.01, thresholdmax, by = 0.001) #use by=0.001 due to low risks and more homogenous case mix; otherwise use 0.01
    
    NB_f <- lapply(thresholds, function(pt) {
      
      #NB treat all
      m_all <- 1 - summary(survfit(Surv(time, event) ~ 1, data = data), times = tmax)$surv #BC deaths
      NB_all <- m_all - (1-m_all) * (pt/(1-pt))
      
      #NB for model
      prop_pred <- nrow(subset(data, predicted_risk >= pt))/nrow(data)  #proportion predicted high risk from model; alternatively, mean(predrisk>=pt)
      
      surv <- try(
        summary(survfit(Surv(time, event) ~ 1, data = data[data$predicted_risk >= pt, ]), 
                times = tmax), silent = TRUE)
      
      if (class(surv) == "try-error") {
        TP <- 0
        FP <- 0
        NB <- 0 #no observations above threshold
      } else {
        m_model <- 1 - surv$surv
        TP <- m_model * prop_pred
        FP <- (1 - m_model) * prop_pred
        NB <- TP - FP * (pt/(1-pt))
      }
      
      NBres <- data.frame("threshold" = pt,
                          "NB_all" = NB_all,
                          "TP_all" = m_all,
                          "FP_all" = (1-m_all),
                          "NB" = NB,
                          "TP" = TP,
                          "FP" = FP)
    })
    
    #Bind results
    NB_res <- do.call(rbind, NB_f)
    
    return(NB_res)
  }
  
  nb_fn_TN <- function(data, tmax, thresholdmax = 1) { #note: labels can be confusing, refer to formula in manuscript
    
    thresholds <- seq(0.01, thresholdmax, by = 0.001)
    
    NB_f <- lapply(thresholds, function(pt) {
      
      #NB treat none
      m_none <- 1 - summary(survfit(Surv(time, event) ~ 1, data = data), times = tmax)$surv #BC deaths
      NB_none <- (1-m_none) - m_none * ((1-pt)/pt) #(1-m_none) is just surv prob estimate at tmax, weight is reversed
      
      #NB for model
      prop_pred <- nrow(subset(data, predicted_risk < pt))/nrow(data)  #proportion predicted low risks/below threshold
      
      surv <- try(
        summary(survfit(Surv(time, event) ~ 1, data = data[data$predicted_risk < pt, ]), 
                times = tmax), silent = TRUE)
      
      if (class(surv) == "try-error") {
        TN <- 1
        FN <- 1
        NB <- 1 #no observations below threshold
      } else {
        m_model <- surv$surv
        TN <- m_model * prop_pred
        FN <- (1 - m_model) * prop_pred
        NB <- TN - FN * ((1-pt)/pt)
      }
      
      NBres <- data.frame("threshold" = pt,
                          "NB_none" = NB_none,
                          "TN_none" = (1-m_none),
                          "FN_none" = m_none,
                          "NB" = NB,
                          "TN" = TN,
                          "FN" = FN)
    })
    
    #Bind results
    NB_res <- do.call(rbind, NB_f)
    
    return(NB_res)
  }