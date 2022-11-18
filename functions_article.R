
#  Plots
#---------


#Function for my km curves using ggsurvplot

my_ggkmplot <- function(myfit, title, labels, conf.int=F, fun="event", xlim=c(15,50), isFive=T, break.time.by=10, legend="right"){
  
  palette5 = c("#92C5DE","#0571B0","#756BB1","#CA0020", "#F4A582")
  palette3 = c("#0571B0","#756BB1","#CA0020")
  if(isFive){palette=palette5
  }else{palette=palette3}
  
  
  ggsurvplot(myfit,
             fun = fun, conf.int = conf.int, censor = F,
             title = title,  size = 0.5, xlim=xlim, break.time.by=break.time.by,
             palette = palette, ggtheme = theme_bw(), 
             legend = legend, legend.title = "PRS",  legend.labs=labels)
}




# Cox plot using ggsurvplot from survminer package

my_ggcoxplot <- function(fit, new_df, title, labels, conf.int=F, xlim=c(15,45), ylim=c(0,0.2), legend="right", 
                         ylab="Cumulative event", xlab="Years", break.time.by=10, break.y.by=10, isFive=T, legend.title="PRS"){
  
  labels=labels
  palette5 = c("#92C5DE","#0571B0","#756BB1","#CA0020", "#F4A582")
  palette3 = c("#0571B0","#756BB1","#CA0020")
  if(isFive){palette=palette5
  }else{palette=palette3}
  
  ggsurvplot(fit, data = new_df,
             fun = "event", conf.int = conf.int, censor = F, 
             xlim=xlim, break.time.by = break.time.by, ylim=ylim, break.y.by=break.y.by,
             title = title,  size = 0.5, 
             palette = palette, ggtheme = theme_bw(), 
             legend = legend, legend.title = legend.title,  legend.labs=labels,
             xlab = xlab, ylab = ylab, 
             font.tickslab=9.5, font.legend=9.5, surv.scale="percent")
  
}



my_expand_covs <- function(df, batch) {
  
  expand.grid(
    batch = batch,
    BL_YEAR = round(mean(df$BL_YEAR)),
    PC1 = mean(df$PC1),
    PC2 = mean(df$PC2),
    PC3 = mean(df$PC3),
    PC4 = mean(df$PC4),
    PC5 = mean(df$PC5),
    PC6 = mean(df$PC6),
    PC7 = mean(df$PC7),
    PC8 = mean(df$PC8),
    PC9 = mean(df$PC9),
    PC10 = mean(df$PC10)
  ) 
}


# Tables
#-------

# kableExtra wrapper

my.kable <- function(my.table) {
  my.table %>%
    knitr::kable() %>%
    kable_classic() %>%
    row_spec(0,bold=TRUE)
  
}



# Extract coefficients and ci's

my_extr_coef <- function(cx, title, select){
  summary(cx)$coefficients[,c(2,5)] %>%
    cbind(row.names(.),.,summary(cx)$conf[,3:4]) %>%
    as.tibble() %>%
    rename(names=V1, est="exp(coef)", pval="Pr(>|z|)",  lower="lower .95", upper="upper .95") %>%
    filter(str_detect(names, select)) %>%
    mutate(name=title)
}




# Extracts number of cases and controls

extract_ns2 <- function(endpoint, score_name, df) {
  
  df %>%
    select(all_of(score_name), all_of(endpoint)) %>%
    rename_(quantile = "score_name", endpoint = "endpoint") %>%
    drop_na %>%
    group_by(quantile, endpoint) %>%
    summarise(n=n()) %>%
    spread(endpoint, n) %>%
    rename(cases = "1", controls="0") %>%
    ungroup() %>%
    mutate(endpoint = endpoint,
           score_name = score_name) %>%
    select(endpoint, score_name, quantile, controls, cases) #%>%
    #add_row(endpoint=endpoint, score_name=score_name, quantile=title, controls=NA, cases=NA, .before = 1)

}




# Combines estimates and ci's to same column, rounds some values, replaces NA's etc. 
#   Input table 'res.df' must have columns "est", "lower", "upper", "pval"

my_tidy_table <- function(res.df, est_repl = " - ", p_repl = " - ", est_dec = 2){
 
  res.df %>%
    mutate_at(c("est", "lower", "upper", "pval"),  as.numeric) %>%
    mutate_at(c("est", "lower", "upper"), round, est_dec) %>%
    mutate_at("pval", signif, 2) %>%
    mutate(est = str_glue("{est} ({lower}-{upper})"),
           est = str_replace(est, "^N.*", est_repl),
           pval = if_else(pval < 1e-300, 0, pval), 
           pval = as.character(pval),
           pval = replace_na(pval, p_repl),
           pval = str_replace(pval, "^0$", "<1e-300")
    )%>%
    mutate_at(c("est", "lower", "upper", "pval"),  as.character) %>%
    mutate_at(c("est", "lower", "upper", "pval"),  replace_na, " - ") %>%
    select(-"lower",-"upper")
  
}



#my_cxlist_to_hrtable2(cxs.cs, list.df.p, ep_title="Endpoint", select="SBP_") 
my_cxlist_to_hrtable2 <- function(cxlist, df.list, ep_title="Endpoint", select="SBP_"){
  
  #Input: list of cx-results
  #Requires my_exr_coef() and my_tidy_table()
  
  names <- as.list(names(cxlist))
  
  #  tmp_list <- lapply(names, function(name){
  #    my_extr_coef(cxlist[[name]], title=name , select=paste0(name,"_") )})
  tmp_list <- lapply(names, function(name){
    my_extr_coef(cxlist[[name]], title=name , select=select)}) %>%
    setNames(names)
  
  #n's picked
  ns_list <- lapply(names, function(name){
    
    ns_tmp <- df.list[[name]] %>% 
      group_by(get(name)) %>% 
      summarise(n=n()) %>% 
      t() %>%  as_tibble %>%
      slice(-1)
    ns_tmp <- ns_tmp[,1:2]
    names(ns_tmp) <- c("Controls", "Cases")
    ns_tmp
    #rename(Controls = `1`, Cases = `2`) %>%   #Old version
    #select("Controls", "Cases")               #Stoped working at rmarkdown... 
  }) %>% setNames(names)                       #It has "V1" and "V2" as colnames, while at actual R they are `1` and `2`
  
  tmp_list2 <- lapply(names, function(name){
    bind_cols(tmp_list[[name]], ns_list[[name]])
  })%>% setNames(names)
  
  bind_rows(tmp_list2) %>%
    my_tidy_table() %>%
    select(name, est, pval, Cases, Controls) %>%
    rename(!!ep_title := name, `HR (95% CI)` = est, `P-value` = pval)
}

#Function from Felix; prettifies estimates

prettify_estimate <- function(estimate, stderr, r) {
  rounded_estimate <- signif(estimate, r)
  rounded_ci_low   <- signif(estimate - 1.96 * stderr, r)
  rounded_ci_high  <- signif(estimate + 1.96 * stderr, r)
  pretty_string    <- paste0(rounded_estimate, " (", rounded_ci_low, "-", rounded_ci_high, ")")
  return(pretty_string)
}

# Extracts coefficients from glm results

extract_glm_table <- function(myfit){ 
  
  cbind(summary(myfit)$coef, confint(myfit)) %>%
    as.data.frame() %>%
    rownames_to_column(var = "Variable") %>%
    rename(Low = "2.5 %", High = "97.5 %", "P-value" = "Pr(>|z|)") %>%
    mutate_at(vars(Estimate, Low, High), ~round(exp(.),2)) %>%
    mutate(`P-value` = format(signif(`P-value`,2), scientific = TRUE),
           Odds = str_glue("{Estimate} ({Low}-{High})")) %>% 
    filter(Variable != "(Intercept)") %>%
    select(Variable, Odds, `P-value`)
  
}



#Self written functions to generate 'table one'
#There would be plenty of packages, but this was written when packages where limited

#This is unnecessarily complicated when no vector, but only one value, but idea
#is to keep it systematically similar in both cases.

binvar_to_chartable <- function(bin_vars, case, df) {
  
  tmp <-  df %>%
    select(eval(case), all_of(bin_vars)) %>%
    pivot_longer(cols = all_of(bin_vars), names_to = "Characteristics") %>%
    group_by_at(c(case, "Characteristics")) %>%
    summarise(n=sum(value, na.rm=T), perc=mean(value, na.rm=T)*100) %>%
    pivot_wider(names_from = case, values_from = c(n, perc)) %>%
    mutate(`Relative difference, (%)` = ((perc_1 - perc_0)/perc_0)*100 )
  
  as_tibble(bin_vars) %>%   #To keep given order of variables
    rename(Characteristics = value) %>%
    left_join(tmp, .by="Characteristics") %>%
    mutate_if(is.numeric, round, 2) %>%
    mutate(Case = str_glue("{n_1} ({perc_1})"),
           Control = str_glue("{n_0} ({perc_0})")) %>%
    select(Characteristics, Case, Control) %>%
    mutate(Characteristics = str_glue("{Characteristics}, n (%)"))
  
}

contvar_to_chartable <- function(cont_vars, case,  df) {
  
  tmp <- df %>% 
    select(eval(case), all_of(cont_vars)) %>%
    pivot_longer(cols = all_of(cont_vars), names_to = "Characteristics") %>%
    group_by_at(c(case, "Characteristics")) %>%
    summarise(mean=mean(value, na.rm=T), sd=sd(value, na.rm=T)) %>%
    pivot_wider(names_from = case, values_from = c(mean, sd)) %>%
    mutate(`Relative difference, (%)` = ((mean_1 - mean_0)/mean_0)*100 )
  
  as_tibble(cont_vars) %>%   #To keep given order of variables
    rename(Characteristics = value) %>%
    left_join(tmp, .by="Characteristics") %>%
    mutate_if(is.numeric, round, 1) %>%
    mutate(Case = str_glue("{mean_1} +/- {sd_1}"), 
           Control = str_glue("{mean_0} +/- {sd_0}")) %>%
    select(Characteristics, Case, Control)  %>%
    mutate(Characteristics = str_glue("{Characteristics}, mean +/- SD"))
  
}


n_to_chartable <- function(case, df) {
  
  df %>%
    group_by_at(c(case)) %>% 
    summarise(n=n()) %>%
    ungroup()%>%
    pivot_wider(names_from = case, values_from = n) %>%
    mutate(Characteristics = str_glue("N")) %>%
    rename(Case =  `1`, Control = `0`) %>%
    mutate(Case = str_glue("{Case}"), Control = str_glue("{Control}")) %>%
    select(Characteristics, Case, Control)
  
}                



#Other functions
#---------------


my.roc.test <- function(roc1,roc2, name=""){
  
  tmp <- roc.test(roc1, roc2)
  
  cbind(cindex = tmp$estimate[1],
        est = tmp$estimate[1]-tmp$estimate[2], 
        low = tmp$conf.int[1], 
        high = tmp$conf.int[2], 
        pval = tmp$p.value ) %>%
    as.tibble() %>%
    mutate_at(vars(est, low, high), ~format(round(.,4), scientific = F) ) %>%
    #mutate_at(vars(est), ~format(round(.,4), scientific = F) ) %>%
    mutate(Name = name,
           `C-index` = round(cindex,4),
           `P-value` = format(signif(pval,2), nsmall = 1),
           Increment = str_glue("{est} ({low}-{high})")) %>%
           #Increment = str_glue("{est}")) %>% 
    select(Name, `C-index`, Increment, `P-value`) %>%
    as.data.frame()
  
}


my.reclas.wrapper <- function(reclas, name=""){
  
  cbind(nri_est = reclas$nri.cat$est,
        nri_low = reclas$nri.cat$low, 
        nri_high = reclas$nri.cat$high, 
        nri_pval = reclas$nri.cat$pval,
        idi_est = reclas$idi$est,
        idi_low = reclas$idi$low, 
        idi_high = reclas$idi$high, 
        idi_pval = reclas$idi$pval,
        recl_contr = reclas$nri.table.absent[2,3],
        recl_case = reclas$nri.table.present[1,3]) %>%
    as.tibble() %>%
    mutate_at(vars(nri_est, nri_low, nri_high, idi_est, idi_low, idi_high), ~format(round(.,4), scientific = F) ) %>%
    mutate(Name = name,
           `NRI P-value` = format(signif(nri_pval,2), nsmall = 1),
           `IDI P-value` = format(signif(idi_pval,2), nsmall = 1),
           `NRI (95% CI)` = str_glue("{nri_est} ({nri_low}-{nri_high})"),
           `Correctly reclassified cases (%)` = recl_case,
           `Correctly reclassified controls (%)` = recl_contr,
           `IDI (95% CI)` = str_glue("{idi_est} ({idi_low}-{idi_high})")) %>% 
    select(Name, `NRI (95% CI)`,`NRI P-value`, `Correctly reclassified cases (%)`, `Correctly reclassified controls (%)`, `IDI (95% CI)`, `IDI P-value`) %>%
    as.data.frame()
  
}  


#This is modified from 'predictABEL' function 'reclassification.
#Differences are in output format.
#requires 'predictABEL'-package

my.reclassification <- function (data, cOutcome, predrisk1, predrisk2, cutoff) 

  {
  c1 <- cut(predrisk1, breaks = cutoff, include.lowest = TRUE, 
            right = FALSE)
  c2 <- cut(predrisk2, breaks = cutoff, include.lowest = TRUE, 
            right = FALSE)
  tabReclas <- table(`Initial Model` = c1, `Updated Model` = c2)
  cat(" _________________________________________\n")
  cat(" \n     Reclassification table    \n")
  cat(" _________________________________________\n")
  ta <- table(c1, c2, data[, cOutcome])
  cat("\n Outcome: absent \n  \n")
  TabAbs <- ta[, , 1]
  tab1 <- cbind(TabAbs, ` % reclassified` = round((rowSums(TabAbs) - 
                                                     diag(TabAbs))/rowSums(TabAbs), 3) * 100)
  names(dimnames(tab1)) <- c("Initial Model", "Updated Model")
  print(tab1)
  cat("\n \n Outcome: present \n  \n")
  TabPre <- ta[, , 2]
  tab2 <- cbind(TabPre, ` % reclassified` = round((rowSums(TabPre) - 
                                                     diag(TabPre))/rowSums(TabPre), 3) * 100)
  names(dimnames(tab2)) <- c("Initial Model", "Updated Model")
  print(tab2)
  cat("\n \n Combined Data \n  \n")
  Tab <- tabReclas
  tab <- cbind(Tab, ` % reclassified` = round((rowSums(Tab) - 
                                                 diag(Tab))/rowSums(Tab), 3) * 100)
  names(dimnames(tab)) <- c("Initial Model", "Updated Model")
  print(tab)
  cat(" _________________________________________\n")
  c11 <- factor(c1, levels = levels(c1), labels = c(1:length(levels(c1))))
  c22 <- factor(c2, levels = levels(c2), labels = c(1:length(levels(c2))))
  
  
  
  x <- improveProb(x1 = as.numeric(c11) * (1/(length(levels(c11)))), 
                   x2 = as.numeric(c22) * (1/(length(levels(c22)))), y = data[, 
                                                                              cOutcome])
  y <- improveProb(x1 = predrisk1, x2 = predrisk2, y = data[, 
                                                            cOutcome])
  cat("\n NRI(Categorical) [95% CI]:", round(x$nri, 4), "[", 
      round(x$nri - 1.96 * x$se.nri, 4), "-", round(x$nri + 
                                                      1.96 * x$se.nri, 4), "]", "; p-value:", round(2 * 
                                                                                                      pnorm(-abs(x$z.nri)), 5), "\n")
  cat(" NRI(Continuous) [95% CI]:", round(y$nri, 4), "[", round(y$nri - 
                                                                  1.96 * y$se.nri, 4), "-", round(y$nri + 1.96 * y$se.nri, 
                                                                                                  4), "]", "; p-value:", round(2 * pnorm(-abs(y$z.nri)), 
                                                                                                                               5), "\n")
  cat(" IDI [95% CI]:", round(y$idi, 4), "[", round(y$idi - 
                                                      1.96 * y$se.idi, 4), "-", round(y$idi + 1.96 * y$se.idi, 
                                                                                      4), "]", "; p-value:", round(2 * pnorm(-abs(y$z.idi)), 
  
  #CHANGE: Here I output everything in sensible format.                                                                                                                 
                                                                                                                                                                                                                                    5), "\n")
  nri.cat <- list(est = round(x$nri, 4), 
                  low = round(x$nri - 1.96 * x$se.nri, 4),
                  high = round(x$nri + 1.96 * x$se.nri, 4),
                  pval = signif(2 * pnorm(-abs(x$z.nri)), 3),
                  cutoff = cutoff)
  
  nri.cont <- list(est = round(y$nri, 4), 
                   low = round(y$nri - 1.96 * y$se.nri, 4),
                   high = round(y$nri + 1.96 * y$se.nri, 4),
                   pval = signif(2 * pnorm(-abs(y$z.nri)), 3))
  
  idi <- list(est = round(y$idi, 4), 
              low = round(y$idi - 1.96 * y$se.idi, 4),
              high = round(y$idi + 1.96 * y$se.idi, 4),
              pval = signif(2 * pnorm(-abs(y$z.idi)), 3))
  
  list(nri.cat = nri.cat, nri.cont = nri.cont, idi = idi, nri.table.absent = tab1, nri.table.present = tab2, nri.table = tab)
  
}


