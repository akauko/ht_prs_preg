
#  Plots
#---------


# Covariate values are set for the survival curve

my_expand <- function(my_cox, var, df, batch) {
  
  new_df <- expand.grid(
    tmp = levels(get(var, df)),
    #SEX_IMPUTED = levels(df$SEX_IMPUTED),
    SEX_IMPUTED = 0.5,
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
  ) %>%
    rename_(.dots = setNames("tmp", var)) 
}




# Cox plot using ggsurvplot from survminer package

my_ggcoxplot <- function(fit, new_df, title, conf.int=T, ylim=c(0,1), xlim=c(0,80), legend="right", ylab="Cumulative disease rate", vjust=-12, break.time=20){
  
  palette = c("#92C5DE","#0571B0","#756BB1","#CA0020", "#F4A582")
  labels= c("<2.5%","2.5-20%","20-80%","80-97.5%",">97.5%")
  
 p <- ggsurvplot(fit, data = new_df,
             fun = "event", conf.int = conf.int, censor = F, 
             ylim = ylim, xlim=xlim, break.time.by = break.time,
             title = title,  size = 0.5, 
             palette = palette, ggtheme = theme_bw(), 
             legend = legend, legend.title = "PRS",  legend.labs=labels,
             xlab = "Age, Years", ylab = ylab,
             font.tickslab=9.5, font.legend=9.5)
 
 p$plot +
   theme(plot.title = element_text(vjust = vjust, hjust = 0.054, size=12))  +
   guides(color = guide_legend(reverse = TRUE), fill = guide_legend(reverse = TRUE))

 
 
}



# Kaplan Maier plot using ggsurvplot from survminer package

my_ggkmplot <- function(myfit, title){
  
  palette = c("#92C5DE","#0571B0","#756BB1","#CA0020", "#F4A582")
  labels= c("<2.5%","2.5-20%","20-80%","80-97.5%",">97.5%")

  ggsurvplot(myfit,
             fun = "event", conf.int = TRUE, censor = F,
             title = title,  size = 0.5,
             palette = palette, ggtheme = theme_bw(), 
             legend = "right", legend.title = "PRS",  legend.labs=labels)
}


# Restricted mean survival times plotted using forestplot from forestplot package

my_rmst_plot <- function(table.text, table.means, clip, xticks, zero, xlab="Age at disease onset"){

  forestplot(table.text, table.means,
           graph.pos = 2,
           is.summary=c(TRUE,TRUE,rep(FALSE,6),TRUE,rep(FALSE,6),TRUE,rep(FALSE,6)),
           clip=clip, 
           boxsize = 0.27,
           lwd.ci = 1, col=fpColors(lines="black"),  
           #line.margin = 0.5,
           xlab=xlab,
           xticks = xticks,
           txt_gp = fpTxtGp(ticks=gpar(cex=0.9), 
                            xlab=gpar(cex=0.9), 
                            summary=gpar(cex=1),
                            label=gpar(cex=1)),
           grid=T, zero=zero
           #title="RMST"
  )

}


# Tables
#-------



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

my_tidy_table <- function(res.df, est_repl = "-", est_dec = 2){
 
  res.df %>%
    mutate_at(c("est", "lower", "upper", "pval"),  as.numeric) %>%
    mutate_at(c("est", "lower", "upper"), round, est_dec) %>%
    mutate_at("pval", signif, 1) %>%
    mutate(est = str_glue("{est} ({lower}, {upper})"),
           est = str_replace(est, "^N.*", est_repl),
           pval = if_else(pval < 1e-300, 0, pval), 
           pval = replace_na(pval, "-"),
           pval = str_replace(pval, "^0$", "<1e-300")
           #pval = str_replace(pval, "e", "x10^")
           
    )%>%
    mutate_at(c("est", "lower", "upper", "pval"),  replace_na, "-") %>%
    select(-"lower",-"upper")
  
}


# Creates hr table for given scores
#   usage: my_hr_table_cont2(list(cx.sbp.cs, cx.dbp.cs), c("SBP", "DBP"), "_SCALED")
#
# Requires:
#   'my_extr_coef()' for extracting and combining coefficients and ci's
#   'my_tidy_table()' for cleaning the table

my_hr_table_cont2  <- function(cxs, titles, select="_SCALED"){
  
  mapply(my_extr_coef, cxs, titles, select) %>%
  t %>% as.tibble %>%
  my_tidy_table("1 (reference)", 2) %>% 
  mutate(name=as.character(name)) %>%
  select(name, est, pval) %>%
  rename("HR (95% CI)"=est, "P-value"=pval)
}

#cx1 <- cx.sbp.f
#cx2 <- cx.sbp.m
#endpoint <- "I9_HYPTENS"
#dfs  <- list(df_f, df_m)
#score_names <-  c("SBP_CAT", "SBP_CAT") 
#titles <- c("<b>SBP female</b>", "<b>SBP male</b>")

# Creates hr table for two risk scores of interest
#  
# Requires:
#   'my_extr_coef()' for extracting and combining coefficients and ci's
#   'my_tidy_table()' for cleaning the table
#   'extract_ns2()' to extract number of cases and controls 

my_hr_table2 <- function(cx1, cx2, endpoint="I9_HYPTENS", dfs=list(df,df), score_names = c("SBP_CAT", "DBP_CAT"), titles = c("<b>SBP</b>", "<b>DBP</b>")) {
  
  #First the number of cases and controls:
  n.1 <- extract_ns2(endpoint, score_names[1], dfs[[1]])
  n.2 <- extract_ns2(endpoint, score_names[2], dfs[[2]])

  n.all <- bind_rows(n.1, n.2) %>%
    mutate(quantile = as.character(quantile))
  
  #coefficients are extracted
  select = "_CAT"
  cx.1 <- my_extr_coef(cx1, score_names[1], select)
  cx.2 <- my_extr_coef(cx2, score_names[2], select)
  
  cx.all <- bind_rows(cx.1, cx.2) %>%
    separate(names, c("tmp", "quantile"), sep = "_CAT") %>%
    rename("score_name" = "name")
  
  #everythig combined
  n.all %>%
    left_join(cx.all, by = c("score_name", "quantile")) %>%
    my_tidy_table("1 (reference)", 2) %>% 
    mutate(ns = str_glue("{cases} / {controls}"))%>%
    select(quantile, est, pval, ns) %>%
    add_row(quantile=titles[1], est="", pval="", ns="", .before = 1) %>%
    add_row(quantile=titles[2], est="", pval="", ns="", .before = 7) %>%
    rename("PRS"=quantile, "HR (95% CI)"=est, "P-value"=pval, "Cases / Controls" = ns)
}


#cx1 <- cx.sbp.f
#cx2 <- cx.sbp.y.f
#cx3 <- cx.sbp.o.f
#endpoints <- c("I9_HYPTENS","HT_EARLY","HT_LATE")
#dfs <- list(df_f, df_f, df_f)
#score_names <-  c("SBP_CAT", "SBP_CAT", "SBP_CAT") 
#titles <-  c("<b>All</b>", "<b>Early</b>", "<b>Late</b>")

# Creates hr table for three risk scores of interest
#  
# Requires:
#   'my_extr_coef()' for extracting and combining coefficients and ci's
#   'my_tidy_table()' for cleaning the table
#   'extract_ns2()' to extract number of cases and controls 

my_hr_table3 <- function(cx1, cx2, cx3, endpoints=c("I9_HYPTENS","HT_EARLY","HT_LATE"), dfs=list(df,df,df), 
                         score_names = c("SBP_CAT", "SBP_CAT", "SBP_CAT"), 
                         titles = c("<b>All</b>", "<b>Early</b>", "<b>Late</b>")) {
  
  #First the number of cases and controls:
  n.1 <- extract_ns2(endpoints[1], score_names[1], dfs[[1]])
  n.2 <- extract_ns2(endpoints[2], score_names[2], dfs[[2]])
  n.3 <- extract_ns2(endpoints[3], score_names[3], dfs[[3]])  
  
  n.all <- bind_rows(n.1, n.2, n.3) %>%
  mutate(quantile = as.character(quantile))
  
  #coefficients are extracted
  select = "_CAT"
  cx.1 <- my_extr_coef(cx1, score_names[1], select) %>% mutate(endpoint = endpoints[1])
  cx.2 <- my_extr_coef(cx2, score_names[2], select) %>% mutate(endpoint = endpoints[2])
  cx.3 <- my_extr_coef(cx3, score_names[3], select) %>% mutate(endpoint = endpoints[3])
  
  cx.all <- bind_rows(cx.1, cx.2, cx.3) %>%
    separate(names, c("tmp", "quantile"), sep = "_CAT") %>%
    rename("score_name" = "name")

  
  #everythig combined
  test <- 
    n.all %>%
    left_join(cx.all, by = c("score_name", "quantile", "endpoint")) %>%
    my_tidy_table("1 (reference)", 2) %>% 
    mutate(ns = str_glue("{cases} / {controls}"))%>%
    select(quantile, est, pval, ns) %>%
    add_row(quantile=titles[1], est="", pval="", ns="", .before = 1) %>%
    add_row(quantile=titles[2], est="", pval="", ns="", .before = 7) %>%
    add_row(quantile=titles[3], est="", pval="", ns="", .before = 13) %>%
    rename("PRS"=quantile, "HR (95% CI)"=est, "P-value"=pval, "Cases / Controls" = ns)
}

#cx1 <- cx.sbp.i
#cx2 <- cx.sbp.y.i
#cx3 <-cx.sbp.o.i
#dfs=list(df, df, df)
#select=":SEX"

#sub("e(-\\d+)", "x10<sup>\\1</sup>", `P-value`, perl=T))

my_hr_table3_notitle <- function(cx1, cx2, cx3, endpoints=c("I9_HYPTENS","HT_EARLY","HT_LATE"), dfs=list(df,df,df), 
                         score_names = c("SBP_CAT", "SBP_CAT", "SBP_CAT"),  select = "_CAT") {
  
  #First the number of cases and controls:
  n.1 <- extract_ns2(endpoints[1], score_names[1], dfs[[1]])
  n.2 <- extract_ns2(endpoints[2], score_names[2], dfs[[2]])
  n.3 <- extract_ns2(endpoints[3], score_names[3], dfs[[3]])  
  
  n.all <- bind_rows(n.1, n.2, n.3) %>%
    mutate(quantile = as.character(quantile))
  
  #coefficients are extracted
  select = select
  cx.1 <- my_extr_coef(cx1, score_names[1], select) %>% mutate(endpoint = endpoints[1])
  cx.2 <- my_extr_coef(cx2, score_names[2], select) %>% mutate(endpoint = endpoints[2])
  cx.3 <- my_extr_coef(cx3, score_names[3], select) %>% mutate(endpoint = endpoints[3])
  
  cx.all <-
    bind_rows(cx.1, cx.2, cx.3) %>%
    mutate(quantile = str_replace(names, "^.*_CAT(.*%).*$", "\\1")) %>%
    rename("score_name" = "name")
  
  #everythig combined
  test <- 
    n.all %>%
    left_join(cx.all, by = c("score_name", "quantile", "endpoint")) %>%
    my_tidy_table("1 (reference)", 2) %>% 
    mutate(ns = str_glue("{cases} / {controls}"))%>%
    select(quantile, est, pval, ns) %>%
    rename("PRS"=quantile, "HR (95% CI)"=est, "P-value"=pval, "Cases / Controls" = ns)
}

my_rmst_coefs <- function(rmst){

  #Mean age differences extracted  
  diff.table <- lapply(rmst, function(x) x$unadjusted.result[1,]) %>%
    do.call(rbind, .) %>%
    cbind(q,.) %>%
    as.tibble() %>%
    rename(est="Est.", lower="lower .95", upper="upper .95", pval=p)
  
  #This takes average of reference values in paiwise runs
  mean.2080 <- mean(sapply(rmst, function(x) x$RMST.arm0$rmst[1]))
  
  #Mean ages extracted and all data combined and processed
  mean.table <- 
    lapply(rmst, function(x) x$RMST.arm1$rmst) %>%
    do.call(rbind, .) %>%
    cbind(q,.) %>%
    as.tibble() %>%
    rename(mean="Est.", lower_mean="lower .95", upper_mean="upper .95") %>%
    select(-se) %>%
    mutate_at(c("mean", "lower_mean", "upper_mean"), as.numeric) %>%
    add_row(q = "20-80%", mean = mean.2080, lower_mean = mean.2080, upper_mean=mean.2080, .before = 3) %>%
    left_join(diff.table, by=c("q"))

}




#Other functions
#---------------

#Pairs quantile 20-80 with other quantiles for pairwise restricted mean survival time runs

my_extr_cats2 <- function(quant, cat, df, endpoints=c("I9_HYPTENS_AGE", "I9_HYPTENS")){
  
  tmp <- df %>% 
    rename_(CAT = cat) %>%
    filter(CAT == "20-80%") %>%
    mutate(CAT = 0)
  
  df %>%
    rename_(CAT = cat) %>%
    filter(CAT == quant) %>%
    mutate(CAT = 1) %>%
    bind_rows(tmp) %>%
    select(all_of(endpoints), CAT) %>%
    droplevels() %>%
    drop_na() 
  
} 

