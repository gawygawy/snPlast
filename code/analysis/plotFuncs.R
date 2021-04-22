library(lemon) # to be able to use facet_rep_wrap
library(egg)
library(ggpubr)

cond.col.scheme <- c("#7570B3", "#D95F02", "#1B9E77")
cond.vecs <- c("GFP", "OFF", "ON")

plotGpPerf <- function(.gp.data){

  ggplot(.gp.data) + 
    geom_point(aes(x=day, y=success.p, colour=condition),  size=4, shape=19) + 
    geom_line(aes(x=day, y=success.p, colour=condition)) +
    geom_hline(yintercept = 0.8, linetype = "dashed") + 
    geom_errorbar(aes(x=day, ymin=success.p-sem, ymax=success.p+sem, colour=condition), width=.5) + 
    ylim(0, 1) +
    ylab("Successful trials/day") +
    xlab("Day") +
    #scale_y_continuous(labels = scales::number_format(accuracy = 0.1), breaks=seq(0, 1, 0.2)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks=seq(0, 1, 0.2)) +
    scale_x_continuous(labels = scales::number_format(accuracy = 1), breaks=seq(0, 12, 2)) +
    scale_color_manual(values=c(cond.col.scheme)) + 
    facet_rep_grid(~phase, scales="free_x", space="free_x") + 
    theme_classic(base_line_size = 1, base_size=20)#+ 
}

plotLogisticRegression <- function(re.df, fix.df){

ggplot() + 
  geom_line(aes(x=trial, y=predicted, colour=condition), size=1, data = as.tibble(fix.df) %>% mutate(condition=factor(condition, levels=c("GFP", "OFF", "ON")))) +
  geom_line(aes(x=trial, y=indiv.pred, group=sub, colour=condition), alpha=0.3, data=as.tibble(re.df) %>% mutate(condition=factor(condition, levels=c("GFP", "OFF", "ON")))) +
  facet_rep_grid(condition~phase,scales = "free_x", space = "free_x") + 
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks=seq(0, 1, 0.2)) +
    scale_x_continuous(labels = scales::number_format(accuracy = 1), breaks=seq(0, 12, 2)) +
  scale_color_manual(values = cond.col.scheme) + 
  theme_classic(base_size=22)
}

plotFixEfPreds <- function(fix.df){
  ggplot() + 
  geom_line(aes(x=trial, y=predicted, colour=condition), size=1, data = as.tibble(fix.df)) +
  facet_rep_grid(~phase,scales = "free_x", space = "free_x") + 
  scale_x_continuous(breaks=seq(1,12,1)) + 
  scale_y_continuous(labels=scales::percent) + 
  scale_color_manual(values = cond.col.scheme) + 
  my.theme
}

plotRandomEffects <- function(sub.re){
  ggplot(sub.re, aes(y=sub,x=condval, colour=condition)) +
    geom_point() + 
    facet_rep_grid(condition~term, scales="free") +
    geom_errorbarh(aes(xmin=condval -condsd,
                      xmax=condval +condsd), height=0) + 
    theme_classic() + 
    scale_color_manual(values=cond.col.scheme, guide=FALSE) +
    xlab("Coefficient value (deviation from average)") + 
    ylab("Mouse")
}

 plotHit80 <- function(hit.80.dat){
   ggplot(hit.80.dat, aes(x=condition, y=hit.80, fill=condition)) + 
     geom_violin(alpha = 0.25) + 
     geom_dotplot(binaxis="y", stackdir="center", dotsize = 0.6) +
     facet_rep_wrap(~phase) +
     ylab("days to attain 80% performance") + 
     scale_y_continuous(breaks=seq(2,12,2)) +
     scale_fill_manual(values=cond.col.scheme) +
     theme_classic(base_line_size = 1, base_size = 20)
 }

plotHit80Sig <- function(.hit80dat){

  compare_gp <- list(c("GFP", "OFF"), c("GFP", "ON"), c("ON", "OFF"))
  ggviolin(.hit80dat, x="condition", y="hit.80", fill="condition", add="dotplot", add.params=list(size=0.6, fill="condition"),
    palette=cond.col.scheme, facet.by="phase", trim=TRUE, alpha=0.4, 
    ylab = "days to reach 80%", xlab="group", order=c("GFP", "OFF", "ON")) + 
  scale_y_continuous(breaks=seq(2,12,2)) + 
  stat_compare_means(comparisons = compare_gp, label = "p.signif", size=5) + 
  theme_pubr(base_size=16)
 }

plotLearnCats <- function(.mhit80){

p <- ggbarplot(.mhit80 %>% group_by(cat, condition) %>% summarize(count=n()), 
    "condition", "count", 
    fill="cat", palette=c("white", "black", "brown", "grey"), alpha=0.5) + 
  theme_pubr(base_size=20)

ggpar(p,
   xlab = "group", ylab = "count", legend="right")
}

plotEstParams <- function(data, x, y, filename) {
  x <- enquo(x) # condition
  y <- enquo(y) # ACh or ACh.DA

  p <- ggplot(data, aes(x=!!x, y=!!y, fill=!!x)) + 
    geom_dotplot(binaxis="y", stackdir="center") + 
    geom_violin(alpha = 0.25) +
    scale_fill_manual(values=cond.col.scheme) +
    theme_classic() 

  if(!missing(filename)){
    ggsave(filename=filename, plot=p, device=cairo_pdf, dpi = 600)

  }

  return(p)
}

plotParamsPerf <- function(est.params, m.hit.80, y.param, learn.cat, filename){
  data <- left_join(est.params, m.hit.80)
  learn.cat <- enquo(learn.cat) # slow.learn or slow.rev? 
  y.param <- enquo(y.param) # ACh or ACh.DA?

  p <- ggplot(data) + 
    geom_dotplot(aes(x= factor(!!learn.cat), y= !!y.param, fill=condition), binaxis="y", stackdir="center") +
    facet_rep_wrap(~ condition) + 
    scale_fill_manual(values=cond.col.scheme) + 
    my.theme

  return(p)
}
#
# plotIndivFits <- function(est.day.performance, exp.dat, m.success.perc, which.con){

#   est.day.performance <- est.day.performance %>% rename(mod.sem = sem)

#   plt.dat <- left_join(est.day.performance, exp.dat %>% unnest(), 
#     by = c("sub", "condition", "phase", "day")) %>% 
#     left_join(., m.success.perc %>% unnest(), by = c("sub", "phase", "day")) %>% 
#     filter(condition == which.con)

#   con.col <- cond.col.scheme[cond.vecs == which.con]

#   ggplot(plt.dat) + 
#     geom_point(aes(x=trial, y=rew.found), colour=con.col, size=0.4) + 
#     geom_line(aes(x=day, y=predicted), colour=con.col) + 
#     geom_point(aes(x=day, y=success.p), colour=con.col, shape=1, size=2.5) + 
#     geom_ribbon(aes(x=day, y=predicted, ymin=predicted-mod.sem, ymax=predicted+mod.sem), fill = "grey70", alpha=0.5) +
#     geom_errorbar(aes(x=day, ymin=success.p-sem, ymax=success.p+sem), colour=con.col, width=.1) + 
#     facet_grid(sub~phase) +
#     ylim(0, 1) +
#     ylab("Proportion of successful trials/day") +
#     xlab("Day") +
#     scale_y_continuous(breaks=seq(0, 1)) +
#     scale_x_continuous(labels = scales::number_format(accuracy = 1), breaks=seq(2, 12, 2)) +
#     facet_rep_grid(sub~phase, scales="free_x", space="free_x") + 
#     theme_classic() +
#     theme(legend.position="none", 
#         panel.spacing.y = unit(0.3, "lines"), 
#         strip.background = element_rect(color=NA, fill=alpha("grey", 0.5)))

# }

# plotIndivSelect <- function(est.day.performance, exp.dat, m.success.perc, subs){

#   est.day.performance <- est.day.performance %>% rename(mod.sem = sem)

#   plt.dat <- left_join(est.day.performance, exp.dat %>% unnest(), 
#     by = c("sub", "condition", "phase", "day")) %>% 
#     left_join(., m.success.perc %>% unnest(), by = c("sub", "phase", "day")) %>% 
#     filter(sub %in% subs) %>% 
#     arrange(condition)

#   ggplot(plt.dat, aes(colour=condition)) + 
#     geom_point(aes(x=trial, y=rew.found), size=0.4) + 
#     geom_line(aes(x=day, y=predicted)) + 
#     geom_point(aes(x=day, y=success.p), shape=1, size=2.5) + 
#     #geom_ribbon(aes(x=day, y=predicted, ymin=predicted-mod.sem, ymax=predicted+mod.sem), fill = "grey70", alpha=0.5) +
#     geom_errorbar(aes(x=day, ymin=success.p-sem, ymax=success.p+sem), width=.1) + 
#     facet_grid(sub~phase) +
#     scale_color_manual(values=cond.col.scheme) + 
#     ylim(0, 1) +
#     ylab("Proportion of successful trials/day") +
#     xlab("Day") +
#     scale_y_continuous(breaks=seq(0, 1)) +
#     scale_x_continuous(labels = scales::number_format(accuracy = 1), breaks=seq(2, 12, 2)) +
#     facet_rep_grid(sub~phase, scales="free_x", space="free_x") + 
#     theme_classic() +
#     theme(legend.position="none", 
#         panel.spacing.y = unit(0.3, "lines"), 
#         strip.background = element_rect(color=NA, fill=alpha("grey", 0.5)))
# }


plotCV <- function(params0, params1, all, varX, varY){
  params0 %>% rename(ACh0 = ACh, DA0 = DA, ACh.DA0 = ACh.DA)
  params1 %>% rename(ACh1 = ACh, DA1 = DA, ACh.DA1 = ACh.DA)

  two.phase <- left_join(params0, params1, by = c("sub", "condition")) %>% 
    left_join(all, by = c("sub", "condition"))

  x <- enquo(varX)
  y <- enquo(varY)

  ggplot(two.phase, aes(x=!! x, y = !!y))

}

plotGpEstDayPerf <- function(est.day.performance, which.con){

group.avg <- est.day.performance %>% 
  group_by(condition, phase, day) %>% 
  summarize(success.p = mean(predicted), sem=sd(predicted)/sqrt(n()))


con.col <- cond.col.scheme[cond.vecs %in% which.con]
# predicted group performance 
  ggplot(group.avg, aes(x=day, y = success.p, colour = condition, group = condition)) + 
  geom_line() + 
  geom_errorbar(aes(x=day, ymin=success.p-sem, ymax=success.p+sem, colour=condition), width=.3) + 
  facet_rep_grid(~phase, scales="free_x", space="free_x") +
  scale_colour_manual(values=con.col) + 
  ylim(0, 1) +
  ylab("Successful trials/day") +
  xlab("Day") +
  #scale_y_continuous(labels = scales::number_format(accuracy = 0.1), breaks=seq(0, 1, 0.2)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks=seq(0, 1, 0.2)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 1), breaks=seq(2, 12, 2)) +
  theme_classic(base_size=20)
}

plotSubEstDayPerf <- function(est.day.performance){
# predicted day performance 
  ggplot(est.day.performance, aes(x=day, y = predicted, colour = condition, group = condition)) + 
    geom_line(aes(group = sub), alpha = 0.4) + 
    #stat_summary(fun.y = mean, geom = "line") + 
    scale_colour_manual(values = cond.col.scheme) +  
    ylim(0, 1) +
    ylab("successful trials/day") +
    xlab("day") +
    #scale_y_continuous(labels = scales::number_format(accuracy = 0.1), breaks=seq(0, 1, 0.2)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks=seq(0, 1, 0.2)) +
    scale_x_continuous(labels = scales::number_format(accuracy = 1), breaks=seq(1, 12, 1)) +
    facet_rep_grid(condition~phase, scales="free_x", space="free_x") + 
    my.theme
}

plotSubOnLearningRate <- function(sim.learning.rate, est.params, sub.cis){
  ggplot() + 
    geom_tile(aes(x=da, y=ach.da, fill=trial), data = sim.learning.rate, alpha = 0.6) + 

    geom_jitter(aes(x = da, y = ach.da), colour="black", height = 0, width = 0.00005, data = est.params, size=0.3) + 
    facet_rep_grid(condition~phase) + 
    scale_fill_gradient2(low = "blue", mid="yellow", high = "red", midpoint = 0.5) + 
    scale_colour_manual(values=cond.col.scheme) + 
    theme_classic() + 
    ylab("ach:da") +
    xlab("da")
}

plotSubOnHit80 <- function(.hit80dat, .phase, .params, .midpt){
  #midpt : for colour scale, 50 for 0, 70 for 1

  ggplot() + 
  geom_tile(aes(x=DA, y=ACh.DA, fill=hit.80), data = .hit80dat %>% filter(phase == .phase)) +  
  scale_fill_gradient2(low = "yellow", mid="red", high = "black", midpoint=.midpt, na.value = "black") + 
  geom_errorbar(aes(x=(da.low + da.hi)/2, ymin=achda.low, ymax=achda.hi), width=0, size=0.4, data=.params) + 
  geom_errorbarh(aes(y=(achda.low + achda.hi)/2, xmin=da.low, xmax=da.hi), height=0, size=0.4, data=.params) + 
  facet_wrap(~ condition, ncol=1) + 
  theme_classic(base_size=16)
}

plotModHit80 <- function(.hit80dat, .phase,.midpt){
  #midpt : for colour scale, 50 for 0, 70 for 1

  ggplot() + 
  geom_tile(aes(x=DA, y=ACh.DA, fill=hit.80), data = .hit80dat %>% filter(phase == .phase)) +  
  scale_fill_gradient2(low = "yellow", mid="red", high = "black", midpoint=.midpt, na.value = "black") + 
  theme_classic(base_size=24)
}
  
# ggplot(est.hit80, aes(x=condition, y=hit.80, fill=condition)) + 
#     geom_dotplot(binaxis="y", stackdir="down", position=position_nudge(x=-0.025, y=0), dotsize=0.5, binwidth=0.5) +
#     geom_flat_violin(scale="count", trim=FALSE) + 
#     facet_wrap(~phase)


plotGpParams <- function(.gp.params, y.var, facet.var=NULL){

  my_comparisons <- list( c("GFP", "OFF"), c("GFP", "ON"), c("ON", "OFF"))

  ggboxplot(.gp.params, x="sub", y=y.var, facet.by=facet.var)+ 
    stat_compare_means(comparisons = my_comparisons, label = "p.signif", size=3)
}

plotParamsByCat <- function(data, .x1, .x2, .cat){
  x1 <- enquo(.x1)
  x2 <- enquo(.x2)
  cat <- enquo(.cat)

  dat <- data %>% 
    select(sub, !!x1, !!x2, !!cat, condition) %>% 
    gather(phase, new_x, -sub, -!!cat, -condition)

  p <- ggplot(dat, aes(x=phase, y=new_x, group=sub)) + 
    geom_line(size=0.3) +
    geom_point() + 
    facet_wrap(vars(!!cat)) + 
    theme_classic()

  return(p)
}

daVals <- c(0.75, 1.75, 2.75, 0.75, 1.75, 2.75)
xVals <- c(0.00075, 0.00175, 0.00275, 0.00315, 0.00415, 0.00515)

plotParamsByPhase <- function(data){
  
  ggplot(dat3, aes(x=DA, y=ACh, group=sub)) + 
    geom_rect(xmin=0.00075, xmax=0.00275, ymin = -Inf, ymax=Inf, alpha=0.2, fill="grey") + 
    geom_vline(xintercept = c(0.00175, 0.00415), linetype="dashed", size=0.3) +
    geom_vline(xintercept = c(0.00075, 0.00275, 0.00315, 0.00515)) +
    geom_line() + 
    geom_point() + 
    annotate(geom="text", x=0.00175, y=0.0015, label="Learning", size=5) + 
    annotate(geom="text", x=0.00415, y=0.0015, label="Reversal", size=5) + 
    facet_grid(slow.rev ~ condition) +
    scale_color_brewer(palette="Dark2") + 
    xlim(0.00075, 0.00555) + 
    scale_x_continuous(breaks=xVals, labels=daVals) + 
    theme(axis.text.x = element_text(angle=45, size=rel(1.5)), 
      axis.text.y = element_text(size = rel(1.5)),
      panel.background = element_blank(), 
      panel.border = element_blank())
  
}

plotRecoveredParams <- function(.recovered.dat){
  ggscatter(.recovered.dat, x = "groundTruth", y = "recovered", add = "reg.line", facet.by="param", scales="free") +
  stat_cor() + 
  geom_abline(xintercept=0, linetype = "dashed")
}

subLabeller <- function(facet.label){
  return(lapply(facet.label,function(x) str_split(x, "_", simplify=TRUE)[, 2]))
}

plotMainSubFits <- function(.data){
  main.plot <- 
  ggplot(.data, aes(group=orderSub)) + 
    geom_point(aes(x=day, y=success.p), size=0.8) + 
    geom_line(aes(x=day, y=predicted))+
    geom_errorbar(aes(x=day, ymin=predicted-sem, ymax=predicted+sem), width=.5) +  
    facet_wrap(~orderSub, ncol=4, labeller=subLabeller) + 
    theme_bw(base_size=9) + 
    ylab("successful trials/day") + 
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks=seq(0, 1, 0.2)) +
      scale_x_continuous(labels = scales::number_format(accuracy = 1), breaks=seq(2, 20, 2)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.position = "bottom",
          legend.title = element_blank())
    
    return(main.plot)
}

get_inset <- function(df){
  label <- data.frame(
  ACh = c(0.00155, 0.00155), 
  DA = c(0.00175, 0.00415), 
  label = c("init.", "rev.")
  )
  
p <- ggplot(data=df, aes(x=DA, y=ACh)) + 
      geom_rect(xmin=0.00075, xmax=0.00275, ymin = -Inf, ymax=Inf, alpha=0.2, fill="grey") + 
      geom_vline(xintercept = c(0.00175, 0.00415), linetype="dashed", size=0.3) +
      geom_vline(xintercept = c(0.00075, 0.00275, 0.00315, 0.00515)) +
      geom_text(data=label, aes(label = label), size=1.75) +
      geom_line()+ 
      geom_point(size=0.35)+ 
      ylim(0, 0.0016) +
      xlim(0.00075, 0.00555) + 
      ylab("ACh") + 
      scale_color_brewer(palette="Dark2") + 
      scale_x_continuous(breaks=xVals, labels=daVals) + 
      theme_bw(base_size=6) + 
          guides(fill=FALSE) +
        theme(panel.background = element_rect(fill="white"),  ## white plot background 
              axis.title.y = element_text(size=rel(0.8), vjust=0.1),
              axis.title.x = element_text(size=rel(0.8)),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.x=element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              legend.position = "none", 
              plot.background = element_blank())
  return(p)
}

## This function allows us to specify which facet to annotate
annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
{
  layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob, 
                                          xmin = xmin, xmax = xmax, 
                                          ymin = ymin, ymax = ymax))
}

plotInsetSubFits <- function(.data){
  
 insets <- .data %>% 
  split(f = .$orderSub) %>% 
  map(~annotation_custom2(
    grob=ggplotGrob(get_inset(.)), 
    data = data.frame(orderSub=unique(.$orderSub)),
    #ymin=-0.095, ymax=0.45, xmin = 11, xmax=22) # OFF settings 
    ymin=-0.125, ymax=0.45, xmin = 11, xmax=22) # ON settings 
    #ymin=-0.1, ymax=0.5, xmin = 11, xmax=22) # GFP settings 
  ) 

  return(insets)
}

plotMainSubFits2 <- function(.data){
 main.plot2 <- 
  ggplot(.data, aes(group=sub)) + 
    geom_point(aes(x=day, y=success.p, colour=condition)) + 
    geom_line(aes(x=day, y=predicted, colour=condition))+
    geom_errorbar(aes(x=day, ymin=predicted-sem, ymax=predicted+sem, colour=condition), width=.5) +  
    facet_rep_grid(sub ~ phase, scales="free_x", space="free_x") + 
    theme_bw(base_size=14) + 
    scale_color_manual(values=cond.col.scheme) + 
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks=seq(0, 1, 0.2)) +
      scale_x_continuous(labels = scales::number_format(accuracy = 1), breaks=seq(2, 20, 2)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          legend.position = "bottom",
          legend.title = element_blank()) 

  return(main.plot2)
}

get_inset2 <- function(df){
  p <- ggplot(data=df, aes(x=DA, y=ACh)) + 
      geom_vline(xintercept = 0.00175, linetype="dashed", size=0.3) +
      geom_vline(xintercept = c(0.00075, 0.00275)) +
      geom_point(size=0.2)+ 
      geom_errorbar(aes(ymin=ach.low, ymax=ach.hi), width=0, size=0.1) + 
      geom_errorbarh(aes(xmin=da.low, xmax=da.hi), height=0, size=0.1) + 
      ylab("ACh(1e-03)") +
      xlab("DA(1e-03)") +
      scale_color_brewer(palette="Dark2") + 
      scale_x_continuous(breaks=xVals, labels=daVals) + 
      scale_y_continuous(labels = scales::number_format(scale=1000, accuracy=0.1), limits=c(0, 0.0016)) + 
      theme_bw(base_size=9) + 
          guides(fill=FALSE) +  
        theme(panel.background = element_rect(fill="white"),  ## white plot background 
              axis.title.y = element_text(size=rel(0.8)),
              axis.title.x = element_text(size=rel(0.8)),
              axis.text.x = element_text(size=rel(0.8)),
              axis.text.y = element_text(size=rel(0.8)),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              legend.position = "none", 
              plot.background = element_blank())
  return(p)
}

plotInsetSubFits2 <- function(.data){
  
 insets <- .data %>% 
  split(f = .$sub) %>% 
  map(~annotation_custom2(
    grob=ggplotGrob(get_inset2(.)), 
    data = data.frame(sub=unique(.$sub)),
    ymin=-0.08, ymax=0.6, xmin = 7, xmax=12.5)
  ) 

  return(insets)
}

plotRecoveredParams <- function(.data){
  ggscatter(.data, x = "groundTruth", y = "recovered", add = "reg.line", facet.by="param", scales="free") +
    stat_cor() + 
    geom_abline(linetype = "dashed")
}
