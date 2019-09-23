library(magrittr)
#devtools::install_github("caleblareau/BuenColors")
library(BuenColors)
library(xtable)

#load truth
load("frontdoor-truth.rda")

#load results
correct <- readRDS("data/out_correct.rds")
correct2 <- readRDS("data/out_correct2.rds")
miss1 <- readRDS("data/out_miss1.rds")
miss2 <- readRDS("data/out_miss2.rds")

#save results in R space
mle.correct <- round(colMeans(correct[[1]]),3)
mle.correct2 <- round(colMeans(correct2[[1]]),3)
mle.miss1 <- round(colMeans(miss1[[1]]),3)
mle.miss2 <- round(colMeans(miss2[[1]]),3)

sp.1.correct <- round(colMeans(correct[[2]]),3)
sp.1.correct2 <- round(colMeans(correct2[[2]]),3)
sp.1.miss1 <- round(colMeans(miss1[[2]]),3)
sp.1.miss2 <- round(colMeans(miss2[[2]]),3)

sp.2.correct <- round(colMeans(correct[[3]]),3)
sp.2.correct2 <- round(colMeans(correct2[[3]]),3)
sp.2.miss1 <- round(colMeans(miss1[[3]]),3)
sp.2.miss2 <- round(colMeans(miss2[[3]]),3)

sp.correct <- round(colMeans(correct[[4]]),3)
sp.correct2 <- round(colMeans(correct2[[4]]),3)
sp.miss1 <- round(colMeans(miss1[[4]]),3)
sp.miss2 <- round(colMeans(miss2[[4]]),3)

#TABLE#
results <- rbind(mle.correct,sp.1.correct,sp.2.correct,sp.correct,
                 mle.correct2,sp.1.correct2,sp.2.correct2,sp.correct2,
                 mle.miss1,sp.1.miss1,sp.2.miss1,sp.miss1,
                 mle.miss2,sp.1.miss2,sp.2.miss2,sp.miss2)
colnames(results) <- c("PSI","PIIE","VAR","BIAS","COVERAGE")

xtable(results)

#PLOT#


results.figure <- data.frame(cbind(rbind( cbind(1,c(correct[[1]][,2],correct[[2]][,2],correct[[3]][,2],correct[[4]][,2])),
                                          cbind(2,c(correct2[[1]][,2],correct2[[2]][,2],correct2[[3]][,2],correct2[[4]][,2])),
                                          cbind(3,c(miss1[[1]][,2],miss1[[2]][,2],miss1[[3]][,2],miss1[[4]][,2])),
                                          cbind(4,c(miss2[[1]][,2],miss2[[2]][,2],miss2[[3]][,2],miss2[[4]][,2]))),rep(sort(rep(1:4,10000)),4)))

colnames(results.figure) <- c("Model","PIIE","Estimator")
       
plot <- ggplot(results.figure, aes(x=factor(Model), y=PIIE, fill=factor(Estimator))) 
p2 <- plot + geom_boxplot() + 
  scale_x_discrete("\n Model specification",labels=c("1" = "(a)", "2" = "(b)", "3" = "(c)", "4" = "(d)")) + 
  ylab("population intervention indirect effect (PIIE) \n") +
  geom_hline(yintercept=truth.est[3],color="black") + 
  scale_fill_manual(values = c("white","gray87","gray51","gray29"),name="Estimator Type",labels=c("MLE","SP 1","SP 2","SP DR")) +
  pretty_plot() +
  theme(legend.key.height = unit(1, "cm"),legend.key.width = unit(1, "cm"),legend.text=element_text(size=20),axis.title=element_text(size=24),text = element_text(size=20),axis.text.x = element_text(size=20, angle=0),axis.text.y = element_text(size=20, angle=0)) 

ggsave(p2, file = "sim-results-notitle.png", dpi = 300, width = 15, height = 10)
