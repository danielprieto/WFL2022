##Power Analysis for experimental design 
#Funding agency: Wings for Life
#Proposal:
##v1.0
#Author: Daniel Prieto dprieto(at)fcien.edu.uy
#Date:2022-10-26
##
##############################################
##Effect size determination of existing data##
##############################################
#>>>EdU proliferation assay<<<
##Data source: Preliminary data from our lab
#source DOI: unpublished
#Load data
data.prel <- read.csv("./c26fl_prolifer.csv")
data1 <- as.numeric(data.prel$C26Fl)#group A: Cx26fl/fl
data2 <- as.numeric(data.prel$Control)#group B: FoxJ1>tdT
data1 <- data2[1:23]#Discard NA values
##Check normality
shapiro.test(data1)#Shapiro-Wilk normality test
shapiro.test(data2)
#p>0.05 means we cannot reject the null hypothesis (that the data are normally distributed).
#
##Determine size of the effect
yA <- mean(data1)#mean of the measures for group A
yB <- mean(data2)#mean of the measures for group A
sd.pool <- sqrt(((sd(data1))^2+(sd(data2))^2)/2)#pooled standard deviation
d <- (yB-yA)/sd.pool #ES - Bausell & Li (2002)  DOI:10.1017/CBO9780511541933
d #Show effect size (Cohen's "d")
#
##################################################################
##A priori ower analysis of unknown data (effect size estimated)##
##################################################################
##One-way ANOVA
library(pwr)#Load library
#size of the effect (f)
#small=0.1; medium=0.25; large=0.4 #from Cohen's classification
#the number of groups in the proposed experiment (k)
#desired statistical power = 0.8
#returns "n", the number of samples required to achieve that power
pwr.anova.test(k=2, f=0.68, sig.level=0.05, power=0.8) 
#Plot Power=f(sample size)
##Plot data
library(ggplot2)
n <- c(seq(2,10,by=1),seq(12,20,by=2),seq(25,50,by=5))#define n intervals
p <- pwr.anova.test(k=3, f=0.68, sig.level=0.05, power=NULL, n=n) #calculate for interval
q<- qplot(p$n,p$power)+#create dot plot
  geom_point(size = 2, alpha=0.7)+ 
  theme(axis.text.x = element_text(hjust = 1, size = 16, color = "black"), 
        axis.text.y = element_text(hjust = 1, size = 16, color = "black"),
        axis.title=element_text(size=14,face="bold"),legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1))+
  labs(x = "Sample size", y = "Statistical power") +
  labs(title = "Power analysis", subtitle = "Balanced one-way analysis of variance power calculation") +
  theme(plot.title = element_text(size = 14, color = "black", face = "bold"), plot.subtitle=element_text(size = 12))+
  geom_ribbon(stat="smooth", color = "grey", size = 0.3, aes(outfit=fit<<-..y..), alpha=0.3, method = loess)
#
##Make model explicit
modl <- loess(p$n ~ p$power)#fit with Loess model of local polynomial adjustment
dPow <- 0.8#desired power
sSize <- predict(modl, dPow)#Sample size at value
sSize <- round(sSize, digits = 0)
#
#Put it into the plot
q <- q + geom_vline(xintercept = sSize, linetype="dotted", color = "black", size=1)# Add a vertical line segment
q <- q + geom_hline(yintercept = dPow, linetype="dotted", color = "grey", size=1)# Add horizontal line segment
q <- q + scale_y_continuous(breaks = sort(c(dPow, seq(0,1,0.25))), 
                                    labels=sort(c(dPow,seq(0,1,0.25))))
q + scale_x_continuous(breaks = sort(c(sSize, seq(0,100,25))), 
                          labels=sort(c(sSize,seq(0,100,25)))) +
theme(axis.text.x = element_text(angle = 45))#Show plot

########################################
##Bootstrap expansion of pilot studies##
########################################
##Data source:
#Load data here or use same as used for ES determination
#Bootstrapping is an approach to power analysis (Efron & Tibshirani, 1998)
#that assumes that the available data represents the population. 
#Random samples of a fixed size n with replacement from the existing dataset are drawn.
#The assumed population difference between the groups is simulated within each sample. 
#Power is the proportion of samples where a significant group difference was found (Walters, 2004).
#Code below based on code from Elizabeth Colantuoni:
#http://www.biostat.jhsph.edu/~ejohnson/regression/Sample%20Size%20Power%20Considerations.pdf
#Define 1000 pseudo replicates
power = function(sample1, sample2, reps=1000, size=10) {
  results <- sapply(1:reps, function(r) {
    resample1 <- sample(sample1, size=size, replace=TRUE)
    resample2 <- sample(sample2, size=size, replace=TRUE)
    test <- wilcox.test(resample1, resample2, paired=FALSE)
    test$p.value
  })
  sum(results<0.05)/reps
}
#Find power for sample sizes 3-100
pow <- c()
size <- c()
for (j in 3:100){
  pow[j] <- power(data1, data2, reps=1000, size=j)
  size[j] <- j
}
bPower <- as.data.frame(cbind(pow, size))
#Plot stuff
library(ggplot2)
Pplot <- ggplot(bPower, aes(y = pow, x = size))+
  geom_point(size = 2, alpha=0.7)+ 
  theme(axis.text.x = element_text(hjust = 1, size = 16, color = "black"), 
        axis.text.y = element_text(hjust = 1, size = 16, color = "black"),
        axis.title=element_text(size=14,face="bold"),legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size = 1))+
  labs(x = "Sample size", y = "Statistical power") +
  labs(title = "Power analysis", subtitle = "Bootstrap simulation") +
  theme(plot.title = element_text(size = 14, color = "black", face = "bold"), plot.subtitle=element_text(size = 12))+
  geom_ribbon(stat="smooth", color = "grey", size = 0.3, aes(outfit=fit<<-..y..), alpha=0.3, method = loess)
##Now make model explicit
mdl <- loess(size ~ pow)#fit with Loess model of local polynomial adjustment
#This is the sample size needed for the specified Power:
desPow <- 0.8#desired power
sampSize <- predict(mdl, desPow)#Sample size at value
sampSize <- round(sampSize, digits = 0)
##
#Put it into the plot
Pplot <- Pplot + geom_vline(xintercept = sampSize, linetype="dotted", color = "black", size=1)# Add a vertical line segment
Pplot <- Pplot + geom_hline(yintercept = desPow, linetype="dotted", color = "grey", size=1)# Add horizontal line segment
Pplot <- Pplot + scale_y_continuous(breaks = sort(c(desPow, seq(0,1,0.25))), 
                                    labels=sort(c(desPow,seq(0,1,0.25))))
Pplot+ scale_x_continuous(breaks = sort(c(sampSize, seq(0,100,25))), 
                          labels=sort(c(sampSize,seq(0,100,25))))+ 
  theme(axis.text.x = element_text(angle = 45))
#EOF