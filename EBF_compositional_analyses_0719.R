library(ggplot2)
library(ggsignif)
library(phyloseq)
library(PMA)
library(GGally)
install.packages("glmnet")
source("/Users/bpb/Documents/GitHub/GutMicrobes_analyses/CoDA_functions_053119.R") #available @ https://github.com/itsmisterbrown/GutMicrobes_analyses


#load data, pseudocount of +1 already added for logratio analysis
feed <- phyloseq::import_biom("/Users/bpb/Documents/GitHub/GutMicrobes_analyses/Filtered_Feeding_053119_ASV_table_w_tax_md.biom")

#update colnames
colnames(tax_table(feed)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#change theme
theme_set(theme_bw())

#set aliases
phy.obj <- feed
asv.tab <- t(otu_table(phy.obj))
sampledf <- sample_data(phy.obj)
tax.tab <- tax_table(phy.obj)

#convert sample data to numeric
sampledf$TNF <- as.numeric(paste0(sampledf$TNF))
sampledf$HLADR <- as.numeric(paste0(sampledf$HLADR))
sampledf$CCR5 <- as.numeric(paste0(sampledf$CCR5))
sampledf$CD25hiCD39plus <- as.numeric(paste0(sampledf$CD25hiCD39plus))
sampledf$CCR5plusCD25plus <- as.numeric(paste0(sampledf$CCR5plusCD25plus))
sampledf$CD25hiCD39 <- as.numeric(paste0(sampledf$CD25hiCD39))
sampledf$CD25plus <- as.numeric(paste0(sampledf$CD25plus))
sampledf$HLADRplusCD25plus <- as.numeric(paste0(sampledf$HLADRplusCD25plus))

#perform unweighted CLR transform and sparse principal component analysis via L1 PMD
#CLR transform
asv.tab.clr <- clr(asv.tab)
#SPC via L1 PMD
spc.orth <- PMA::SPC(unclass(asv.tab.clr), sumabsv=2.5, K=90, orth=TRUE, center = FALSE, trace = F) #values used in manuscript analysis

#update colnames and generate unweighted ILR basis from L1 PMD output
rownames(spc.orth$v) <- colnames(asv.tab.clr)
#create ILR part weights vector
ilr.weights <- rep(1, ncol(asv.tab.clr))
names(ilr.weights) <- colnames(asv.tab.clr)
#build basis for ILR transform (balances)
ilr.basis.df <- balance.basis(spc.orth$v, p=ilr.weights)
colnames(ilr.basis.df) <- paste("SPB", 1:ncol(ilr.basis.df), sep = ".")

#Perform the ILR using the basis identified via L1 PMD (set of balances)
#perform the ILR
ilr.asv.tab <- ilr(asv.tab, V = ilr.basis.df, p = ilr.weights)
ilr.asv.tab[is.nan(ilr.asv.tab)] <- 0 #set any balances without numerator and denominator parts to 0 to work with glmnet


#Perform penalized logistic regression to identify balances associated with feeding practice
#penalized logistic regression, binomial
l1mod <- glmnet(ilr.asv.tab, sampledf$EBF, alpha=1, family = "binomial")

ebf.bals <- as.matrix(coefficients(l1mod, s=0.1))
ebf.bals <- rownames(ebf.bals)[which(ebf.bals != 0)]
(ebf.bals <- ebf.bals[2:length(ebf.bals)]) #the first entry is just the intercept, which we don't want

gpair <- GGally::ggpairs(cbind.data.frame(ilr.asv.tab[,ebf.bals], sampledf$EBF), 
                 upper = list(combo = "box", continuous = "cor"), lower = list(continuous = "smooth")) #SPB3 looks good (referred to as Balance 1 in manuscript)

#examine SPB3
eb <- extract.balances(W = spc.orth$v, taxonomy = tax.tab)
spb3.tax <- data.frame(eb[3])

#examine asv tab when subset to taxa in Balance 1
asv.tab.spb3 <- asv.tab[,rownames(data.frame(spb3.tax))]

#subset Balance 1 to ASVs present in > 33% of samples
#determine number of samples that each ASV is present in (prevalence)
asv.spb3.prev <- apply(X = data.frame(asv.tab.spb3), 2, FUN = function(x) length(x[x>2])) #>2 since we added a pseudocount of 1

#subset to ASVs with prevalence > 33%
asv.spb3.hiprev <- asv.spb3.prev[which(asv.spb3.prev> 42)] #42 sample minimum corresponds to ~33% prevalence cutoff

#examine taxonomy of SPB3 (Balance 1) with ASVs with high prevalence
spb3.tax[names(asv.spb3.hiprev),] #Balance 1

#create vector of values from subset of SPB3 (Balance 1)
balance1 <- create.balance(df = asv.tab, num.tax = c("ASV2", "ASV6", "ASV13"), den.tax = c("ASV1", "ASV3", "ASV7", "ASV8", "ASV10"), weighted = F)
colnames(balance1) <- "Balance1"

#merge sample data and values from reduced balance
com.df <- cbind.data.frame(sampledf, balance1)

#adjust factor levels for plotting
com.df$TimePoint <- factor(com.df$TimePoint, levels = c("Birth", "Week 6", "Week 14"))

#wilcox tests since sample sizes vary
wilcox.test(formula=Balance1~as.factor(EBF), data=com.df, subset = c(TimePoint %in% c("Week 14")))

#plot balances longitudinally
ggplot(as.data.frame(com.df), aes(x=TimePoint, y=Balance1, color=EBF)) +
  geom_boxplot(size=2) +
  geom_hline(yintercept = 0, lty=2, size=1.5) +
  theme(axis.text.y = element_text(size = 25, colour = "black"),
        axis.title.y = element_text(size = 15, colour = "black"),
        axis.text.x = element_text(size = 25, colour = "black"),
        axis.title.x = element_text(size = 0, colour = "black"),
        legend.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 0),
        legend.position = c(0.1, 0.9),
        strip.background =element_rect(fill="white"), 
        strip.text = element_text(size = 35, colour = "black")) +
  scale_color_manual(values = rev(fcols)) +
  geom_signif(y_position=c(8.25), xmin=c(2.8), xmax=c(3.2),
              annotation=c("*"), textsize = 10, size=2, color="black", vjust = 0.5) +
  labs(y=expression(paste(log(italic(over("B.breve"%*%"B.gallicum"%*%"B.bifidum",
                                          "E.coli"%*%"B.dorei"%*%"V.dispar"%*%"B.vulgatus"%*%"R.gnavus"))))))

#examine interactions between Balance 1 and flow cytometry data
flownames <- c("CCR5", "CCR5plusCD25plus", "CD25plus", "HLADR", "HLADRplusCD25plus", "Balance1")
gpair2 <- GGally::ggpairs(com.df, columns = flownames, mapping = aes(color=EBF),
                         upper = list(combo = "box", continuous = "cor"), lower = list(continuous = "smooth")) + scale_y_log10()

#specific flow models with potentially significant interactions
lmt1 <- glm(formula = log10(CD25plus)~Balance1, data = com.df, family = gaussian(link = identity))
lmt2 <- glm(formula = log10(CCR5plusCD25plus)~Balance1, data = com.df, family = gaussian(link = identity))
lmt3 <- glm(formula = log10(HLADRplusCD25plus)~Balance1, data = com.df, family = gaussian(link = identity))

lmt1
summary(lmt3)

#plot model fit
#create labels for strip text
com.df$hladrlab <- c(rep("HLA-DR+CD25+", nrow(com.df)))
com.df$CCR5lab <- c(rep("CCR5+CD25+", nrow(com.df)))
com.df$cd25pluslab <- c(rep("CD25+", nrow(com.df)))

#generate plots
plot1 <- ggplot(com.df, aes(x=Balance1, y=CD25plus)) +
  geom_point(size=5, aes(color=EBF)) +
  geom_point(color = "grey90", size = 2) +
  scale_color_manual(values = rev(fcols)) +
  theme(axis.text.y = element_text(size = 25, colour = "black"),
        axis.title.y = element_text(size = 25, colour = "black"),
        axis.text.x = element_text(size = 25, colour = "black"),
        axis.title.x = element_text(size = 25, colour = "black"),
        legend.text = element_text(size = 20, colour = "black"),
        legend.title = element_text(size = 0),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
        legend.position=c(0.15, 0.15),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 25)) +
  scale_y_log10(limits=c(1,100)) + annotation_logticks(sides = "l") +
  lims(x=c(-9,7.5)) +
  facet_wrap(~cd25pluslab) +
  geom_smooth(method = "glm", data = com.df, formula = y~x, method.args = list(family = gaussian(link = identity)), se = TRUE) + 
  labs(y="Expression", x = "Balance 1")


plot2 <- ggplot(com.df, aes(x=Balance1, y=CCR5plusCD25plus)) +
  geom_point(size=5, aes(color=EBF)) +
  geom_point(color = "grey90", size = 2) +
  scale_color_manual(values = rev(fcols)) +
  theme(axis.text.y = element_text(size = 25, colour = "black"),
        axis.title.y = element_text(size = 25, colour = "black"),
        axis.text.x = element_text(size = 25, colour = "black"),
        axis.title.x = element_text(size = 25, colour = "black"),
        legend.text = element_text(size = 20, colour = "black"),
        legend.title = element_text(size = 0),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
        legend.position="none",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 25)) +
  scale_y_log10(limits = c(1, 100)) + annotation_logticks(sides = "l") +
  lims(x=c(-9,7.5)) +
  facet_wrap(~CCR5lab) +
  geom_smooth(method = "glm", data = com.df, formula = y~x, method.args = list(family = gaussian(link = identity)), se = TRUE) + 
  labs(x="Balance 1", y="Expression")


plot3 <- ggplot(com.df, aes(x=Balance1, y=HLADRplusCD25plus)) +
  geom_point(size=5, aes(color=EBF)) +
  geom_point(color = "grey90", size = 2) +
  scale_color_manual(values = rev(fcols)) +
  theme(axis.text.y = element_text(size = 25, colour = "black"),
        axis.title.y = element_text(size = 25, colour = "black"),
        axis.text.x = element_text(size = 25, colour = "black"),
        axis.title.x = element_text(size = 25, colour = "black"),
        legend.text = element_text(size = 20, colour = "black"),
        legend.title = element_text(size = 0),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
        legend.position="none",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 25)) +
  scale_y_log10(limits=c(0.1,100)) + annotation_logticks(sides = "l") +
  lims(x=c(-9,7.5)) +
  facet_wrap(~hladrlab) +
  geom_smooth(method = "glm", data = com.df, formula = y~x, method.args = list(family = gaussian(link = identity)), se = TRUE) + 
  labs(x="Balance 1", y="Expression")

#put it all together
flow.plot <- cowplot::plot_grid(plot1, plot2, plot3, align = "h", nrow=2, labels = "AUTO", 
                                label_size = 20)
flow.plot

View(com.df)

#summary stats by group
com.df %>% 
  group_by(EBF, TimePoint) %>% 
  summarise(median = median(Balance1), mean = mean(Balance1))
