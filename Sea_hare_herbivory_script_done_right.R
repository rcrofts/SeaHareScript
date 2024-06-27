Algal.herbivory.trials <- read.csv("C:/Users/Rory/Documents/Uni/Masters/Research/Stats/Algal.herbivory.trials.csv")
attach(Algal.herbivory.trials)

library(tidyverse)
library(dplyr)
library(car)
library(MASS)
library(rcompanion)
library(ggpubr)
library(broom)

##replace negative values for algae eaten with zeros
algae_no_negs<-ifelse(algae.eaten.ml<0, 0, algae.eaten.ml)
hist(sqrt(algae_no_negs))
length(algae_no_negs)

Algal.herbivory.trials <- cbind(Algal.herbivory.trials, algae_no_negs)

##store data under this name

seahare_data <- Algal.herbivory.trials

##check a random sample

set.seed(1234)
dplyr::sample_n(seahare_data, 10)

## check the structure

str(seahare_data)

## convert treatment into factor

seahare_data$treatment <- as.factor(seahare_data$treatment)


## print head of the file

head(seahare_data)

## show group levels

levels(seahare_data$treatment)

## reorder so control is first

seahare_data$treatment <- ordered(seahare_data$treatment,
                                  levels = c("Control", "Aplysia", "Dolabella.a", "Dolabella.u"))

## compute summary statistics by group

group_by(seahare_data, treatment) %>%
  summarise(
    count = n(),
    mean = mean(algae_no_negs, na.rm = TRUE),
    sd = sd(algae_no_negs, na.rm = TRUE),
    median = median(algae_no_negs, na.rm = TRUE),
    IQR = IQR(algae_no_negs, na.rm = TRUE)
  )

## visualize data in a box plot

ggboxplot(seahare_data, x = "treatment", y = "algae_no_negs", 
          color = "treatment", palette = c("#00AFBB", "#E7B800", "#FC4E07", "#adc178"),
          order = c("Control", "Aplysia", "Dolabella.a", "Dolabella.u"),
          ylab = "Algae eaten", xlab = "Treatment")

## visualize data in a line graph

ggline(seahare_data, x = "treatment", y = "algae_no_negs", 
       add = c("mean_se", "jitter"), 
       order = c("Control", "Aplysia", "Dolabella.a", "Dolabella.u"),
       ylab = "Algae eaten", xlab = "Treatment")

## compute Kruskal-Wallis test

kruskal.test(algae_no_negs ~ treatment, data = seahare_data)

## p < 0.05 so significant differences between groups

## multiple pairwise comparison between groups

pairwise.wilcox.test(seahare_data$algae_no_negs, seahare_data$treatment,
                     p.adjust.method = "BH")






## analysis of algae eaten vs sea hare size
plot(algae_no_negs, collection.size.ml)
size.vs.herb.model <- lm(collection.size.ml ~ algae_no_negs)
abline(size.vs.herb.model)
summary(size.vs.herb.model)

##select only Dolabella auricularia to investigate relationship
algae.eaten.dola <- subset(Algal.herbivory.trials, treatment == "Dolabella.a")

##select only Dolabella ? to investigate relationship
algae.eaten.dola.u <- subset(Algal.herbivory.trials, treatment == "Dolabella.u")

##select only Aplysias due to the range of sizes
proportional.eaten.aplysia <- subset(Algal.herbivory.trials, treatment == "Aplysia")

## algae eaten vs size by species

plot(algae.eaten.dola$algae_no_negs, algae.eaten.dola$collection.size.ml, xlab = "Algae eaten per 24h (ml)", ylab = "Collection size (ml)")
size.vs.herb.model.dola <- lm(algae.eaten.dola$collection.size.ml ~ algae.eaten.dola$algae_no_negs)
abline(size.vs.herb.model.dola)
summary(size.vs.herb.model.dola)

plot(algae.eaten.dola.u$algae_no_negs, algae.eaten.dola.u$collection.size.ml, xlab = "Algae eaten per 24h (ml)", ylab = "Collection size (ml)")
size.vs.herb.model.dola.u <- lm(algae.eaten.dola.u$collection.size.ml ~ algae.eaten.dola.u$algae_no_negs)
abline(size.vs.herb.model.dola.u)
summary(size.vs.herb.model.dola.u)

plot(proportional.eaten.aplysia$algae_no_negs, proportional.eaten.aplysia$collection.size.ml, xlab = "Algae eaten per 24h (ml)", ylab = "Collection size (ml)")
size.vs.herb.model.aplys <- lm(proportional.eaten.aplysia$collection.size.ml ~ proportional.eaten.aplysia$algae_no_negs)
abline(size.vs.herb.model.aplys)
summary(size.vs.herb.model.aplys)

## proportional algae eaten compared to sea hare size (Aplysia only)
proportional.algae.eaten <- proportional.eaten.aplysia$algae_no_negs/proportional.eaten.aplysia$collection.size.ml
proportional.eaten.aplysia <- cbind(proportional.eaten.aplysia, proportional.algae.eaten)


## looking at relationship between proportional algae eaten and size
proportional.eaten.aplysia <- na.omit(proportional.eaten.aplysia)
plot((proportional.eaten.aplysia$proportional.algae.eaten) ~ proportional.eaten.aplysia$collection.size.ml, xlab = "Collection size (ml)", ylab = "Algae eaten as a proportion of body size (ml/ml 24h)")
summary(proportional.eaten.aplysia$proportional.algae.eaten)

boxplot(proportional.eaten.aplysia$proportional.algae.eaten)


## fit an exponential decay model to this curve

nlm_prop_alg <- nls(proportional.algae.eaten ~ SSasymp(collection.size.ml, Asym, R0, lrc), data = proportional.eaten.aplysia)

nlm_prop_alg
summary(nlm_prop_alg)

theme_set(theme_bw(12))


qplot(collection.size.ml,proportional.algae.eaten, data = augment(nlm_prop_alg), xlab = "Collection size (ml)", ylab = "Algae eaten as a proportion of body size (ml/ml 24h)") + geom_line(aes(y = .fitted))

ggplot(proportional.eaten.aplysia, aes(x = proportional.eaten.aplysia$collection.size.ml, y = proportional.eaten.aplysia$proportional.algae.eaten)) +
  geom_point() +
  geom_smooth(method = "nls", formula = y ~ 0.45792 + (23.03373 - 0.45792)*exp(-exp(-0.99304)*x))


xvals <- seq(0,600, by = 5)
yvals <- 0.45792 + (23.03373 - 0.45792)*exp(-exp(-0.99304)*xvals) 
plot(xvals, yvals, type = "l")



nlm_prop_alg2 <- predict(nlm_prop_alg,newdata = data.frame(collection.size.ml = seq(0,1000, by = 10)))

nlm_prop_alg2

k = coef(nlm_prop_alg)
k_s_SE = summary(nlm_prop_alg)$coefficients[1,2]
k_s1_LCI = confint(nlm_prop_alg)[1]
k_s1_UCI = confint(nlm_prop_alg)[2]
R2_s = cor(na.omit(proportional.eaten.aplysia$proportional.algae.eaten), predict(nlm_prop_alg))^2
AIC = AIC(nlm_prop_alg)

k

-exp(-0.99304)

## investigating algal composition data
algal.distribution.tabulated <- read.csv("C:/Users/Rory/Documents/Uni/Masters/Research/Stats/algal.distribution.tabulated.csv")

## chi-square test to determine if values are significantly different from expected based on counts
chisq.test(algal.distribution.tabulated$No.Quadrats.present)
## they are

## investigate relationship of algae displacement values vs dry weight
algae.dry.weights <- read.csv("~/Documents/Uni/Masters/Research/Stats/algae.dry.weights.csv")
plot(algae.dry.weights$Water.displaced..ml., algae.dry.weights$Dry.weight..g., xlab = "Water displaced (ml)", ylab = "Dry weight (g)")
alg.dry.weight.model <- lm(algae.dry.weights$Dry.weight..g. ~ algae.dry.weights$Water.displaced..ml.)
abline(alg.dry.weight.model)
summary(alg.dry.weight.model)
coef(alg.dry.weight.model)

## is there any difference between groups of algae bundles for this relationship
boxplot(algae.dry.weights$Water.displaced..ml. ~ algae.dry.weights$Group)

## looks very similar
dry.weight.rem.1 <- filter(algae.dry.weights, Group == "Rem.1")
dry.weight.rem.2 <- filter(algae.dry.weights, Group == "Rem.2")
dry.weight.pure <- filter(algae.dry.weights, Group == "Pure")

## repeat analysis for each group to see which is most accurate

plot(dry.weight.rem.1$Water.displaced..ml., dry.weight.rem.1$Dry.weight..g.)
dry.weight.rem1.model <- lm(dry.weight.rem.1$Dry.weight..g. ~ dry.weight.rem.1$Water.displaced..ml.)
abline(dry.weight.rem1.model)
summary(dry.weight.rem1.model)
coef(dry.weight.rem1.model)

plot(dry.weight.rem.2$Water.displaced..ml., dry.weight.rem.2$Dry.weight..g.)
dry.weight.rem2.model <- lm(dry.weight.rem.2$Dry.weight..g. ~ dry.weight.rem.2$Water.displaced..ml.)
abline(dry.weight.rem2.model)
summary(dry.weight.rem2.model)
coef(dry.weight.rem2.model)

plot(dry.weight.pure$Water.displaced..ml., dry.weight.pure$Dry.weight..g., xlab = "Water displaced (ml)", ylab = "Dry weight (g)")
dry.weight.pure.model <- lm(dry.weight.pure$Dry.weight..g. ~ dry.weight.pure$Water.displaced..ml.)
abline(dry.weight.pure.model)
summary(dry.weight.pure.model)
coef(dry.weight.pure.model)

## Pure has the highest R squared, may be the best fit.

##finding mean algae eaten per 24h for Aplysia
mean(proportional.eaten.aplysia$algae.eaten.ml)
summary(proportional.eaten.aplysia$algae.eaten.ml)

## translating to grams (using overall formula)
(0.61281*111.6667) + 3.95630

##translating to grams using Pure formula
(0.69634*111.667) - 0.86457

76.89363/24

## each Aplysia eats on average 3.203901 g of Laurencia per hour

