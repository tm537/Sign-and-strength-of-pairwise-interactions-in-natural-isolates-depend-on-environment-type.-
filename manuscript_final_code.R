library(tidyverse)
library(pbkrtest)
library(lmerTest)
library(DHARMa)
library(glmmTMB)
library(Matrix)
library(lme4)
library(emmeans)
library(gllvm)
library(glmmTMB)
library(car)
library(effects)
library(MuMIn)
library(rstatix)
library(ggpubr)
library(broom)
library(igraph)
library(viridis)
library(ggraph)
library(car)
library(ggplot2)
library(wesanderson)
library(mosaic)
library(cowplot)
library(lme4)
library(nlme)
library(multcompView)
library(multcomp)
library(ggrepel)
library(grid)
library(dglm)
library(vegan)
library(tidyverse)
library(emmeans)
library(magrittr)
library(broom.mixed)
library(viridis)
library(dplyr)
library(lmerTest)
library(caret)
library(igraph)
library(ggraph)
library(gridExtra)
library(lattice)
library(glmm)
library(bbmle)
library(TMB)
library(glmmTMB)
library(DHARMa)
library(dotwhisker)
library(parameters)
library(smplot2)
library(VGAM)
library(nnet)
library(GGally)
library(hrbrthemes)
library(ggsci)
library(see)
library(cowplot)
library(patchwork)

getwd()

different <- 'green'
general<- 'black'
local <- 'skyblue'

m9 <- 'royalblue'
tsb <-  'grey'
sb <- 'orange'

new_clone<- read.csv("final_data.csv",header=T)

clone<- new_clone %>%
  filter(!is.na(coculture1)) %>%
  mutate(monoculture= (monoculture*monoculture_density)) %>%
  mutate(coculture1 = (coculture1*coculture1_density)) %>%
  mutate(coculture2 = (coculture2*coculture2_density)) %>%
  mutate(coculture3 = (coculture3*coculture3_density)) %>%
  mutate(innoculation = innoculation*innoc_density)


## make longer so that reps included for pairs

clone<- clone %>%
  pivot_longer(cols = c(coculture1,coculture2,coculture3),
               names_to= 'Reps',
               values_to = 'Coculture')

## index for interaction proxy

clone <- clone %>%
  mutate(index= (Coculture-monoculture)/(Coculture+monoculture))

## density change 

clone <- clone %>%
  mutate(den_ch = (monoculture/innoculation)) %>%
  mutate(den_ch = (log(den_ch)/48))

##

model_clone<- clone %>%
  filter(!is.na(index))

model_clone$clone_sum <- paste(model_clone$clone.ID,model_clone$media)

model_clone$pair.ID<-as.factor(model_clone$pair.ID)
model_clone$clone.ID<-as.factor(model_clone$clone.ID)

## create variable for vial.id 

model_clone$vial.id <- paste(model_clone$pair.ID,model_clone$Reps,model_clone$media)
model_clone$vial.id <- as.factor(model_clone$vial.id)


pair_interaction<-model_clone %>%
  mutate(sign = ifelse(index>0,'+','-')) %>%
  group_by(vial.id,media,Reps,pair.ID,locality) %>%
  summarise(combo_sign = paste(sort(sign), collapse = "/"),
            combined_index = sum(index)) %>%
  ungroup()

pair_interaction$media_pair<-paste(pair_interaction$media,pair_interaction$pair.ID)

## figure 2

p1<- ggplot(data = clone, aes(x = media, y = den_ch, fill = media)) +
  sm_slope(group = factor(clone.ID),
           labels = c("Minimal Media","TSB","Soil Wash"),
           show_mean = TRUE,
           show_err = F) +
  #scale_fill_manual(values=c(m9,tsb,sb)) +
  geom_hline(yintercept = 0) +
  ylab("Density change (ln)") +
  scale_fill_manual(values=c(m9,tsb,sb))

p1
ggsave()


### S4

p4.1 <- ggplot(data = mean_int, aes(x = media, y = int_mean, fill = media)) +
  sm_slope(
    group = factor(clone_number),
    labels = c("M9", "TSB", "SB"),
    show_mean = TRUE,
    show_err = FALSE,
    point.params = list(
      size = 3.5,
      shape = 21
    )
  ) +
  geom_hline(yintercept = 0) +
  ylab("RII Index") +
  scale_fill_manual(values = c(m9, tsb, sb)) +
  theme(
    text = element_text(size = 25),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 25),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22)
  )

p4.1

### S3 


den_ch_av <-model_clone %>%
  group_by(clone.ID,media,den_ch,pair.ID) %>%
  summarise(mean_index = mean(index))

den_ch_av <- den_ch_av %>%
  filter(!is.na(den_ch))


den_ch_av$pred <- predict(d_c_i)

den_ch_ind <- ggplot(data=den_ch_av, aes(x=den_ch,y=mean_index,colour = media)) +
  geom_point()+
  geom_line(aes(y = pred), size = 1) +
  scale_colour_manual(values = c(m9,tsb,sb),,
                      labels = c('M9', 'TSB', 'Soil Wash')) +
  theme_bw() +
  ylab('Average RII index per clone') +
  xlab('Monoculture Growth (ln)')

den_ch_ind

den_ch_av$media <- as.factor(den_ch_av$media)

d_c_i <- lm(mean_index~den_ch*media,data = den_ch_av)
summary(d_c_i)
dci2<-update(d_c_i, .~. -den_ch:media)
anova(d_c_i,dci2)
simulateResiduals(dci2,plot = T)
summary(dci2)
dci3<-update(dci2, .~. -media)
summary(dci3)
anova(dci2,dci3)
simulateResiduals(dci3,plot=T)
summary(dci3)
dci4<- update(dci3, .~. -den_ch)
null_mod <- lm(mean_index ~1,data=den_ch_av)
summary(null_mod)






supp.labs <- c("Minimial Media", "TSB", 'Soil Wash')
names(supp.labs) <- c("A", "B",'C')


### take average for each pair rep


mean_int<-model_clone %>%
  group_by(clone_sum,locality,media,pair.ID,clone.ID) %>%
  summarise(int_mean = mean(index),
            sd_mean = sd(index)) %>%
  ungroup()

## figure 3

p2 <- ggplot(data=mean_int,aes(x=factor(locality),y=int_mean,colour = factor(locality))) +
  geom_boxplot(outlier.alpha=0) +
  geom_jitter(aes(size=sd_mean))+ #,
  #position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.775)) +
  scale_colour_manual(values=c(different,general,local)) +
  scale_size_continuous(range=c(1,3.5)) +
  theme_bw() +
  xlab('Locality') +
  ylab('RII Index') +
  labs(colour = "Locality")+
  facet_wrap(~media,labeller = labeller(media = supp.labs)) +
  theme(legend.position = 'none') +
  theme(text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

p2

## figure 6 


p6<- ggplot(data=mean_int, aes(x=media , y=abs(int_mean),colour=media)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter() +
  theme_bw() +
  scale_color_manual(values=c(m9,tsb,sb)) +
  theme(legend.position = 'none') +
  scale_x_discrete(labels=c('Minimal Media','TSB','Soil Wash')) +
  ylab('Mean Isolate RII Index (absolute values)') + xlab('Media') +
  theme(text = element_text(size=20))


p6



#####
##
#####


## modelling how media affects desnity change 

den_ch_lm <- lm(den_ch~media,data=clone)
summary(den_ch_lm)
simulateResiduals(den_ch_lm,plot=T)
den_ch_lm1<- lm(den_ch~1,data=clone)
anova(den_ch_lm,den_ch_lm1)

p_ch<- emmeans(den_ch_lm,specs = 'media',adjust='tukey')
pairs(p_ch)
emmeans(den_ch_lm, ~media,adjust='tukey')

###
###


### model on RII index

model2<- lmer(index~media*locality+(1|media:pair.ID)+(1|clone.ID),data = model_clone)
summary(model2)
DHARMa::simulateResiduals(model2,plot = T)

model3 <- lmer(index~media+locality+(1|media:pair.ID)+(1|clone.ID),data = model_clone)
anova(model2,model3) ## interaction not signif 

r.squaredGLMM(model2) 
r.squaredGLMM(model3) ## 37% variance explained by random variables 

model4 <- lmer(index~media+(1|media:pair.ID)+(1|clone.ID),data = model_clone)
anova(model3,model4) # locality not signif

r.squaredGLMM(model4) ## 31% variance explained by random variables 
summary(model4)

model5<- lmer(index~1+(1|media:pair.ID)+(1|clone.ID),data = model_clone)
anova(model4,model5) ## media v signif 

null_model<- lmer(index~1+(1|media:pair.ID)+(1|clone.ID),data = model_clone)
summary(null_model)
r.squaredGLMM(null_model)

null_model2<- lmer(index~1+(1|clone.ID),data = model_clone)
summary(null_model2)
r.squaredGLMM(null_model2)

ctr<-emmeans(model4,'media')
ctr
pair_ctr<-pairs(ctr)
pair_ctr

plot1<- emmip(model4, ~media, type="response", CIs=TRUE)
plot1

levels(ctr)

## Supplementary figure S2 

c_plot<- plot +
  scale_colour_manual(values=c(tsb,sb,m9)) +
  scale_y_discrete(labels=c('TSB','Soil Wash','Minimal Media')) +
  ylab('Media') +
  theme_bw()

c_plot

#### code for plot showing random effect

X <-model_clone%>%
  mutate(fit.m = predict(model4, re.form = NA),
         fit.c = predict(model4, re.form = NULL))

p_marginal <- X %>% 
  ggplot(aes(x = pair.ID, y= index)) +
  geom_point(pch =16, col = 'grey') +
  geom_line(aes(y=fit.m), col=1, size = 1) +
  coord_cartesian(ylim = c(-1,1)) +
  facet_wrap(vars(media))

p_marginal

pl<-model_clone %>%
  count(clone.ID, media)


##S1

p_conditional<- X %>%
  ggplot(aes(x = media, y = index)) +
  geom_point(pch = 1) +
  geom_line(aes(y = fit.c,colour = media), size = 2) +
  facet_wrap(vars(pair.ID)) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_colour_manual(values=c(m9,tsb,sb), labels = c ('Minimal Media','TSB','Soil Wash')) +
  theme_bw() 


## model of absolute values for strength of media 

var_index <- lmer(abs(index)~media*locality+(1|media:pair.ID)+(1|clone.ID),data=model_clone)
summary(var_index)

var_index2 <- update(var_index, .~. -media:locality)
anova(var_index,var_index2)

var_index3 <- update(var_index2, .~. -locality)
anova(var_index2,var_index3)

var_index4 <- update(var_index3, .~. -media)
anova(var_index3,var_index4)

summary(var_index3)
r.squaredGLMM(var_index3)

null_model<- lmer(index~1+(1|media:pair.ID)+(1|clone.ID),data = model_clone)
summary(null_model)
r.squaredGLMM(null_model)

null_model2<- lmer(abs(index)~1+(1|clone.ID),data = model_clone)
summary(null_model2)
r.squaredGLMM(null_model2)

null_model3<- lmer(abs(index)~1+(1|media:pair.ID)+(1|clone.ID),data = model_clone)
summary(null_model3)
r.squaredGLMM(null_model3)


#### multinomial on interaction type

pair_interaction$media<- as.factor(pair_interaction$media)
pair_interaction$combo_sign<- as.factor(pair_interaction$combo_sign)

pair_interaction <- pair_interaction[-c(55,58,61,91,97,138,141),]
pair_interaction <- pair_interaction[-c(90,130),]

pair_interaction$combo_sign <- droplevels(pair_interaction$combo_sign)
levels(pair_interaction$combo_sign)

### running multinomial, where our comparison category is -/- as this is the most common category 

# Run a "only intercept" model

pair_interaction$locality<-as.factor(pair_interaction$locality)

library(forcats)

pair_interaction <- pair_interaction %>%
  mutate(media = fct_recode(media,
                            "Minimal Media" = "A",
                            "TSB" = "B",
                            "Soil Wash" = "C"))



OIM <- multinom(combo_sign ~ 1, data = pair_interaction)
summary(OIM)

fit_full <- multinom(combo_sign ~ media+locality, data = pair_interaction)
summary(fit_full)

fit_full1<- multinom(combo_sign ~ media, data = pair_interaction)

fit_full2<- multinom(combo_sign ~ 1, data = pair_interaction)

anova(fit_full1,fit_full2)

anova(fit_full,fit_full1)
summary(fit_full1)


z <- summary(fit_full1)$coefficients/summary(fit_full1)$standard.errors
z

p_value <- (1 - pnorm(abs(z), 0, 1)) * 2
p_value

anova(OIM,fit_full1)

library(ggeffects)

p<-ggeffect(fit_full1, terms = "media") %>%
  plot()
p


response.labs <- c("Competitive (-/-)", "Exploitative (+/-)", "Cooperative (+/+)")
names(response.labs) <- c("X...", "X....1",
                          "X....2")


p<-p+ ylab('Predicted probability') +
  ggtitle('') +
  xlab('Media')

p<- p +
  scale_colour_manual(values = c(m9,tsb,sb)) +
  facet_wrap(~ response.level,
             labeller = labeller(response.level = response.labs)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) 

p



df <- ggeffect(fit_full1, terms = "media") |> 
  as.data.frame()


### Final figure 5 


p_2 <- ggplot(df, aes(x = x, y = predicted, colour = x)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width = 0.15
  ) +
  facet_wrap(
    ~ response.level,
    labeller = labeller(response.level = response.labs)
  ) +
  scale_colour_manual(values = c(m9, tsb, sb)) +
  ylab("Predicted probability") +
  xlab("Media") +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    legend.position = "none"
  ) +
  theme_bw()+
  theme(legend.position = "none")

p_2




