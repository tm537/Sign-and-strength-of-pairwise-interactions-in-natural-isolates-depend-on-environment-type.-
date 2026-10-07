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

new_clone <- new_clone[-c(15,16,61,62,63,64,65,66,103,104,105,106,107,108), ] # remove isolates where monoculture not available 

new_clone %>% count(media)



str(new_clone)

new_clone$monoculture <- as.numeric(new_clone$monoculture)
new_clone$monoculture_density <- as.numeric(new_clone$monoculture_density)
new_clone$mono_plating <- as.numeric(new_clone$mono_plating)

## cfu/mL transformation

clone<- new_clone %>%
  mutate(monoculture= (monoculture*monoculture_density*mono_plating)) %>%
  mutate(coculture1 = (coculture1*coculture1_density*40)) %>%
  mutate(coculture2 = (coculture2*coculture2_density*40)) %>%
  mutate(coculture3 = (coculture3*coculture3_density*40)) %>%
  mutate(innoculation = innoculation*innoc_density*innoc_plating)


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
  mutate(den_ch = log10(den_ch))

clone <- clone %>%
  mutate(raw_dench = (monoculture/innoculation)) %>%
  mutate(den_ch_ln = log(raw_dench))

model_clone<- clone %>%
  filter(!is.na(index))

model_clone$clone_sum <- paste(model_clone$clone.ID,model_clone$media)
model_clone$pair.ID<-as.factor(model_clone$pair.ID)
model_clone$clone.ID<-as.factor(model_clone$clone.ID)

## create variable for vial.id 

model_clone$vial.id <- paste(model_clone$pair.ID,model_clone$Reps,model_clone$media)
model_clone$vial.id <- as.factor(model_clone$vial.id)


# amount of positive and negative interactions


model_clone %>%
  summarise(
    positive = sum(index > 0, na.rm = TRUE),
    negative = sum(index < 0, na.rm = TRUE)
  )

# create datarframe with the pair interaction (+/= , +/- , -/- and the summed interaction index for that pair)

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
  ylab("Density change log(10)") +
  scale_fill_manual(values=c(m9,tsb,sb))

p1
ggsave()

ggsave('figure2_final.png',plot=get_last_plot(),width=6,height=6)


### figure 3 

p4.1 <- ggplot(data = mean_int, aes(x = media, y = int_mean, fill = media)) +
  sm_slope(
    group = factor(clone.ID),
    labels = c("Minimal Media", "TSB", "Soil Wash"),
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

mean_int %>%
  count(media)

ggsave('Figure3.png',plot = get_last_plot(),height = 8, width = 8)


### S3 -- supplementary figure 3, showing graphical representation of modelled relationship 


den_ch_av <-model_clone %>%
  group_by(clone.ID,media,den_ch,pair.ID) %>%
  summarise(mean_index = mean(index)) %>%
  ungroup()

den_ch_av <- den_ch_av %>%
  filter(!is.na(den_ch))


model_clone$pred <- predict(d_c_i)

den_ch_ind <- ggplot(data=den_ch_av, aes(x=den_ch,y=mean_index,colour = media)) +
  geom_point()+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1,alpha=0.2) +
  scale_colour_manual(values = c(m9,tsb,sb),,
                      labels = c('M9', 'TSB', 'Soil Wash')) +
  theme_bw() +
  ylab('Average RII index per clone') +
  xlab('Monoculture Density Change (log10)')

den_ch_ind

ggsave('quad_plot.png', plot = get_last_plot(), height = 5, width = 7)


## modelliing above realtionship, using raw data - not averaged index

sum(is.na(model_clone$den_ch))
sum(is.na(model_clone$media))
sum(is.na(model_clone$index))

d_c_i <- lmer (index~poly(den_ch, 2)*media+ (1|clone.ID:media) ,data = model_clone)
summary(d_c_i)
simulationOutput <- simulateResiduals(fittedModel = d_c_i, plot = FALSE)
plotResiduals(simulationOutput, form = model_clone$den_ch)
plotResiduals(simulationOutput, form = model_clone$media)
plotResiduals(simulationOutput, form = model_clone$clone_sum)
plotResiduals(simulationOutput, form = model_clone$locality)


## model overpredicting majority of residuals, curevd relationship 

plot(model_clone$den_ch, model_clone$index, 
     col = as.factor(model_clone$media), pch = 16,
     xlab = "den_ch", ylab = "index")
legend("topright", legend = levels(as.factor(model_clone$media)), 
       col = 1:length(levels(as.factor(model_clone$media))), pch = 16)


d_c_i2 <- lmer(index ~ poly(den_ch, 2) * media + (1|clone.ID:media), data = model_clone)
simulationOutput <- simulateResiduals(fittedModel = d_c_i2, plot = FALSE)
plotResiduals(simulationOutput, form = model_clone$den_ch)
plotResiduals(simulationOutput, form = model_clone$media)
plotResiduals(simulationOutput, form = model_clone$clone.ID)
testDispersion(simulationOutput)
simulateResiduals(d_c_i2,plot=T)
outliers <- outliers(simulationOutput)
outliers
model_clone[outliers, ]
row_id <- outliers(simulationOutput)
predict(d_c_i2)[row_id]      # expected
model_clone$index[row_id]    # observed

dci2.1 <- lmer(index ~ poly(den_ch, 2) + media + (1|clone.ID:media), data = model_clone)
simulateResiduals(dci2.1,plot=T)


anova(d_c_i2,dci2.1) # interaction between quadratic term and media is significant 

joint_tests(d_c_i2,by = 'media')
trends <- emtrends(d_c_i2, ~ media, var = "den_ch")
trends

boot_conf <- confint(d_c_i2, method="boot", nsim=1000)
boot_conf

preds <- ggpredict(dci2.1, terms = c("den_ch [all]",'media'))
plot(preds)


## run analysis without line 254 outlier 

model_clone_254 <- model_clone[-254, ]

d_c_i2_254 <- lmer(index ~ poly(den_ch, 2) * media + (1|clone.ID), data = model_clone_254)
simulationOutput <- simulateResiduals(fittedModel = d_c_i2_254, plot = FALSE)
plotResiduals(simulationOutput, form = model_clone$den_ch)
plotResiduals(simulationOutput, form = model_clone$media)
simulateResiduals(d_c_i2_254,plot=T)
outliers <- outliers(simulationOutput)











supp.labs <- c("Minimial Media", "TSB", 'Soil Wash')
names(supp.labs) <- c("A", "B",'C')


### take average for each pair rep


mean_int<-model_clone %>%
  group_by(clone_sum,locality,media,pair.ID,clone.ID) %>%
  summarise(int_mean = mean(index),
            sd_mean = sd(index)) %>%
  ungroup()

## supplementary figure 4

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

ggsave('supp_figure.png', plot = get_last_plot(), height = 5 , width = 7)

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

ggsave('abs_plot.png',plot=get_last_plot() , height=6,width=6)


#####
##
#####


## modelling how media affects desnity change 

den_ch_lm <- lm(den_ch~media,data=clone)
summary(den_ch_lm)
simulateResiduals(den_ch_lm,plot=T)
simout <- simulateResiduals(fittedModel = den_ch_lm,plot=FALSE)
plotResiduals(simout,form = clone$media)
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
r.squaredGLMM(model3) ## 21% variance explained by random variables 

model4 <- lmer(index~media+(1|media:pair.ID)+(1|clone.ID),data = model_clone)
anova(model3,model4) # locality not signif

r.squaredGLMM(model4) ## 20.5% variance explained by random variables 
summary(model4)

model5<- lmer(index~1+(1|media:pair.ID)+(1|clone.ID),data = model_clone)
anova(model4,model5) ## media v signif 

r.squaredGLMM(model5)

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
p_conditional

ggsave('supp_2.png',plot=get_last_plot(),height=6,width=6)


## model of absolute values for strength of media 

var_index <- lmer(abs(index)~media*locality+(1|media:pair.ID)+(1|clone.ID),data=model_clone)
summary(var_index)
simulateResiduals(var_index,plot=T)

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


emms <- emmeans(var_index3, ~ media)
emms
pairs(emms)


#### multinomial on interaction type

pair_interaction$media<- as.factor(pair_interaction$media)
pair_interaction$combo_sign<- as.factor(pair_interaction$combo_sign)


pair_interaction$combo_sign <- droplevels(pair_interaction$combo_sign)
levels(pair_interaction$combo_sign)

### running multinomial, where our comparison category is -/- in TSB as this is the most common category 

# Run a "only intercept" model

pair_interaction %>%
  count(media,combo_sign)


pair_interaction$locality<-as.factor(pair_interaction$locality)

library(forcats)

pair_interaction <- pair_interaction %>%
  mutate(media = fct_recode(media,
                            "Minimal Media" = "A",
                            "TSB" = "B",
                            "Soil Wash" = "C"))


pair_interaction$media <- relevel(pair_interaction$media, ref = "TSB")
fit_full1 <- multinom(combo_sign ~ media, data = pair_interaction)


OIM <- multinom(combo_sign ~ 1, data = pair_interaction)
summary(OIM)

fit_full <- multinom(combo_sign ~ media+locality, data = pair_interaction)
summary(fit_full)

fit_full1<- multinom(combo_sign ~ media, data = pair_interaction)

fit_full2<- multinom(combo_sign ~ 1, data = pair_interaction)

anova(fit_full,fit_full1,fit_full2)

summary(fit_full1)


summary(fit_full1)$coefficients

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

ggsave("LRT.png",plot=get_last_plot(),height=5,width=7)





