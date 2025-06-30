## ----setup, include=FALSE-----------------------------------------------------

knitr::knit_hooks$set(purl = knitr::hook_purl)

library(ggplot2)
library(lmerTest) 
library(psych)
library(tidyverse)
library(table1)
library(survival)
library(dplyr)
library(stats)
library(emmeans)
library(ggpubr)
library(parameters)
library(ggthemes)
library(reactable)
library(car)
library(jtools)
library(rstatix)
library(interactions)


## ----Create demographics table------------------------------------------------

summary(total_df)

demo <- subset(total_df, select = c('user_id', 'age', 'pronoun_1', 'indigenous_status', 'treatment_stage'))

demo[demo == "NA"] <- "UNKNOWN"
demo[demo == "NULL"] <- "UNKNOWN"
demo[demo == ""] <- "UNKNOWN"
demo[demo == "NotStated"] <- "UNKNOWN"


demo$pronoun_1 <- factor(demo$pronoun_1, levels =c('female','male','neutral', 'else'), labels =c('She/Her', 'He/Him', 'They/Them', 'Else'))
demo$indigenous_status <- factor(demo$indigenous_status, levels = c('NonIndigenous', 'Aboriginal',  'UNKNOWN', 'TorresStraitIslander', 'AboriginalAndTorresStraitIslander'),
                                labels = c('Non-Indigenous', 'Aboriginal',  'Unknown', 'Torres Strait Islander', 'Aboriginal and Torres Strait Islander'))
demo$treatment_stage <- factor(demo$treatment_stage, levels = c('WAITING_F2F','RECEIVING_F2F', 'UNKNOWN', 'APPROACHING_DISCHARGE','SUBTHRESHOLD_DISCHARGE','DISCHARGED'),
                                labels = c('Waiting','Receiving', 'Unknown', 'Approaching Discharge','Subthreshold Discharge','Discharged'))

label(demo$age) <- "Age"
label(demo$pronoun_1) <- "Pronouns"
label(demo$indigenous_status) <- "Indigenous status"
label(demo$treatment_stage) <- "Treatment stage"


units(demo$age) <- "years"



labels <- list(
  variables = list(age="Age (years)",
                   pronoun_1="Pronouns, number (%)",
                   indigenous_status="Indigenous status, number (%)",
                   treatment_stage="Treatment stage, number (%)"

))

my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), c("", "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}

my.render.cat <- function(x) {
  c("",sapply(stats.default(x),function(y) with(y,
                                                sprintf("%d (%0.0f %%)", FREQ, PCT))))
}

strata <- list(Total=demo)

table1(strata, labels, render.missing=NULL,
       render.continuous = my.render.cont, render.categorical = my.render.cat)



## ----Summary statistics-------------------------------------------------------

total_df |> summarise(mean_var1 = mean(engaged_duration, na.rm=TRUE), sd_var1 = sd(engaged_duration, na.rm=TRUE), mean_var2 = mean(active_days, na.rm=TRUE), sd_var2 = sd(active_days, na.rm=TRUE), mean_var3 = mean(duration, na.rm=TRUE), sd_var3 = sd(duration, na.rm=TRUE), mean_var4 = mean(sessions, na.rm=TRUE), sd_var4 = sd(sessions, na.rm=TRUE), mean_var5 = mean(reactions, na.rm=TRUE), sd_var5 = sd(reactions, na.rm=TRUE), mean_var6 = mean(message_number, na.rm=TRUE), sd_var6 = sd(message_number, na.rm=TRUE))


k10_1 <- k10_baseline |> 
  subset(select = c("user_id", "k10_total")) |>
  dplyr::rename(baseline = k10_total)
summary(k10_1)

k10_2 <- k10_week6 |> 
  subset(select=c("user_id","k10_total")) |>
  dplyr::rename(week_6 = k10_total) 
summary(k10_2)

k10_3 <- k10_week12 |> 
  subset(select=c("user_id","k10_total")) |>
  dplyr::rename(week_12 = k10_total) 
summary(k10_3)

k10_list <- list(k10_1, k10_2, k10_3)
K10 <- k10_list %>% reduce(full_join, by='user_id')

K10 <- K10 |> 
  gather(key = "time", value = "k10_total", baseline, week_6, week_12) |>
  convert_as_factor(user_id,time)

k10.aov <- anova_test(data = K10, dv = k10_total, wid = user_id, within = time)
get_anova_table(k10.aov)


baseline_mean <- k10_1 |>  dplyr::summarise(mean=mean(baseline), sd=sd(baseline), freq=n())
week6_mean <- k10_2 |>  dplyr::summarise(mean=mean(week_6), sd=sd(week_6), freq=n())
week12_mean <- k10_3 |>  dplyr::summarise(mean=mean(week_12), sd=sd(week_12), freq=n())
k10_mean <- bind_rows(baseline_mean,week6_mean,week12_mean, .id = "time") 


ggplot(k10_mean, aes(x=time, y=mean, fill=time)) +         
  geom_bar(stat="identity",color="black",position=position_dodge()) +
  geom_errorbar( aes(x=time, ymin=mean-sd, ymax=mean+sd), width=0.2, colour="black", alpha=0.9, size=1.3) +
  scale_x_discrete(labels = c('Onboarding','Week 6','Week 12')) +
  labs(title=" ", x=" ", y = "K10 score") + theme_classic() + 
  theme(axis.title.y = element_text(size = rel(2), angle = 90, margin=margin(0,20,0,0)), axis.title.x = element_text(size = rel(3), angle = 00, margin=margin(20,0,0,0))) +
  theme(plot.title = element_text(size = rel(3), angle = 00)) +
  theme(axis.text.x = element_text(angle = 00, size=18), axis.text.y = element_text(angle = 00, size=20)) +
  theme(legend.position = "none") +
  theme(axis.line.x = element_line(color="black", linewidth = 1.2), axis.line.y = element_line(color="black", linewidth = 1.2))+
  theme(axis.ticks.y=element_line(size=(1.5)), axis.ticks.x=element_line(size=(1.5)), axis.ticks.length=unit(0.4, "cm")) +
  scale_fill_manual(values=c("#669d62", "#3c7c3d", "#1f5b25"))

k10_wide <- reshape(K10, idvar = "user_id", timevar = "time", direction = "wide")



## ----Engagement metrics-------------------------------------------------------

mean_weeklyTime <- time |> group_by(week_number) |>
  rename(duration = time_on_platform_min) |>
  summarise(mean=mean(duration), sd=sd(duration)) 
mean_weeklySessions <- session_number |> group_by(week_number) |>
  rename(sessions = session.Count.per.week) |>
  summarise(mean=mean(sessions), sd=sd(sessions)) 
mean_weeklyReactions <- reactions |> group_by(week_number) |>
  dplyr::rename(reactions = number.of.reactions.made.by.YP) |>
  summarise(mean=mean(reactions), sd=sd(reactions)) 
mean_weeklyClinical <- message_number |> group_by(week_number) |>
  summarise(mean=mean(message_number), sd=sd(message_number)) 


mean_timePlot <- ggplot(mean_weeklyTime, aes(x=week_number, y=mean)) +
  geom_line(color="#7f793c",linewidth=1.3,linetype = "solid") + geom_point(color="#7f793c",size=3) +
  labs(x="Week",y="Time on platform (minutes)") +
  scale_x_continuous(limits = c(1,12), breaks = c(2,4,6,8,10,12)) +
  scale_y_continuous(limits = c(0,30), breaks = c(0,5,10,15,20,25,30)) +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 15, angle = 90, margin=margin(0,20,0,0)), axis.title.x = element_text(size = 20, angle = 00, margin=margin(20,0,0,0))) +
  theme(plot.title = element_text(size = rel(3), angle = 00)) +
  theme(axis.text.x = element_text(angle = 00, size=20), axis.text.y = element_text(angle = 00, size=20)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
  theme(axis.line.x = element_line(color="black", linewidth = 1.2), axis.line.y = element_line(color="black", linewidth = 1.2)) +
  theme(axis.ticks.y=element_line(linewidth=(1.5)), axis.ticks.x=element_line(size=(1.5)), axis.ticks.length=unit(0.4, "cm")) 



mean_sessionsPlot <- ggplot(mean_weeklySessions, aes(x=week_number, y=mean)) +
  geom_line(color="#565c33",size=1.3, linetype = "dashed") + geom_point(color="#565c33",size=3) +
  labs(x="Week",y="Sessions on platform") +
  scale_x_continuous(limits = c(1,12), breaks = c(2,4,6,8,10,12)) +
  scale_y_continuous(limits = c(0,6), breaks = c(0,1,2,3,4,5,6)) +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 15, angle = 90, margin=margin(0,20,0,0)), axis.title.x = element_text(size = 20, angle = 00, margin=margin(20,0,0,0))) +
  theme(plot.title = element_text(size = rel(3), angle = 00)) +
  theme(axis.text.x = element_text(angle = 00, size=20), axis.text.y = element_text(angle = 00, size=20)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
  theme(axis.line.x = element_line(color="black", linewidth = 1.2), axis.line.y = element_line(color="black", linewidth = 1.2)) +
  theme(axis.ticks.y=element_line(size=(1.5)), axis.ticks.x=element_line(size=(1.5)), axis.ticks.length=unit(0.4, "cm"))



mean_reactionsPlot <- ggplot(mean_weeklyReactions, aes(x=week_number, y=mean)) +
  geom_line(color="#184948",size=1.3, linetype = "twodash") + geom_point(color="#184948",size=3) +
  labs(x="Week",y="Reactions on social network") +
  scale_x_continuous(limits = c(1,12), breaks = c(2,4,6,8,10,12)) +
  scale_y_continuous(limits = c(0,4), breaks = c(0,1,2,3,4)) +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 15, angle = 90, margin=margin(0,20,0,0)), axis.title.x = element_text(size = 20, angle = 00, margin=margin(20,0,0,0))) +
  theme(plot.title = element_text(size = rel(3), angle = 00)) +
  theme(axis.text.x = element_text(angle = 00, size=20), axis.text.y = element_text(angle = 00, size=20)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
  theme(axis.line.x = element_line(color="black", linewidth = 1.2), axis.line.y = element_line(color="black", linewidth = 1.2)) +
  theme(axis.ticks.y=element_line(size=(1.5)), axis.ticks.x=element_line(size=(1.5)), axis.ticks.length=unit(0.4, "cm"))


mean_clinicalPlot <- ggplot(mean_weeklyClinical, aes(x=week_number, y=mean)) +
  geom_line(color="#022a2a",size=1.3, linetype = "twodash") + geom_point(color="#022a2a",size=3) +
  labs(x="Week",y="Messages to clinical team") +
  scale_x_continuous(limits = c(1,12), breaks = c(2,4,6,8,10,12)) +
  scale_y_continuous(limits = c(0,4), breaks = c(0,1,2,3,4)) +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 15, angle = 90, margin=margin(0,20,0,0)), axis.title.x = element_text(size = 20, angle = 00, margin=margin(20,0,0,0))) +
  theme(plot.title = element_text(size = rel(3), angle = 00)) +
  theme(axis.text.x = element_text(angle = 00, size=20), axis.text.y = element_text(angle = 00, size=20)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
  theme(axis.line.x = element_line(color="black", linewidth = 1.2), axis.line.y = element_line(color="black", linewidth = 1.2)) +
  theme(axis.ticks.y=element_line(size=(1.5)), axis.ticks.x=element_line(size=(1.5)), axis.ticks.length=unit(0.4, "cm"))


ggarrange(mean_timePlot,mean_sessionsPlot,mean_reactionsPlot,mean_clinicalPlot,
          labels = c("A", "B", "C", "D"),
          ncol = 1, nrow = 4)



## ----Character strengths -- prepare data--------------------------------------
strengths_ID <- total_df[, 10:69]



kindness <- strengths_ID[,1:3] |> mutate_if(is.character, as.numeric)
kindness$kindness_total = rowSums(kindness)
kindness = tibble::rownames_to_column(kindness)

teamwork <- strengths_ID[,4:6] |> mutate_if(is.character, as.numeric)
teamwork$teamwork_total = rowSums(teamwork)
teamwork = tibble::rownames_to_column(teamwork)

creativity <- strengths_ID[,7:9] |> mutate_if(is.character, as.numeric)
creativity$creativity_total = rowSums(creativity)
creativity = tibble::rownames_to_column(creativity)

discretion <- strengths_ID[,10:12] |> mutate_if(is.character, as.numeric)
discretion$discretion_total = rowSums(discretion)
discretion = tibble::rownames_to_column(discretion)

humour <- strengths_ID[,13:15] |> mutate_if(is.character, as.numeric)
humour$humour_total = rowSums(humour)
humour = tibble::rownames_to_column(humour)

forgiveness <- strengths_ID[,16:18] |> mutate_if(is.character, as.numeric)
forgiveness$forgiveness_total = rowSums(forgiveness)
forgiveness = tibble::rownames_to_column(forgiveness)

honesty <- strengths_ID[,19:21] |> mutate_if(is.character, as.numeric)
honesty$honesty_total = rowSums(honesty)
honesty = tibble::rownames_to_column(honesty)

love_of_learning <- strengths_ID[,22:24] |> mutate_if(is.character, as.numeric)
love_of_learning$love_of_learning_total = rowSums(love_of_learning)
love_of_learning = tibble::rownames_to_column(love_of_learning)

inspiration <- strengths_ID[,25:27] |> mutate_if(is.character, as.numeric)
inspiration$inspiration_total = rowSums(inspiration)
inspiration = tibble::rownames_to_column(inspiration)

leadership <- strengths_ID[,28:30] |> mutate_if(is.character, as.numeric)
leadership$leadership_total = rowSums(leadership)
leadership = tibble::rownames_to_column(leadership)

perspective <- strengths_ID[,31:33] |> mutate_if(is.character, as.numeric)
perspective$perspective_total = rowSums(perspective)
perspective = tibble::rownames_to_column(perspective)

gratitude <- strengths_ID[,34:36] |> mutate_if(is.character, as.numeric)
gratitude$gratitude_total = rowSums(gratitude)
gratitude = tibble::rownames_to_column(gratitude)

courage <- strengths_ID[,37:39] |> mutate_if(is.character, as.numeric)
courage$courage_total = rowSums(courage)
courage = tibble::rownames_to_column(courage)

curiosity <- strengths_ID[,40:42] |> mutate_if(is.character, as.numeric)
curiosity$curiosity_total = rowSums(curiosity)
curiosity = tibble::rownames_to_column(curiosity)

social_intelligence <- strengths_ID[,43:45] |> mutate_if(is.character, as.numeric)
social_intelligence$social_intelligence_total = rowSums(social_intelligence)
social_intelligence = tibble::rownames_to_column(social_intelligence)

fairness <- strengths_ID[,46:48] |> mutate_if(is.character, as.numeric)
fairness$fairness_total = rowSums(fairness)
fairness = tibble::rownames_to_column(fairness)

hope <- strengths_ID[,49:51] |> mutate_if(is.character, as.numeric)
hope$hope_total = rowSums(hope)
hope = tibble::rownames_to_column(hope)

enthusiasm <- strengths_ID[,52:54] |> mutate_if(is.character, as.numeric)
enthusiasm$enthusiasm_total = rowSums(enthusiasm)
enthusiasm = tibble::rownames_to_column(enthusiasm)

love <- strengths_ID[,55:57] |> mutate_if(is.character, as.numeric)
love$love_total = rowSums(love)
love = tibble::rownames_to_column(love)

perseverance <- strengths_ID[,58:60] |> mutate_if(is.character, as.numeric)
perseverance$perseverance_total = rowSums(perseverance)
perseverance = tibble::rownames_to_column(perseverance)

# Grid arrange

par(mfrow = c(5,4))
hist(courage$courage_total, col = "#9a133d", xlab = "", ylab = "Frequency", main = "Courage", ylim = c(0, 1200))
axis(side=1, at=seq(3, 15, 1), labels=seq(3, 15, 1))
hist(creativity$creativity_total, col = "#b93961", xlab = "", ylab = "Frequency", main = "Creativity", ylim = c(0, 1200))
axis(side=1, at=seq(3, 15, 1), labels=seq(3, 15, 1))
hist(curiosity$curiosity_total, col = "#d8527c", xlab = "", ylab = "Frequency", main = "Curiosity", ylim = c(0, 1200))
axis(side=1, at=seq(3, 15, 1), labels=seq(3, 15, 1))
hist(discretion$discretion_total, col = "#f28aaa", xlab = "", ylab = "Frequency", main = "Discretion", ylim = c(0, 1200))
axis(side=1, at=seq(3, 15, 1), labels=seq(3, 15, 1))
hist(enthusiasm$enthusiasm_total, col = "#f9b4c9", xlab = "", ylab = "Frequency", main = "Enthusiasm", ylim = c(0, 1200))
axis(side=1, at=seq(3, 15, 1), labels=seq(3, 15, 1))
hist(fairness$fairness_total, col = "#f9e0e8", xlab = "", ylab = "Frequency", main = "Fairness", ylim = c(0, 1200))
axis(side=1, at=seq(3, 15, 1), labels=seq(3, 15, 1))
hist(forgiveness$forgiveness_total, col = "#ffffff", xlab = "", ylab = "Frequency", main = "Forgiveness", ylim = c(0, 1200))
axis(side=1, at=seq(3, 15, 1), labels=seq(3, 15, 1))
hist(gratitude$gratitude_total, col = "#eaf3ff", xlab = "", ylab = "Frequency", main = "Gratitude", ylim = c(0, 1200))
axis(side=1, at=seq(3, 15, 1), labels=seq(3, 15, 1))
hist(honesty$honesty_total, col =  "#c5daf6", xlab = "", ylab = "Frequency", main = "Honesty", ylim = c(0, 1200))
axis(side=1, at=seq(3, 15, 1), labels=seq(3, 15, 1))
hist(hope$hope_total, col = "#a1c2ed", xlab = "", ylab = "Frequency", main = "Hope", ylim = c(0, 1200))
axis(side=1, at=seq(3, 15, 1), labels=seq(3, 15, 1))
hist(humour$humour_total, col = "#6996e3", xlab = "", ylab = "Frequency", main = "Humour", ylim = c(0, 1200))
axis(side=1, at=seq(3, 15, 1), labels=seq(3, 15, 1))
hist(inspiration$inspiration_total, col = "#4060c8", xlab = "", ylab = "Frequency", main = "Inspiration", ylim = c(0, 1200))
axis(side=1, at=seq(3, 15, 1), labels=seq(3, 15, 1))
hist(kindness$kindness_total, col = "#1a318b", xlab = "", ylab = "Frequency", main = "Kindness", ylim = c(0, 1200))
axis(side=1, at=seq(3, 15, 1), labels=seq(3, 15, 1))
hist(leadership$leadership_total, col = "#e7e5cc", xlab = "", ylab = "Frequency", main = "Leadership", ylim = c(0, 1200))
axis(side=1, at=seq(3, 15, 1), labels=seq(3, 15, 1))
hist(love$love_total, col = "#c2d6a4", xlab = "", ylab = "Frequency", main = "Love", ylim = c(0, 1200))
axis(side=1, at=seq(3, 15, 1), labels=seq(3, 15, 1))
hist(love_of_learning$love_of_learning_total, col = "#9cc184", xlab = "", ylab = "Frequency", main = "Love of Learning", ylim = c(0, 1200))
axis(side=1, at=seq(3, 15, 1), labels=seq(3, 15, 1))
hist(perseverance$perseverance_total, col = "#669d62", xlab = "", ylab = "Frequency", main = "Perseverance", ylim = c(0, 1200))
axis(side=1, at=seq(3, 15, 1), labels=seq(3, 15, 1))
hist(perspective$perspective_total, col = "#3c7c3d", xlab = "", ylab = "Frequency", main = "Perspective", ylim = c(0, 1200))
axis(side=1, at=seq(3, 15, 1), labels=seq(3, 15, 1))
hist(social_intelligence$social_intelligence_total, col = "#1f5b25", xlab = "", ylab = "Frequency", main = "Social Intelligence", ylim = c(0, 1200))
axis(side=1, at=seq(3, 15, 1), labels=seq(3, 15, 1))
hist(teamwork$teamwork_total, col = "#1e3d14", xlab = "", ylab = "Frequency", main = "Teamwork", ylim = c(0, 1200))
axis(side=1, at=seq(3, 15, 1), labels=seq(3, 15, 1))

## ----Character strengths -- Factor analysis-----------------------------------

# Factorability of the data
# Overall KMO MSA is 0.95
KMO(r=cor(strengths_ID))

# Bartlett’s Test of Sphericity
cor_matrix <- cor(strengths_ID)
#check_sphericity_bartlett(cor_matrix, n=nrow(strengths_ID))

# Plot
fafitfree <- fa(strengths_ID, nfactors = ncol(strengths_ID), rotate = "oblimin")
n_factors <- length(fafitfree$e.values)
scree     <- data.frame(
  Factor_n =  as.factor(1:n_factors), 
  Eigenvalue = fafitfree$e.values)

ggplot(scree, aes(x = Factor_n, y = Eigenvalue, group = 1)) + 
  geom_point(size=2.5) + geom_line(size=1) +
  xlab("Number of factors") +
  ylab("Eigenvalue") +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 20, angle = 90, margin=margin(0,20,0,0)), axis.title.x = element_text(size = 20, angle = 00, margin=margin(20,0,0,0))) +
  theme(plot.title = element_text(size = rel(3), angle = 00)) +
  theme(axis.text.x = element_text(angle = 00, size=10), axis.text.y = element_text(angle = 00, size=15)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
  theme(axis.line.x = element_line(color="black", linewidth = 1.2), axis.line.y = element_line(color="black", linewidth = 1.2)) +
  theme(axis.ticks.y=element_line(size=(1.5)), axis.ticks.x=element_line(size=(1.5)), axis.ticks.length=unit(0.4, "cm")) 

# Decide on rotation
Nfacts <- 3
#fit <- factanal(strengths_ID, Nfacts, rotation = "oblimin")
#print(fit, digits=2, cutoff=0.3, sort=TRUE) # If enough correlations are >0.3, use oblique rotation. If not, change to orthogonal

# Orthogonal rotation
# Run final factor analysis
fa <- fa(strengths_ID, nfactors = 3, n.obs = nrow(strengths_ID), rotate = "equamax", fm="ml", scores="regression")
loads <- fa$loadings
fa.diagram(loads)

# Save loadings
fa_loadings <- data.frame(unclass(fa$loadings))
print(fa$loadings,cutoff = 0.3)

# Save factor scores
fa.scores <- factor.scores(x=strengths_ID, method="Anderson", f=fa)
scores = data.frame(fa.scores$scores)
colnames(scores) <- c("Factor1", "Factor2", "Factor3")

character_dataset <- dplyr::bind_cols(total_df, scores)
character_factors <- subset(character_dataset, select = c("user_id", "Factor1", "Factor2", "Factor3"))

# Plot loadings

fa_loadings$Subscale <- row.names(fa_loadings)
fa_loadings$Subscale[1:3] <- "Kindness"
fa_loadings$Subscale[4:6] <- "Teamwork"
fa_loadings$Subscale[7:9] <- "Creativity"
fa_loadings$Subscale[10:12] <- "Discretion"
fa_loadings$Subscale[13:15] <- "Humour"
fa_loadings$Subscale[16:18] <- "Forgiveness"
fa_loadings$Subscale[19:21] <- "Honesty"
fa_loadings$Subscale[22:24] <- "Love of Learning"
fa_loadings$Subscale[25:27] <- "Inspiration"
fa_loadings$Subscale[28:30] <- "Leadership"
fa_loadings$Subscale[31:33] <- "Perspective"
fa_loadings$Subscale[34:36] <- "Gratitude"
fa_loadings$Subscale[37:39] <- "Courage"
fa_loadings$Subscale[40:42] <- "Curiosity"
fa_loadings$Subscale[43:45] <- "Social Intelligence"
fa_loadings$Subscale[46:48] <- "Fairness"
fa_loadings$Subscale[49:51] <- "Hope"
fa_loadings$Subscale[52:54] <- "Enthusiasm"
fa_loadings$Subscale[55:57] <- "Love"
fa_loadings$Subscale[58:60] <- "Perseverance"
fa_loadings$qn_num <- row.names(fa_loadings)

F1 <-ggplot(fa_loadings, aes(x = qn_num, y = ML1, fill = Subscale), position = position_dodge2(preserve = "single")) + 
  geom_col(width = 1.5) +
  labs(color = NULL) +
  ylab("Loadings") +
  ggtitle("Factor 1") +
  scale_y_continuous(breaks=c(-1.0,-0.5, 0, 0.5, 1.0)) +
  theme_hc()+ scale_colour_hc() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_manual(values=c("#9a133d", "#b93961", "#d8527c", "#f28aaa", "#f9b4c9", "#f9e0e8", "#ffffff", "#eaf3ff", "#c5daf6", "#a1c2ed", "#6996e3", "#4060c8", "#1a318b", "#e7e5cc", "#c2d6a4", "#9cc184", "#669d62", "#3c7c3d", "#1f5b25", "#1e3d14"))

F2 <-ggplot(fa_loadings, aes(x = qn_num, y = ML2, fill = Subscale), position = position_dodge2(preserve = "single")) + 
  geom_col(width = 1.5) +
  labs(color = NULL) +
  ylab("Loadings") +
  ggtitle("Factor 2") +
  scale_y_continuous(breaks=c(-1.0,-0.5, 0, 0.5, 1.0)) +
  theme_hc()+ scale_colour_hc() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_manual(values=c("#9a133d", "#b93961", "#d8527c", "#f28aaa", "#f9b4c9", "#f9e0e8", "#ffffff", "#eaf3ff", "#c5daf6", "#a1c2ed", "#6996e3", "#4060c8", "#1a318b", "#e7e5cc", "#c2d6a4", "#9cc184", "#669d62", "#3c7c3d", "#1f5b25","#1e3d14"))


F3 <-ggplot(fa_loadings, aes(x = qn_num, y = ML3, fill = Subscale), position = position_dodge2(preserve = "single")) + 
  geom_col(width = 1.5) +
  labs(color = NULL) +
  ylab("Loadings") +
  ggtitle("Factor 3") +
  scale_y_continuous(breaks=c(-1.0,-0.5, 0, 0.5, 1.0)) +
  theme_hc()+ scale_colour_hc() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_fill_manual(values=c("#9a133d", "#b93961", "#d8527c", "#f28aaa", "#f9b4c9", "#f9e0e8", "#ffffff", "#eaf3ff", "#c5daf6", "#a1c2ed", "#6996e3", "#4060c8", "#1a318b", "#e7e5cc", "#c2d6a4", "#9cc184", "#669d62", "#3c7c3d", "#1f5b25","#1e3d14"))


ggarrange(F1, F2, F3 + rremove("x.text"), 
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3, 
          common.legend = TRUE, 
          legend = "right")




## -----------------------------------------------------------------------------

mean_load <- data.frame(matrix(ncol = 3, nrow = 0))
x <- c("ML1", "ML2", "ML3")
colnames(mean_load) <- x

columns <- list("ML1", "ML2", "ML3")

for (i in columns) {
  mean_load[1, i] <- mean(fa_loadings[1:3, i])
  mean_load[2, i] <- mean(fa_loadings[4:6, i])
  mean_load[3, i] <- mean(fa_loadings[7:9, i])
  mean_load[4, i] <- mean(fa_loadings[10:12, i])
  mean_load[5, i] <- mean(fa_loadings[13:15, i])
  mean_load[6, i] <- mean(fa_loadings[16:18, i])
  mean_load[7, i] <- mean(fa_loadings[19:21, i])
  mean_load[8, i] <- mean(fa_loadings[22:24, i])
  mean_load[9, i] <- mean(fa_loadings[25:27, i])
  mean_load[10, i] <- mean(fa_loadings[28:30, i])
  mean_load[11, i] <- mean(fa_loadings[31:33, i])
  mean_load[12, i] <- mean(fa_loadings[34:36, i])
  mean_load[13, i] <- mean(fa_loadings[37:39, i])
  mean_load[14, i] <- mean(fa_loadings[40:42, i])
  mean_load[15, i] <- mean(fa_loadings[43:45, i])
  mean_load[16, i] <- mean(fa_loadings[46:48, i])
  mean_load[17, i] <- mean(fa_loadings[49:51, i])
  mean_load[18, i] <- mean(fa_loadings[52:54, i])
  mean_load[19, i] <- mean(fa_loadings[55:57, i])
  mean_load[20, i] <- mean(fa_loadings[58:60, i])
}


numbers <- mean_load

mean_load$Subscale <- row.names(mean_load)
mean_load$Subscale[1] <- "Kindness"
mean_load$Subscale[2] <- "Teamwork"
mean_load$Subscale[3] <- "Creativity"
mean_load$Subscale[4] <- "Discretion"
mean_load$Subscale[5] <- "Humour"
mean_load$Subscale[6] <- "Forgiveness"
mean_load$Subscale[7] <- "Honesty"
mean_load$Subscale[8] <- "Love of Learning"
mean_load$Subscale[9] <- "Inspiration"
mean_load$Subscale[10] <- "Leadership"
mean_load$Subscale[11] <- "Perspective"
mean_load$Subscale[12] <- "Gratitude"
mean_load$Subscale[13] <- "Courage"
mean_load$Subscale[14] <- "Curiosity"
mean_load$Subscale[15] <- "Social Intelligence"
mean_load$Subscale[16] <- "Fairness"
mean_load$Subscale[17] <- "Hope"
mean_load$Subscale[18] <- "Enthusiasm"
mean_load$Subscale[19] <- "Love"
mean_load$Subscale[20] <- "Perseverance"


mean_load <- mean_load |>
  dplyr::rename("Factor 1" = "ML1",
                "Factor 2" = "ML2",
                "Factor 3" = "ML3")

mean_load <- mean_load[, c(4, 1, 2, 3)]

is.num <- sapply(mean_load, is.numeric)
mean_load[is.num] <- lapply(mean_load[is.num], round, 2)

BuYlRd <- function(x) rgb(colorRamp(c("#eaf3ff", "#89a6bb", "#454b87"))(x), maxColorValue = 255)


reactable(
  mean_load,
  defaultColDef = colDef(
    style = function(value, index, name) {
      if (is.numeric(value) && value >= 0.3) {
        normalized <- (value - min(numbers)) / (max(numbers) - min(numbers))
        color <- BuYlRd(value)
        list(fontWeight = "bold", background = color)
      } else if (is.numeric(value) && value >= 0.0) {
        normalized <- (value - min(numbers)) / (max(numbers) - min(numbers))
        color <- BuYlRd(value)
        list(background = color)
      }
    },
    format = colFormat(digits = 2),
    minWidth = 50,
    align = "center",
  ),
  columns = list(
    Subscale = colDef(style = list(borderRight = "1px solid rgba(0, 0, 0, 0.1)"))),
  borderless = TRUE,
  style = list(fontFamily = "Work Sans, sans-serif", fontSize = "0.875rem"),
  defaultPageSize = 20
)



## ----Wrangle datasets---------------------------------------------------------

total_factors_df <- merge(total_df, character_factors, by = "user_id")

total_factors_df$pronoun_1 <- as.factor(total_factors_df$pronoun_1)
total_factors_df$treatment_stage <- as.factor(total_factors_df$treatment_stage)

### Duration on platform

weekly_factors_df <- merge(time, total_factors_df, by = "user_id", all = TRUE) |>
                      rename(weekly_duration = time_on_platform_min)

# Check distributions
barplot(table(weekly_factors_df$weekly_duration), ylab = "Frequency", main = "weekly_duration")
barplot(table(weekly_factors_df$Factor1), ylab = "Frequency", main = "Factor1")
barplot(table(weekly_factors_df$Factor2), ylab = "Frequency", main = "Factor2")
barplot(table(weekly_factors_df$Factor3), ylab = "Frequency", main = "Factor3")
barplot(table(weekly_factors_df$age), ylab = "Frequency", main = "age")
barplot(table(weekly_factors_df$k10_total), ylab = "Frequency", main = "k10_total")
barplot(table(weekly_factors_df$pronoun_1), ylab = "Frequency", main = "pronoun_1")
barplot(table(weekly_factors_df$treatment_stage), ylab = "Frequency", main = "treatment_stage")


# Log transform numeric variables
weekly_factors_df$weekly_duration <- log(weekly_factors_df$weekly_duration + 1)
weekly_factors_df$k10_total <- log(weekly_factors_df$k10_total + 1)
weekly_factors_df$age <- log(weekly_factors_df$age + 1)



## ----Survival analysis--------------------------------------------------------

survival_df <- total_factors_df %>%  
  dplyr::filter(account_status != "INACTIVE") %>%
  dplyr::mutate(account_status = dplyr::recode(account_status, ACTIVE = 1, EXPIRED = 2)) %>%
  subset(, select = c("user_id", "active_days", "account_status", "Factor1", "Factor2", "Factor3", 
                      "age", "k10_total", "pronoun_1", "treatment_stage"))

survival <- na.omit(survival_df)

aft_weibull <-
  survreg(Surv(active_days, account_status) ~ Factor1 + Factor2 + Factor3 + age + k10_total + pronoun_1 + treatment_stage,
        data =  survival)
summary(aft_weibull)

aft_exp <-
  survreg(Surv(active_days, account_status) ~ Factor1 + Factor2 + Factor3 + age + k10_total + pronoun_1 + treatment_stage,
          data =  survival, dist = "exponential")
summary(aft_exp)

aft_lnorm <-
  survreg(Surv(active_days, account_status) ~ Factor1 + Factor2 + Factor3 + age + k10_total + pronoun_1 + treatment_stage,
          data =  survival, dist = "lognormal")
summary(aft_lnorm)

aft_llogis <-
  survreg(Surv(active_days, account_status) ~ Factor1 + Factor2 + Factor3 + age + k10_total + pronoun_1 + treatment_stage,
          data =  survival, dist = "loglogistic")
summary(aft_llogis)

unname(head(coef(aft_weibull)))
unname(head(coef(aft_exp)))
unname(head(coef(aft_lnorm)))
unname(head(coef(aft_llogis)))

AIC(aft_weibull,
    aft_exp,
    aft_lnorm,
    aft_llogis) # lnorm appears to be the best model


my_resid <- log(survival$active_days) - aft_lnorm$linear.predictors
qqnorm(my_resid)
qqline(my_resid)

# Log_norm model selected as best fit

## ----Longitudinal analysis----------------------------------------------------

# Check multicolinearity

vif(lm(weekly_duration ~ week_number + Factor1 + Factor2 + Factor3 + k10_total + age + treatment_stage + pronoun_1, data = weekly_factors_df))


time_mixmodel <-
  lmer(
    weekly_duration ~ week_number * Factor1 + week_number * Factor2 + week_number *
      Factor3 + (1 + week_number |
                               user_id),
    data = weekly_factors_df, 
    control = lmerControl(optimizer ="Nelder_Mead")
  )

summ(time_mixmodel)


# Check residual plots
res <- resid(time_mixmodel) 
plot(time_mixmodel) # Plot showing heteroscedacity
# QQ plot
qqnorm(res) 
qqline(res) 
# Density plot
plot(density(res))


## ----Longitudinal analysis with covariates - time on platform-----------------

time_mixmodel_cov <-
  lmer(
    weekly_duration ~ week_number * Factor1 + week_number * Factor2 + week_number *
      Factor3 + k10_total + age + treatment_stage + pronoun_1 + (1 + week_number |
                   user_id),
    data = weekly_factors_df, 
    control = lmerControl(optimizer ="Nelder_Mead")
  )

summary(time_mixmodel_cov)



# Check residual plots
res <- resid(time_mixmodel_cov) 
plot(fitted(time_mixmodel_cov), res) 
abline(0,0) 
qqnorm(res) 
qqline(res) 
plot(density(res))

cooks_distance <- cooks.distance(time_mixmodel_cov)
plot(cooks_distance)
abline(h = 4 / length(cooks_distance), col = "red", lty = 2)


# Compare model fits
AIC(time_mixmodel, time_mixmodel_cov)
BIC(time_mixmodel, time_mixmodel_cov)


fact1_long_time <-
  interact_plot(time_mixmodel_cov, pred = week_number, modx = Factor1,
              x.label = "Week", y.label = "Time on platform (log(minutes))", colors = "Greens",
              legend.main = "Social Harmony Factor") +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 15, angle = 90, margin=margin(0,20,0,0)), axis.title.x = element_text(size = 20, angle = 00, margin=margin(20,0,0,0))) +
  theme(plot.title = element_text(size = rel(3), angle = 00)) +
  theme(axis.text.x = element_text(angle = 00, size=20), axis.text.y = element_text(angle = 00, size=20)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
  theme(axis.line.x = element_line(color="black", linewidth = 1.2), axis.line.y = element_line(color="black", linewidth = 1.2)) +
  theme(axis.ticks.y=element_line(size=(1.5)), axis.ticks.x=element_line(size=(1.5)), axis.ticks.length=unit(0.4, "cm")) +
  scale_x_continuous(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8, 10, 12)) +
  scale_y_continuous(limits = c(0, 2)) 


fact2_long_time <-
  interact_plot(time_mixmodel_cov, pred = week_number, modx = Factor2,
              x.label = "Week", y.label = "Time on platform (log(minutes))", colors = "Reds",
              legend.main = "Positive Determination Factor") +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 15, angle = 90, margin=margin(0,20,0,0)), axis.title.x = element_text(size = 20, angle = 00, margin=margin(20,0,0,0))) +
  theme(plot.title = element_text(size = rel(3), angle = 00)) +
  theme(axis.text.x = element_text(angle = 00, size=20), axis.text.y = element_text(angle = 00, size=20)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
  theme(axis.line.x = element_line(color="black", linewidth = 1.2), axis.line.y = element_line(color="black", linewidth = 1.2)) +
  theme(axis.ticks.y=element_line(size=(1.5)), axis.ticks.x=element_line(size=(1.5)), axis.ticks.length=unit(0.4, "cm")) +
  scale_x_continuous(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8, 10, 12)) +
  scale_y_continuous(limits = c(0, 2)) 


fact3_long_time <-
  interact_plot(time_mixmodel_cov, pred = week_number, modx = Factor3,
              x.label = "Week", y.label = "Time on platform (log(minutes))", colors = "Purples", 
              legend.main = "Courage & Creativity Factor") +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 15, angle = 90, margin=margin(0,20,0,0)), axis.title.x = element_text(size = 20, angle = 00, margin=margin(20,0,0,0))) +
  theme(plot.title = element_text(size = rel(3), angle = 00)) +
  theme(axis.text.x = element_text(angle = 00, size=20), axis.text.y = element_text(angle = 00, size=20)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
  theme(axis.line.x = element_line(color="black", linewidth = 1.2), axis.line.y = element_line(color="black", linewidth = 1.2)) +
  theme(axis.ticks.y=element_line(size=(1.5)), axis.ticks.x=element_line(size=(1.5)), axis.ticks.length=unit(0.4, "cm")) +
  scale_x_continuous(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8, 10, 12)) +
  scale_y_continuous(limits = c(0, 2)) 

ggarrange(fact1_long_time,fact2_long_time,fact3_long_time,
                        labels = c("A", "B", "C"),
                        ncol = 1, nrow = 3, 
                        common.legend = FALSE)



## ----Longitudinal analysis - Number of sessions-------------------------------

sessions_long_df <- merge(session_number, total_factors_df, by = "user_id", all = TRUE) |>    
                      rename(weekly_sessions = session.Count.per.week)

sessions_long_df$weekly_sessions <- log(sessions_long_df$weekly_sessions + 1)
sessions_long_df$k10_total <- log(sessions_long_df$k10_total + 1)
sessions_long_df$age <- log(sessions_long_df$age + 1)

session_mixmodel <-
  lmer(
    weekly_sessions ~ week_number * Factor1 + week_number * Factor2 + week_number *
      Factor3 + k10_total + age + pronoun_1 + treatment_stage + (1 + week_number |
                                                 user_id),
    data = sessions_long_df, 
    control = lmerControl(optimizer ="Nelder_Mead")
  )

summary(session_mixmodel)

cooks_distance <- cooks.distance(session_mixmodel)
plot(cooks_distance)

# Check residual plots
res <- resid(session_mixmodel) 
plot(fitted(session_mixmodel), res) 
abline(0,0) 
qqnorm(res) 
qqline(res) 
plot(density(res))


fact1_long_session <-
  interact_plot(session_mixmodel, pred = week_number, modx = Factor1,
                x.label = "Week", y.label = "Sessions on platform (log)", colors = "Greens",
                legend.main = "Social Harmony Factor") +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 15, angle = 90, margin=margin(0,20,0,0)), axis.title.x = element_text(size = 20, angle = 00, margin=margin(20,0,0,0))) +
  theme(plot.title = element_text(size = rel(3), angle = 00)) +
  theme(axis.text.x = element_text(angle = 00, size=20), axis.text.y = element_text(angle = 00, size=20)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
  theme(axis.line.x = element_line(color="black", linewidth = 1.2), axis.line.y = element_line(color="black", linewidth = 1.2)) +
  theme(axis.ticks.y=element_line(size=(1.5)), axis.ticks.x=element_line(size=(1.5)), axis.ticks.length=unit(0.4, "cm")) +
  scale_x_continuous(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8, 10, 12)) +
  scale_y_continuous(limits = c(0,2))


fact2_long_session <-
  interact_plot(session_mixmodel, pred = week_number, modx = Factor2,
                x.label = "Week", y.label = "Sessions on platform (log)", colors = "Reds",
                legend.main = "Positive Determination Factor") +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 15, angle = 90, margin=margin(0,20,0,0)), axis.title.x = element_text(size = 20, angle = 00, margin=margin(20,0,0,0))) +
  theme(plot.title = element_text(size = rel(3), angle = 00)) +
  theme(axis.text.x = element_text(angle = 00, size=20), axis.text.y = element_text(angle = 00, size=20)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
  theme(axis.line.x = element_line(color="black", linewidth = 1.2), axis.line.y = element_line(color="black", linewidth = 1.2)) +
  theme(axis.ticks.y=element_line(size=(1.5)), axis.ticks.x=element_line(size=(1.5)), axis.ticks.length=unit(0.4, "cm")) +
  scale_x_continuous(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8, 10, 12)) +
  scale_y_continuous(limits = c(0,2))


fact3_long_session <-
  interact_plot(session_mixmodel, pred = week_number, modx = Factor3,
                x.label = "Week", y.label = "Sessions on platform (log)", colors = "Purples", 
                legend.main = "Courage & Creativity Factor") +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 15, angle = 90, margin=margin(0,20,0,0)), axis.title.x = element_text(size = 20, angle = 00, margin=margin(20,0,0,0))) +
  theme(plot.title = element_text(size = rel(3), angle = 00)) +
  theme(axis.text.x = element_text(angle = 00, size=20), axis.text.y = element_text(angle = 00, size=20)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
  theme(axis.line.x = element_line(color="black", linewidth = 1.2), axis.line.y = element_line(color="black", linewidth = 1.2)) +
  theme(axis.ticks.y=element_line(size=(1.5)), axis.ticks.x=element_line(size=(1.5)), axis.ticks.length=unit(0.4, "cm")) +
  scale_x_continuous(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8, 10, 12)) +
  scale_y_continuous(limits = c(0,2))

ggarrange(fact1_long_session,fact2_long_session,fact3_long_session,
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3, 
          common.legend = FALSE)


## ----Longitudinal analysis - Reactions on social network----------------------

reactions_long_df <- merge(reactions, total_factors_df, by = "user_id", all = TRUE) |>
                      rename(weekly_reactions = number.of.reactions.made.by.YP)

reactions_long_df$weekly_reactions <- log(reactions_long_df$weekly_reactions + 1)
reactions_long_df$k10_total <- log(reactions_long_df$k10_total + 1)
reactions_long_df$age <- log(reactions_long_df$age + 1)

reactions_mixmodel_cov <-
  lmer(
    weekly_reactions ~ week_number * Factor1 + week_number * Factor2 + week_number *
      Factor3 + k10_total + age + pronoun_1 + treatment_stage + (1 + week_number |
                                                 user_id),
    data = reactions_long_df, 
    control = lmerControl(optimizer ="Nelder_Mead")
  )

summary(reactions_mixmodel_cov)

cooks_distance <- cooks.distance(reactions_mixmodel_cov)
plot(cooks_distance)

# Check residual plots
res <- resid(reactions_mixmodel_cov) 
plot(fitted(reactions_mixmodel_cov), res) 
abline(0,0) 
qqnorm(res) 
qqline(res) 
plot(density(res))


fact1_long_reactions <-
  interact_plot(reactions_mixmodel_cov, pred = week_number, modx = Factor1,
                x.label = "Week", y.label = "Reactions on social network (log)", colors = "Greens",
                legend.main = "Social Harmony Factor") +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 15, angle = 90, margin=margin(0,20,0,0)), axis.title.x = element_text(size = 20, angle = 00, margin=margin(20,0,0,0))) +
  theme(plot.title = element_text(size = rel(3), angle = 00)) +
  theme(axis.text.x = element_text(angle = 00, size=20), axis.text.y = element_text(angle = 00, size=20)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
  theme(axis.line.x = element_line(color="black", linewidth = 1.2), axis.line.y = element_line(color="black", linewidth = 1.2)) +
  theme(axis.ticks.y=element_line(size=(1.5)), axis.ticks.x=element_line(size=(1.5)), axis.ticks.length=unit(0.4, "cm")) +
  scale_x_continuous(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8, 10, 12)) +
  scale_y_continuous(limits = c(0,2))


fact2_long_reactions <-
  interact_plot(reactions_mixmodel_cov, pred = week_number, modx = Factor2,
                x.label = "Week", y.label = "Reactions on social network (log)", colors = "Reds",
                legend.main = "Positive Determination Factor") +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 15, angle = 90, margin=margin(0,20,0,0)), axis.title.x = element_text(size = 20, angle = 00, margin=margin(20,0,0,0))) +
  theme(plot.title = element_text(size = rel(3), angle = 00)) +
  theme(axis.text.x = element_text(angle = 00, size=20), axis.text.y = element_text(angle = 00, size=20)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
  theme(axis.line.x = element_line(color="black", linewidth = 1.2), axis.line.y = element_line(color="black", linewidth = 1.2)) +
  theme(axis.ticks.y=element_line(size=(1.5)), axis.ticks.x=element_line(size=(1.5)), axis.ticks.length=unit(0.4, "cm")) +
  scale_x_continuous(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8, 10, 12)) +
  scale_y_continuous(limits = c(0,2))


fact3_long_reactions <-
  interact_plot(reactions_mixmodel_cov, pred = week_number, modx = Factor3,
                x.label = "Week", y.label = "Reactions on social network (log)", colors = "Purples", 
                legend.main = "Courage & Creativity Factor") +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 15, angle = 90, margin=margin(0,20,0,0)), axis.title.x = element_text(size = 20, angle = 00, margin=margin(20,0,0,0))) +
  theme(plot.title = element_text(size = rel(3), angle = 00)) +
  theme(axis.text.x = element_text(angle = 00, size=20), axis.text.y = element_text(angle = 00, size=20)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
  theme(axis.line.x = element_line(color="black", linewidth = 1.2), axis.line.y = element_line(color="black", linewidth = 1.2)) +
  theme(axis.ticks.y=element_line(size=(1.5)), axis.ticks.x=element_line(size=(1.5)), axis.ticks.length=unit(0.4, "cm")) +
  scale_x_continuous(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8, 10, 12)) +
  scale_y_continuous(limits = c(0,2))



## ----Longitudinal analysis - Interactions with clinical team------------------

message_long_df <- merge(message_number, total_factors_df, by = "user_id", all = TRUE) |>
                      rename(weekly_messages = message_number.x)

message_long_df$weekly_messages <- log(message_long_df$weekly_messages + 1)
message_long_df$k10_total <- log(message_long_df$k10_total + 1)
message_long_df$age <- log(message_long_df$age + 1)

message_mixmodel_cov <-
  lmer(
    weekly_messages ~ week_number * Factor1 + week_number * Factor2 + week_number *
      Factor3 + k10_total + age + pronoun_1 + treatment_stage + (1 + week_number |
                                                                   user_id),
    data = message_long_df, 
    control = lmerControl(optimizer ="Nelder_Mead")
  )

summary(message_mixmodel_cov)

cooks_distance <- cooks.distance(message_mixmodel_cov)
plot(cooks_distance)

# Check residual plots
res <- resid(message_mixmodel_cov) 
plot(fitted(message_mixmodel_cov), res) 
abline(0,0) 
qqnorm(res) 
qqline(res) 
plot(density(res))


fact1_long_message <-
  interact_plot(message_mixmodel_cov, pred = week_number, modx = Factor1,
                x.label = "Week", y.label = "Messages to clinical team (log)", colors = "Greens",
                legend.main = "Social Harmony Factor") +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 15, angle = 90, margin=margin(0,20,0,0)), axis.title.x = element_text(size = 20, angle = 00, margin=margin(20,0,0,0))) +
  theme(plot.title = element_text(size = rel(3), angle = 00)) +
  theme(axis.text.x = element_text(angle = 00, size=20), axis.text.y = element_text(angle = 00, size=20)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
  theme(axis.line.x = element_line(color="black", linewidth = 1.2), axis.line.y = element_line(color="black", linewidth = 1.2)) +
  theme(axis.ticks.y=element_line(size=(1.5)), axis.ticks.x=element_line(size=(1.5)), axis.ticks.length=unit(0.4, "cm")) +
  scale_x_continuous(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8, 10, 12)) +
  scale_y_continuous(limits = c(0,2))


fact2_long_message <-
  interact_plot(message_mixmodel_cov, pred = week_number, modx = Factor2,
                x.label = "Week", y.label = "Messages to clinical team (log)", colors = "Reds",
                legend.main = "Positive Determination Factor") +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 15, angle = 90, margin=margin(0,20,0,0)), axis.title.x = element_text(size = 20, angle = 00, margin=margin(20,0,0,0))) +
  theme(plot.title = element_text(size = rel(3), angle = 00)) +
  theme(axis.text.x = element_text(angle = 00, size=20), axis.text.y = element_text(angle = 00, size=20)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
  theme(axis.line.x = element_line(color="black", linewidth = 1.2), axis.line.y = element_line(color="black", linewidth = 1.2)) +
  theme(axis.ticks.y=element_line(size=(1.5)), axis.ticks.x=element_line(size=(1.5)), axis.ticks.length=unit(0.4, "cm")) +
  scale_x_continuous(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8, 10, 12)) +
  scale_y_continuous(limits = c(0,2))


fact3_long_message <-
  interact_plot(message_mixmodel_cov, pred = week_number, modx = Factor3,
                x.label = "Week", y.label = "Messages to clinical team (log)", colors = "Purples", 
                legend.main = "Courage & Creativity Factor") +
  theme_classic() + 
  theme(axis.title.y = element_text(size = 15, angle = 90, margin=margin(0,20,0,0)), axis.title.x = element_text(size = 20, angle = 00, margin=margin(20,0,0,0))) +
  theme(plot.title = element_text(size = rel(3), angle = 00)) +
  theme(axis.text.x = element_text(angle = 00, size=20), axis.text.y = element_text(angle = 00, size=20)) +
  theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15)) +
  theme(axis.line.x = element_line(color="black", linewidth = 1.2), axis.line.y = element_line(color="black", linewidth = 1.2)) +
  theme(axis.ticks.y=element_line(size=(1.5)), axis.ticks.x=element_line(size=(1.5)), axis.ticks.length=unit(0.4, "cm")) +
  scale_x_continuous(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8, 10, 12)) +
  scale_y_continuous(limits = c(0,2))

ggarrange(fact1_long_message,fact2_long_message,fact3_long_message,
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3, 
          common.legend = FALSE)

# Single graph

ggarrange(fact1_long_time,fact2_long_time,fact3_long_time,
          fact1_long_session,fact2_long_session,fact3_long_session,
          fact1_long_reactions,fact2_long_reactions,fact3_long_reactions,
          fact1_long_message,fact2_long_message,fact3_long_message,
          labels = c("A", "B", "C",
                     "D", "E", "F",
                     "G", "H", "I",
                     "J", "K", "L"),
          ncol = 3, nrow = 4, 
          common.legend = TRUE, 
          legend = "right")



