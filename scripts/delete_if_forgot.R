Simuls_bfix_beta$com <- as.numeric(Simuls_bfix_beta$cj<4)
alpha_com <- with(Simuls_bfix_beta,aggregate(com,by=list(sub,condition),mean))
names(alpha_com) <- c("sub","condition","com")
alpha_com <- cast(alpha_com,sub~condition)
alpha_com$diff <- alpha_com$plus - alpha_com$minus

sim_diff_com <- beta_com$diff - alpha_com$diff
cor.test(sim_diff_com,com_bic$bic_diff)

m.int.a <- glmer(data=subset(Data,manip=="alpha"),com~condition*difflevel + (1|sub),family='binomial')
m.cond.a <- glmer(data=subset(Data,manip=="alpha"),com~condition*difflevel + (condition|sub),family='binomial',
                  control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=2e5)))
anova(m.int.a,m.cond.a)
Anova(m.cond.a)

m.int.b <- glmer(data=subset(Data,manip=="beta"),com~condition*difflevel + (1|sub),family='binomial')
m.cond.b <- glmer(data=subset(Data,manip=="beta"),com~condition*difflevel + (condition|sub),family='binomial',
                  control = glmerControl(optimizer = "bobyqa",optCtrl=list(maxfun=2e5)))
anova(m.int.b,m.cond.b)
Anova(m.cond.b)
emm <- emmeans(m.cond.b, ~ condition)
pairs(emm)
