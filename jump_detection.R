library(ISwR)
library(splines)
library(data.table)
library(grid)
library(gridExtra)
library(contrast)
library(tidyverse)
library(dplyr)
library(ggplot2)

ep_grp <- matrix(c("ALCOPANCACU",       "K11",
                   "CARMYOPATHALCO",    "I9",
                   "EPISEIZALCO",       "G6",
                   "GEST",              "O15",
                   "HEPATITIS",         "AB1",
                   "PNEUMO",            "AB1",
                   "RESP",              "J10",
                   "RHEUMA",            "M13",
                   "TRAUMBRAIN",        "OTHER",
                   "TUBERCULOSIS",      "AB1", 
                   "ALCOGASTRITIS",     "K11",  
                   "ALCOHOLMENTAL",     "F5", 
                   "ALCOLIVER",         "K11",
                   "ALCONEURODEGEN",    "G6",
                   "ALCOPANCCHRON",     "K11",
                   "ALCOPOLYNEU",       "G6",
                   "CHARCOT",           "OTHER",
                   "CHILDHOOD",         "OTHER",
                   "CHIRHEP",           "K11",
                   "COPD",              "J10",
                   "DENTAL",            "OTHER",
                   "DRY",               "H8",
                   "FTD",               "OTHER",
                   "GIANT",             "OTHER",
                   "ILD",               "J10",
                   "INFLUENZA",         "J10",
                   "LUNG",              "CD3",
                   "MACULA",            "H8",
                   "PANCREATITIS",      "K11",
                   "PDSTRICT",          "G6",
                   "PNEUMOBACTKNOWN",   "J10",
                   "PSOR",              "M13",
                   "SPONDYLOARTHRITIS", "M13",
                   "ST19",              "OTHER",
                   "T2D",               "E4",
                   "THYROIDITIS",       "E4",
                   "WET",               "H8",
                   "SCND",              "M13",
                   "ASTHMA",            "J10",
                   "MACULAR",           "H8",
                   "O15",               "O15",
                   "KRA",               "F5",
                   "AB1",               "AB1",
                   "DM",                "E4",
                   "J10",               "J10",
                   "N14",               "N14",
                   "C3",                "CD2",
                   "D3",                "D3",
                   "H8",                "H8",
                   "I9",                "I9",
                   "K11",               "K11",
                   "F5",                "F5",
                   "G6",                "G6",
                   "CD2",               "CD2",
                   "E4",                "E4",
                   "L12",               "L12",
                   "H7",                "H7",
                   "M13",               "M13"), byrow=T, ncol=2)
colnames(ep_grp) <- c("group","grp_icd10")

ep_grp_text <- matrix(c("OTHER", "Other",
                        "O15",   "Pregnancy, childbirth & puerperium",
                        "AB1",   "Infectious & parasitic",
                        "J10",   "Respiratory system",
                        "N14",   "Genitourinary system",
                        "D3",    "Blood & immune mechanism",
                        "H8",    "Ear & mastoid process",
                        "I9",    "Circulatory system",
                        "K11",   "Digestive system",
                        "F5",    "Mental & behavioural",
                        "G6",    "Nervous system",
                        "CD2",   "Neoplasms",
                        "E4",    "Endocrine, nutritional & metabolic",
                        "L12",   "Skin & subcutaneous tissue",
                        "H7",    "Eye & adnexa",
                        "M13",   "Musculoskeletal system"), byrow=T, ncol=2)
colnames(ep_grp_text) <- c("grp_icd10","grp_text")


# Load data from share folder
df2 <- read.csv("share/complete_endpoints.csv")
alive <- read.csv("share/index_birth_death.csv")

# Narrow down disease scope
list_disease_1 <- df2 %>% 
  mutate(event_yr=as.numeric(event_yr)) %>% 
  filter((event_yr <= 1997) & (event_yr >= 1995)) %>% 
  group_by(event) %>% 
  count() %>% 
  filter(n>5) %>% 
  select(event)
list_disease_2 <- df2 %>% 
  mutate(event_yr=as.numeric(event_yr)) %>% 
  filter((event_yr <= 2001) & (event_yr >= 1999)) %>% 
  group_by(event) %>% 
  count() %>% 
  filter(n>5) %>% 
  select(event)
list_disease <- intersect(list_disease_1$event,list_disease_2$event)


group_id <- sapply(strsplit(list_disease, '_'), '[[', 1)
group_id <- data.frame(group_id, list_disease)
colnames(group_id) <- c('group', 'event')
ep <- left_join(data.frame(ep_grp),data.frame(ep_grp_text))
groups <- select(left_join(group_id, ep),'event','grp_text')

df2 <- df2[complete.cases(df2), ] 
df2 <- select(df2,"ID","event","event_yr","birth_yr","age")
bins = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 120)
df2 <- df2 %>% filter(as.numeric(event_yr)>1985) %>%
  filter(as.numeric(event_yr)<2012) %>%
  mutate(age_2020=2020-birth_yr,
         age_2020_cut=cut(age_2020, breaks=bins, right=FALSE),
         age_cut=cut(age,breaks=bins,right=FALSE),
         event_yr=as.numeric(event_yr))

neventdf <- alive %>%  
  mutate(age_2020=2020-birth_yr,
         age_2020_cut=cut(age_2020, breaks=bins,right=FALSE))  

# Grouped by year bins
neventdf_grouped <- neventdf %>% 
  group_by(death_yr,age_2020_cut) %>% 
  count() %>% 
  filter(death_yr!=9999) %>% 
  mutate(alive_n=nrow(neventdf)-n) %>% 
  select(-c(n))

event_grouped2 <- df2 %>% group_by(event,event_yr,age_2020_cut) %>% count()
event_grouped2$index <- paste0(as.character(event_grouped2$event_yr),as.character(event_grouped2$age_2020_cut))

endpoints = groups %>% select(event)
endpoint_list <- endpoints$event

RES <- NULL

for (dis in endpoint_list)
{
  # Model
  ## Poisson Regression Model: Looking for jump in the data
  
  eventdf_grouped2 <- event_grouped2 %>% filter(event==dis) 
  
  grouped2 <- inner_join(eventdf_grouped2,neventdf_grouped,by = c("event_yr" = "death_yr", "age_2020_cut" = "age_2020_cut"))
  grouped2$event_yr_cut <- factor(ifelse((grouped2$event_yr <= 1994) | (grouped2$event_yr >= 2002),0,ifelse(grouped2$event_yr <= 1997,1,ifelse(grouped2$event_yr >= 1999,2,0))))
  
  fit1 <- glm(n ~ event_yr + event_yr_cut + age_2020_cut +offset(log(alive_n)), family=poisson, data=grouped2,na.action="na.exclude",contrasts = list(event_yr_cut = contr.treatment(n = 3, base = 2)))
  
  # Save coefficients
  coef <- summary(fit1)$coef[4,1]
  se <- summary(fit1)$coef[4,2]
  zscore <- summary(fit1)$coef[4,3]
  pval <- summary(fit1)$coef[4,4]
  
  # Save plot
  grouped2_ <- grouped2 %>% group_by(event_yr) %>% summarize(n=sum(n),pop=sum(alive_n))
  ggplot(data = grouped2_, aes(x = event_yr, y = n/pop)) + 
    geom_point() + xlab("Year") + ylab("Prevalence") + theme_bw() + geom_line()
  image <- ggplot(data = grouped2_, aes(x = event_yr, y = n/pop)) + 
    geom_rect(aes(xmin=1994.5, xmax=1997.5, ymin=min(n/pop), ymax=max(n/pop)), fill="lightgrey", alpha=0.2) + 
    geom_rect(aes(xmin=1998.5, xmax=2001.5, ymin=min(n/pop), ymax=max(n/pop)), fill="lightgrey", alpha=0.2) + 
    geom_point() + xlab("Year") + ylab("Prevalence") + theme_bw() + geom_line()
  ggsave(paste0('plots/',groups[groups$event == dis,]$grp_text,'-',dis,'.png'), image, width = 8, height = 8, dpi = 100)
  
  # Save results
  RES <- rbind(RES,c(groups[groups$event == dis,]$grp_text,dis,
                     coef,se,pval,zscore,
                     grouped2_$event_yr[1],grouped2_$event_yr[length(grouped2_$event_yr)],mean(grouped2_$n)))
}

RESS <- data.frame(RES)
colnames(RESS) <- c('group',"endpoint",
                    "beta","se","pval","zscore","start","end","avgn")
write.csv(RESS,file="results_short_window_poisson.csv", row.names = FALSE)

