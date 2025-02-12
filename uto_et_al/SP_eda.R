## Exploratory data analysis 
#########################
library(ggplot2)

dat_org <- read.table("data/OSCE_data.csv", header=TRUE, sep=",")

str(dat_org)
head(dat_org)

## code from uto et al. to reorganize the data 

setting <- list(K = 4, n_person = 30, n_item = 30, n_rater = 5)
dat <- c()
for(i in 1:length(dat_org[,1])){
  for(j in 3:length(dat_org[1,])){
    if(!is.na(dat_org[i, j])){
      dat <- rbind(dat, cbind(dat_org[i,1:2], j-2, dat_org[i, j]+1))
    }
  }
}
colnames(dat) <- c("Items","Raters","Examinees","Score")

str(dat)


## plot the data  
ggplot(dat, aes(x=Score)) + 
  geom_histogram(aes(y = after_stat(density)), binwidth=1, fill="lightblue", color="black") + 
  facet_grid(Raters ~ Items) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## plot distribution across items of raters score

ggplot(dat, aes(x=Score)) + 
  geom_histogram(aes(y = after_stat(density)), binwidth=1, fill="lightblue", color="black") + 
  facet_grid(~ Raters) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(raters)
install.packages("immer")
