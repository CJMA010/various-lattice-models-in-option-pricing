setwd("~/Documents/projects/c++_project/opt_pricing_binomial_american_put/DerivedData/opt_pricing_binomial_american_put/Build/Products/Debug")
library(ggplot2)


result_a <- read.csv('results_a.csv' , colClasses = c(rep("double", 4), "NULL"),header = FALSE)

result_a$V1 <- -log10(result_a$V1)

ggplot( result_a, aes(V1) ) + 
  geom_line( aes( y = V2 , col='CRR'), size = 0.25) +
  geom_line( aes( y = V3 , col='Jarrow-Rudd' ), size = 0.25 ) +
  geom_line( aes( y = V4 , col='Tian' ), size = 0.25 ) +
  labs( title = "american put option price")+
  xlab('delta t')+
  ylab('price')+
  scale_x_continuous(breaks=-log10(c(1,0.1,0.01,0.002)),labels=c("1","0.1" ,"0.01", "0.002"))





result_b <- read.csv('results_b.csv' , colClasses = c(rep("double", 5), "NULL"),header = FALSE)

ggplot( result_b, aes(V1) ) + 
  geom_line( aes( y = V2 , col='Simple Binomial'), size = 0.25,alpha=0.5) +
  geom_line( aes( y = V3 , col='Binomial BlackScholes'), size = 0.25,alpha=0.5) +
  geom_line( aes( y = V4 , col='BBS with Richardson Extrapolation' ), size = 0.25,alpha=0.5 ) +
  geom_line( aes( y = V5 , col='Binomial Average Method' ), size = 0.25,alpha=0.5 ) +
  labs( title = "american put option price")+
  xlab('delta t')+
  ylab('price')

result_c <- read.csv('results_c.csv' , colClasses = c("double","double" , "NULL"),header = FALSE)

ggplot( result_c, aes(V1) ) + 
  geom_line( aes( y = V2 , col='Trinomial'), size = 0.25,alpha=0.5) +
  labs( title = "american put option price")+
  xlab('delta t')+
  ylab('price')
  #scale_x_continuous(breaks=-log10(c(1,0.1,0.01,0.002)),labels=c("1","0.1" ,"0.01", "0.002"))



comapre <- data.frame(cbind(result_b$V1,result_b$V4,result_c$V2))

ggplot( comapre, aes(X1) ) + 
  geom_line( aes( y = X2 , col='BBSR'), size = 0.25,alpha=0.5) +
  geom_line( aes( y = X3 , col='Trinomial'), size = 0.25,alpha=0.5) +
  labs( title = "american put option price")+
  xlab('delta t')+
  ylab('price')
#scale_x_continuous(breaks=-log10(c(1,0.1,0.01,0.002)),labels=c("1","0.1" ,"0.01", "0.002"))

