#### Header ####

##Desc: Instructions for building the Rényi divergence rates database and plotting two-state Markov chains divergence rates.
##Date: 2019/12/16
##Author: Philippe Regnault (philippe.regnault(at)univ-reims.fr)

####Preamble####

#Defining transition matrices to be considered for the second transition matrix
 Q1 <- matrix(0.5, nrow = 2, ncol =2) #Uniform matrix
 Q2 <- matrix(c(0.2, 0.8, 0.7, 0.3), nrow = 2, byrow = TRUE)
 Q3 <- matrix(c(0.2, 0.8, 0.3, 0.7), nrow = 2, byrow = TRUE)

#KL divergence
kldiv <- function(p, q){
  return(sum(p*log(p/q)))
}

#stationary distribution
sta_dist <- function(p, q){
  res <- eigen(matrix(c(1-p, p, q, 1-q), nrow = 2, byrow = TRUE))$vectors[,1]
  res <- res/sum(res)
  return(res[1])
}


 
#Mixture matrix function
mixture_matrix <- function(P, Q, s){
 return(P^s*Q^(1-s))
}
 
 
 
####Creating Rényi divergence database####

#Discretizing the (0,1]*(0,1] square (non-diagonal coefficients of the transition matrix of first MC)
p <- round(seq(0.01, 0.99, by = 0.01),2)
q <- round(seq(0.01, 0.99, by = 0.01),2)
#A set of values for s
#s <- c(round(seq(0.1, 2, by = 0.1),1), 3:10, 20, 100)
s <- c(0.5, 1, 2)
#A set of values for the transition matric of thr second MC
#3 possiblities
Q <- c('Q1', 'Q2', 'Q3')
#Create a database with all possible 4-tuples of p, q, s and Q
expand.grid(p = p, q = q, s = s, Q = Q, stringsAsFactors = FALSE) -> renyi_div_database
#Processing computation of divergence rates
library('dplyr')
renyi_div_database %>%
  rowwise() %>%
  mutate(lambda = max(Re(eigen(mixture_matrix(matrix(c(1-p, p, q, 1-q), nrow = 2, byrow = TRUE), get(Q), s))$values)),
         sta_dist = sta_dist(p,q),
         div_rate = case_when(s != 1 ~log(lambda)/(s-1), 
                              TRUE ~ sum(c(sta_dist, 1-sta_dist)* apply(matrix(c(1-p, p, q, 1-q), nrow = 2, byrow = TRUE)*log (matrix(c(1-p, p, q, 1-q), nrow = 2, byrow = TRUE) / get(Q)), 1, sum )))
         ) -> renyi_div_database
#Checking database
head(renyi_div_database)
str(renyi_div_database)


####Plotting divergence rates####
library('ggplot2')
renyi_div_database %>% 
  ggplot(mapping = aes(x = p, y = q, z = div_rate)) +
  facet_grid(Q~s) +
  geom_contour(aes(colour = ..level..), bins = 20) +
  scale_color_gradient(low = "blue", high = "red") + 
  coord_fixed() +
  theme_bw() -> plot_div

setwd("/home/philippe/Documents/Conferences/20190706_IEEEISIT")
png(file = "renyi_div_2state.png", width = 960, height = 960)
par(cex = 3)
plot_div + theme(axis.text.y=element_text(size=20),
                 axis.text.x=element_text(size=20, angle = 90),
                 axis.title=element_text(size=20,face="bold"),
                 legend.text = element_text(size = 14),
                 legend.title = element_text(size = 14),
                 strip.text = element_text(size = 14),
                 strip.text.y = element_text(angle = 90))
dev.off()

#### The same, but representing the divergence as a function of the second stochastic matrix ####
P1 <- Q1
P2 <- Q2
P3 <- Q3
P <- c('P1', 'P2', 'P3')
#Create a database with all possible 4-tuples of p, q, s and Q
expand.grid(p = p, q = q, s = s, P = P, stringsAsFactors = FALSE) -> renyi_div_database

library('dplyr')
renyi_div_database %>%
  rowwise() %>%
  mutate(lambda = max(Re(eigen(mixture_matrix(matrix(c(1-p, p, q, 1-q), nrow = 2, byrow = TRUE), get(P),s))$values)),
         sta_dist = sta_dist(get(P)[1,2],get(P)[2,1]),
         div_rate = case_when(s != 1 ~log(lambda)/(s-1), 
                              TRUE ~ sum(c(sta_dist, 1-sta_dist)* apply(get(P)*log(get(P) / matrix(c(1-p, p, q, 1-q), nrow = 2, byrow = TRUE)), 1, sum )))
  ) -> renyi_div_database
#Checking database
head(renyi_div_database)
str(renyi_div_database)


library('ggplot2')
renyi_div_database %>% 
  ggplot(mapping = aes(x = p, y = q, z = div_rate)) +
  facet_grid(P~s,
             labeller = label_both) +
  geom_contour(aes(colour = ..level..), bins = 20) +
  scale_color_gradient(low = "blue", high = "red") + 
  coord_fixed() +
  theme_bw() -> plot_div


setwd("/home/philippe/Seafile/Conferences/20200621_IEEEISIT/Conference_paper/")
# png(file = "renyi_div_2state_ISIT2020.png", width = 960, height = 960)
pdf(file = "renyi_div_2state_ISIT2020.pdf", width = 10, height = 10)
par(cex = 3)
plot_div + 
  theme(axis.text.y=element_text(size=20),
                 axis.text.x=element_text(size=20, angle = 90),
                 axis.title=element_text(size=20,face="bold"),
                 legend.text = element_text(size = 14),
                 legend.title = element_text(size = 14),
                 strip.text = element_text(size = 14),
                 strip.text.y = element_text(angle = 90)) + 
  scale_x_continuous(name = expression(q[12])) + 
  scale_y_continuous(name = expression(q[21]))
dev.off()

