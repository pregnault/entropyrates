#### Header ####

##Desc: Instructions for building the entropy rate database and plotting two-state Markov chain entropy rates and related quantities.
##Date: 2018/08/29
##Author: Philippe Regnault (philippe.regnault(at)univ-reims.fr)

#### Preambule ####

##Definition of functions for computing all involved quantities from the two non-diagonal coefficients of the transition matrix, the parameter s, etc.
#Discriminant of the characteristic polynomial of the perturbed transition matrix
delta <- function(p,q,s){
  res <- ( (1-p)^s - (1-q)^s )^2 + 4*p^s*q^s
  return(res)
}
#Its derivative with respect to s
deltap <- function(p,q,s){
  res <- 2*( log(1-p*(p!=1)) * (1-p)^s  - log(1-q*(q!=1)) * (1-q)^s ) * ( (1-p)^s - (1-q)^s) + 4*log(p*q)*p^s*q^s  
  return(res)
}
#The dominant eigenvalue of the perturbed transition matrix
lambda <- function(p,q,s, delta){
  res <- ( (1-p)^s + (1-q)^s + sqrt(delta)) /2
  return(res)
}
#Its derivative with respect to s
lambdap <- function(p,q,s, delta, deltap) {
  res <- ( log(1-p*(p!=1))*(1-p)^s + log(1-q*(q!=1))*(1-q)^s + deltap / (2*sqrt(delta))) / 2
  return(res)
}
#
eta <- function(p,q,s,delta){
  res <- (sqrt(delta) + (1-q)^s - (1-p)^s) /2
  return(res)
}
#Normalized left eigenvector coordinates
u1 <- function(q,s,eta) {
  return(q^s / (q^s + eta) )
}
u2 <- function(q,s,eta) {
  return(eta / (q^s + eta) )
}
#Normalized right eigenvector coordinates
v1 <- function(p,q,s,eta){
  res <- p^s * (q^s + eta) / (p^s*q^s + eta^2)
}
v2 <- function(p,q,s,eta){
  res <- eta * (q^s + eta) / (p^s*q^s + eta^2)
}
#The residual initial information term 
cs <- function(init1, init2, v1, v2){
  return(init1*v1 + init2*v2)
}

#### Creating entropy rate database ####

#Discretizing the (0,1]*(0,1] square (non-diagonal coefficients of the transition matrix)
p <- round(seq(0.01, 1, by = 0.01),2)
q <- round(seq(0.01, 1, by = 0.01),2)
#A set of values for s
s <- c(round(seq(0.1, 2, by = 0.1),1), 3:10, 20, 100)
#A set of values for the initial distribution
init1 <- round(seq(0, 1, by = 0.1),1)
#Create a database with all possible 4-tuples of p, q, s and init1
expand.grid(p = p, q = q, s = s, init1 = init1) -> ent_database
#Processing computation of spectral elements and entropy rates with the dplyr package
library('dplyr')
#Processing may take a few seconds
#proc.time() -> t
ent_database %>%
  mutate(init2 = 1-init1,
         sta1 = q/(p+q),
         sta2 = p/(p+q),
         delta = delta(p,q,s),
         deltap = deltap(p,q,s),
         lambda = lambda(p,q,s,delta),
         lambdap = lambdap(p,q,s,delta, deltap),
         eta = eta(p,q,s,delta),
         u1 = u1(q,s,eta),
         u2 = u2(q,s,eta),
         v1 = v1(p,q,s,eta),
         v2 = v2(p,q,s,eta),
         c = cs(init1,init2,v1,v2),
         # shannon1 = case_when(s == 1 ~ -lambdap),
         shannon = -sta1 * (p*log(p) + (1-p)*log(1-p*(p!=1))) - sta2*(q*log(q) + (1-q)*log(1-q*(q!=1))) ,
         shannon_type = 2,
         renyi = case_when(s != 1 ~ log(lambda) / (1-s), s == 1 ~ shannon),
         renyi_type = 2,
         tsallis = case_when(s >1 ~ 1/(s-1), s < 1 ~ c / (1-s)),
         tsallis_type = case_when(s >1 ~ 0, s > 1 ~ 3),
         taneja = -c*lambdap/lambda,
         taneja_type = case_when(s >1 ~ -1, s > 1 ~ 3)
  ) -> ent_database
#proc.time()-t

#uncomment the following chunk to write the database in a csv file.
# write.table(x = ent_database, file = "Entropy_rate_database.csv", sep = "\t", dec = ".", 
#             row.names = FALSE, col.names = TRUE,
#             fileEncoding = "UTF-8")


#### Graphical representations of lambda(p,q,s) ####

#As a function of p,q, for different values for s
#Using the ggplot2 package
library('ggplot2')
ent_database %>%
  filter(s %in% c(round(seq(0.2, 2, by = 0.2),2) , 3, 5, 10) & s != 1, init1 == 0.5) %>%
  ggplot(mapping = aes(x = p, y = q, z = lambda)) +
  facet_wrap(~s) +
  geom_contour(aes(colour = ..level..), bins = 20) +
  theme_bw() -> plot_lambdapq
#Uncomment the following chunk to export the figure as a ps file
# postscript(file = 'lambdapq.ps')
# par(cex = 1.4)
# plot_lambdapq + coord_fixed()
# dev.off()

#As a function of s, for different values for p and q
ent_database %>%
  filter(p %in% c(0.3, 0.5, 1) & q %in% c(0.3, 0.5, 1), init1 == 0.5, s!= 100) %>%
  ggplot(mapping = aes(x = s, y = lambda)) +
  facet_grid(q~p) +
  geom_path() +
  theme_bw() -> plot_lambdas
#Uncomment the following chunk to export the figure as a ps file
# postscript(file = 'lambdas.ps')
# par(cex = 1.4)
# plot_lambdas + coord_fixed()
# dev.off()

#### Graphical representations of lambdap ####

#As a function of p,q, for different values for s
library('ggplot2')
ent_database %>%
  filter(s %in% c(round(seq(0.2, 2, by = 0.2),2) , 3, 5, 10) & s != 1, init1 == 0.5) %>%
  ggplot(mapping = aes(x = p, y = q, z = lambdap)) +
  facet_wrap(~s) +
  geom_contour(aes(colour = ..level..), bins = 20) +
  theme_bw() -> plot_lambdappq
#Uncomment the following chunk to export the figure as a ps file
# postscript(file = 'lambdappq.ps')
# par(cex = 1.4)
# plot_lambdappq + coord_fixed()
# dev.off()

#As a function of s, for different values for p and q
ent_database %>%
  filter(p %in% c(0.3, 0.5, 1) & q %in% c(0.3, 0.5, 1), init1 == 0.5, s!= 100) %>%
  ggplot(mapping = aes(x = s, y = lambdap)) +
  facet_grid(q~p) +
  geom_path() +
  theme_bw() -> plot_lambdaps
#Uncomment the following chunk to export the figure as a ps file
# postscript(file = 'lambdaps.ps')
# par(cex = 1.4)
# plot_lambdaps + coord_fixed()
# dev.off()


# #### Graphical representations of lambda, lambdap, c en fonction de s
# ent_database %>%
#   select(p,q,s, init1, lambda, lambdap, c) %>%
#   filter(p == 0.2, q== 0.3, init1 == 0.5) %>%  
#   gather(key = func, value = value, -(1:4)) %>%
#   ggplot(mapping = aes(x = s, y= value, group = func)) +
#   geom_path(mapping = aes(colour = func, linetype = func)) +
#   theme_bw()
# 
# ent_database %>%
#   select(p,q,s, init1, lambda, lambdap, c) %>%
#   filter(p==0.7, q==0.3, init1 == 0.5) %>%  
#   gather(key = func, value = value, -(1:4)) %>%
#   ggplot(mapping = aes(x = s, y= value, group = func)) +
#   geom_path(mapping = aes(colour = func, linetype = func)) +
#   theme_bw() -> plt_functions_2
# plt_functions_2 + theme(axis.text = element_text(size = 12), 
#     axis.text.x = element_text(size = 12), 
#     axis.text.y = element_text(size = 12), 
#     legend.text = element_text(size = 10), 
#     legend.title = element_text(face = "bold"))

#### Graphical representaiton of u, v ####

##Left eigenvector U
#As functions of p,q
ent_database %>%
  select(p,q,s, init1, u1, u2) %>%
  filter(init1 == 0.5, s %in% c(0.5, 2, 5)) %>%  
  gather(key = Eigenvectors, value = value, -(1:4)) %>%
  ggplot(mapping = aes(x = p, y= q, z = value)) +
  facet_grid(Eigenvectors~s) +
  geom_contour(mapping = aes(colour = ..level..), bins = 20) +
  theme_bw() -> plot_Leigenvectorpq
#Uncomment the following chunk to export the figure as a ps file
# postscript(file = "Leigenvectorpq.ps")
# par(cex = 1.4)
# plot_Leigenvectorpq + coord_fixed()
# dev.off()

#Righteigenvector V
ent_database %>%
  select(p,q,s, init1, v1, v2) %>%
  filter(init1 == 0.5, s %in% c(0.5, 2, 5)) %>%  
  gather(key = Eigenvectors, value = value, -(1:4)) %>%
  ggplot(mapping = aes(x = p, y= q, z = value)) +
  facet_grid(Eigenvectors~s) +
  geom_contour(mapping = aes(colour = ..level..), bins = 20) +
  theme_bw() -> plot_Reigenvectorpq
#Uncomment the following chunk to export the figure as a ps file
# postscript(file = "Reigenvectorpq.ps")
# par(cex = 1.4)
# plot_Reigenvectorpq + coord_fixed()
# dev.off()


#as functions of s
ent_database %>%
  select(p,q,s, init1, u1, u2, v1, v2) %>%
  filter(p == 0.2, q== 0.3, init1 == 0.5) %>%  
  gather(key = Eigenvectors, value = value, -(1:4)) %>%
  ggplot(mapping = aes(x = s, y= value, group = Eigenvectors)) +
  geom_path(mapping = aes(colour = Eigenvectors, linetype = Eigenvectors)) +
  theme_bw() -> plot_eigevectorss
#Uncomment the following chunk to export the figure as a ps file
# postscript(file = "eigenvectorss.ps")
# par(cex = 1.4)
# plot_eigenvectorss + coord_fixed()
# dev.off()


#### Graphical representation of c ####

## As a function of p,q
ent_database %>%
  select(p,q,s, init1, c) %>%
  filter(init1 %in% c(0.2, 0.5, 0.8), s %in% c(0.5, 2, 5)) %>%  
  gather(key = Eigenvectors, value = value, -(1:4)) %>%
  ggplot(mapping = aes(x = p, y= q, z = value)) +
  facet_grid(init1~s) +
  geom_contour(mapping = aes(colour = ..level..), bins = 20) +
  theme_bw() -> plot_c
#Uncomment the following chunk to export the figure as a ps file
# postscript(file = 'c.ps')
# par(cex = 1.4)
# plot_c + coord_fixed()
# dev.off()

##As a function of s
ent_database %>%
  select(p,q,s, init1,c) %>%
  filter(p == 0.2, q== 0.3, init1 %in% c(0.2, 0.4, 0.5, 0.6, 0.8,1), s !=100) %>%  
  ggplot(mapping = aes(x = s, y= c)) +
  facet_wrap(~init1, scales='free') +
  geom_path() +
  theme_bw()

#### Graphical representations of entropy functions ####

##As functions of p, q for fixed values of s and init1
ent_database %>%
  select(p,q,s, init1, shannon, renyi, tsallis, taneja) %>%
  mutate(sc = as.character(s)) %>% #Convert s in character because sum(rep(0.1, 12)) == 1.2 returns FALSE !!!
  filter(init1 == 0.5, sc %in% as.character(c(0.2, 0.5, 0.8, 1.2, 2))) %>%  
  gather(key = entropy, value = value, - (1:4), -sc) %>% 
  # filter(entropy != 'tsallis') %>%
  mutate(entropy = factor(entropy, levels = c('shannon', 'renyi', 'tsallis', 'taneja'), ordered = TRUE)) %>%
  ggplot(mapping = aes(x = p, y = q, z = value)) +
  facet_grid(entropy~s, scales = 'fixed') +
  geom_contour(aes(colour = ..level..), bins = 20) +
  theme_bw() -> plot_entropies

#Shannon
ent_database %>%
  select(p,q,s, init1, shannon) %>%
  filter(init1 == 0.5, s == 1) %>%
  select(p,q,shannon) %>%  
  ggplot(mapping = aes(x = p, y = q, z = shannon)) +
  geom_contour(aes(colour = ..level..), bins = 20) +
  theme_bw() -> plot_shannon
plot_shannon

#Renyi for s = 0.5, 2, 10, 100
ent_database %>%
  select(p,q,s, init1, renyi) %>%
  filter(init1 == 0.5, s %in% c(0.5, 1, 2, 10)) %>%
  ggplot(mapping = aes(x = p, y = q, z = renyi)) +
  facet_wrap(~s) + 
  geom_contour(aes(colour = ..level..), bins = 20) +
  theme_bw() -> plot_renyi
#Uncomment the following chunk to export the figure as a ps file
# postscript(file = "renyirates.ps")
# plot_renyi + coord_fixed()
# dev.off()

#Tsallis for s = 0.2, 0.5, 0.8, init1 = 0.2, 0.5, 0.8
ent_database %>%
  select(p,q,s, init1, tsallis) %>%
  filter(init1 %in% c(0.2,0.5, 0.8), s %in% c(0.2, 0.5, 0.8)) %>%
  ggplot(mapping = aes(x = p, y = q, z = tsallis)) +
  facet_grid(s~init1) + 
  geom_contour(aes(colour = ..level..), bins = 20) +
  theme_bw() -> plot_tsallis
plot_tsallis
#Uncomment the following chunk to export the figure as a ps file
# postscript(file = "tsallisrates.ps")
# plot_tsallis + coord_fixed()
# dev.off()

#Taneja
ent_database %>%
  select(p,q,s, init1, taneja) %>%
  filter(init1 %in% c(0.2,0.5, 0.8), s %in% c(0.5, 2, 5)) %>%
  ggplot(mapping = aes(x = p, y = q, z = taneja)) +
  facet_grid(init1~s) + 
  geom_contour(aes(colour = ..level..), bins = 20) +
  theme_bw() -> plot_taneja
#Uncomment the following chunk to export the figure as a ps file
# postscript(file = "tanejarates.ps")
# plot_taneja + coord_fixed()
# dev.off()

# #### Graphical representation of entropy as functions of the initial distribution ####
# ##For p=0.2, q=0.3 and several values of s
# ent_database %>%
#   select(p,q,s,init1, shannon, renyi, tsallis, taneja) %>%
#   mutate(sc = as.character(s)) %>% #Convert s in character because sum(rep(0.1, 12)) == 1.2 returns FALSE !!!
#   filter(p == 0.2, q==0.3, sc %in% as.character(c(0.2, 0.5, 0.8, 1.2, 2, 5))) %>%  
#   gather(key = entropy, value = value, - (1:4), -sc) %>% 
#   mutate(entropy = factor(entropy, levels = c('shannon', 'renyi', 'tsallis', 'taneja'), ordered = TRUE)) %>%
#   ggplot(mapping = aes(x = init1, y = value, group = entropy)) +
#   facet_wrap(~s, scales = 'fixed') +
#   geom_path(aes(colour = entropy, linetype = entropy)) +
#   theme_bw()
# 
# ent_database %>%
#   filter(s == 100, init1 == 0.5) %>%
#   select(p,q,renyi) %>%
#   ggplot() +
#   geom_contour(mapping = aes(x = p, y = q, z = renyi, colour = ..level..), bins = 20) +
#   theme_bw()
