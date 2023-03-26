library(pacman)
p_load(tidyverse)
library(MathBioSim)
library(googledrive)
name = "smoothing"
#dd<-0.2
for(initial_pop in c(1, 10, 50, 100)) {
for(death_r in c(0.6, 0.4)) {
for(b in c(0.6, 0.4)){
for(dd in c(0.2, 0.4)){
for (auto_stop_at_plateau in c(TRUE, FALSE)){
sd = list(sd_b=1, sd_d=1)
for(sd in list(list(0.6, 0.6), list(0.7, 0.2), list(0.2, 0.7))){
print(sd[[1]])
#for (sd_b in c(0.6)){
#b<-0.2
#death_r<-0.02
#auto_stop_at_plateau<-TRUE
epochs_count <- 2000 
plateau_threshold<-8
area_length_x<-100
#initial_pop <- 10
sim<-initialize_simulator(area_length_x = area_length_x, dd=dd,
                          initial_population_x = runif(initial_pop, 1, area_length_x), #c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
                          death_r = death_r,
                          death_y = dnorm(seq(0,5,length.out = 1001), sd = sd[[1]]),  #sd_d),
                          b=b,
                          birth_ircdf_y = qnorm(seq(0.5,1-1e-6,length.out = 101), sd = sd[[2]]), #sd_b
                          realtime_limit = 1000000000)
sim_results <- run_simulation(sim, epochs_count, calculate.pcf = TRUE, plateau_threshold = plateau_threshold, auto_stop_at_plateau = auto_stop_at_plateau)

#ggplot(sim_results$pattern,aes(x=x))+
#  geom_histogram(bins=100,breaks=seq(0,100,length.out = 101))
  
ggplot(sim_results$population,aes(x=time,y=pop))+
  geom_line()

ggplot(sim_results$exp_pop10, aes(x = time, y = exp_pop10)) +
  geom_line()

ggplot(sim_results$exp_pop2,aes(x=time,y=exp_pop2))+
  geom_line()

ggplot(sim_results$pcf,aes(x=r,y=pcf))+
  geom_line()

ggplot(sim_results$K,aes(x=r,y=K_iso))+
  geom_line()

# Запись результатов симуляции в файлы 
setwd("/Users/sokanaid")

folder <- "simulations_tables"
if (file.exists(folder)) {
  cat("The folder already exists")
} else {
  dir.create(folder)
}
setwd(folder)
if (file.exists(name)) {
  cat("The folder already exists")
} else {
  dir.create(name)
}
setwd(name)
folder <- paste(name, "dd_",dd,"_death_r",death_r,"_auto_stop_at_plateau_","_b_" ,b, auto_stop_at_plateau,"_initial_pop_", initial_pop,
                "_sd_b_", sd[[2]], "_sd_d_", sd[[1]],sep = "")
gsub(".", "_", folder)
#folder <- paste(folder, ".csv", sep="")
print("folder:")
print(folder)
if (file.exists(folder)) {
  cat("The folder already exists")
} else {
  dir.create(folder)
}
setwd(folder)
print("after set")
write.csv(sim_results$K, "K.csv")
write.csv(sim_results$population, "population.csv")
write.csv(sim_results$exp_pop, "exp_pop.csv")
write.csv(sim_results$exp_pop2, "exp_pop2.csv")
write.csv(sim_results$exp_pop2, "exp_pop3.csv")
write.csv(sim_results$exp_pop2, "exp_pop4.csv")
write.csv(sim_results$exp_pop2, "exp_pop5.csv")
write.csv(sim_results$exp_pop2, "exp_pop6.csv")
write.csv(sim_results$exp_pop2, "exp_pop7.csv")
write.csv(sim_results$exp_pop2, "exp_pop8.csv")
write.csv(sim_results$exp_pop2, "exp_pop9.csv")
write.csv(sim_results$exp_pop2, "exp_pop10.csv")
fileConn<-file("description.txt")
writeLines(c(paste("dd=",dd, sep="")
             ,paste("death_r=",death_r, sep="")
             , paste("epochs_count=", epochs_count, sep="")
             , paste("plateau_threshold=", plateau_threshold, sep="")
             , paste("area_length_x=",area_length_x, sep="")
             ,paste("initial_pop=",initial_pop, sep="")
             ,paste("auto_stop_at_plateau=",auto_stop_at_plateau, sep="")
             ,paste("b=",b, sep="")
             ,paste("sd_b=",sd[[2]], sep="")
             ,paste("sd_d=",sd[[1]], sep="")
             ,paste("found_plateau=",sim_results$found_plateau, sep="")
             ,paste("realtime_limit_reached=",sim_results$realtime_limit_reached, sep="")
             ), fileConn)
close(fileConn)
#}
}}} }}}
