
#' Simplified simulator runner
#'
#' @param simulator 1-species simulator object (1d, 2d or 3d) 
#' @param epochs amount of population-dependent time units to run
#' @param calculate.pcf if TRUE, adds pcf estimate to output
#' @param pcf_grid grid for pcf calculation
#' @param auto_stop_at_plateau stop automatically when it reaches the plateau
#' @param plateau_rate plateau acceptance threshold 
#' 
#' @return list with results
#' @export
#'
#' @examples
run_simulation<-
  function(simulator, epochs, calculate.pcf=FALSE, pcf_grid, auto_stop_at_plateau = TRUE, plateau_threshold = 10) {
  require(dplyr)
  time<-numeric(epochs+1)
  pop<-numeric(epochs+1)
  exp_pop<-rep(list(numeric(epochs+1)), each = 10)
  realtime_limit_reached<-FALSE
  
  pop[1] = simulator$total_population
  for(i in 1:length(exp_pop)) {
    exp_pop[i][1] = pop[1]
  }
  time[1] = simulator$time
  
  simulation_windows_sizes<-c(max(numeric(epochs+1) %/% 5, 1000),
   max(numeric(epochs + 1) %/% 10, 500))
  last_windows<-c(list(pop[1]), list(pop[1]))
  windows_var<-list(list(), list())
  alpha = 0.25 # константа экспаненциального сглаживания
  found_plateau = FALSE
  print("epochs")
  print(epochs)
  print(simulation_windows_sizes)
  print("simulation_windows_sizes")
  for(j in 2:(epochs + 1)){
    
    simulator$run_events(simulator$total_population)
    pop[j]=simulator$total_population
    time[j]=simulator$time
    for(i in 1:length(exp_pop)) {
      if(i==1){
        exp_pop[[i]][j] =  exp_pop[[i]][j-1]*(1-alpha) + pop[j]*alpha
      }
      else{
        exp_pop[[i]][j] = exp_pop[[i]][j-1]*(1-alpha) + exp_pop[[i-1]][j]*alpha
      }
    }
    if(auto_stop_at_plateau){
      # check the plateau
      for(window_size_ind in seq(1,length(simulation_windows_sizes))){
        if(j %% simulation_windows_sizes[[window_size_ind]]==0){
          var_val<-var(last_windows[[window_size_ind]])
          
          windows_var[[window_size_ind]]<-c(windows_var[[window_size_ind]], var_val)
          if(var_val<0.00000001) {
            # to escape zero division
            var_val<-0.00000001
          }
          prev_ind<-length(windows_var[[window_size_ind]])-1
          if(prev_ind>20) {
              prev_ind<- prev_ind-1
          }
          if(prev_ind >=1 && (windows_var[[window_size_ind]][[prev_ind]]/var_val > plateau_threshold || windows_var[[window_size_ind]][[prev_ind]]/var_val < 1/plateau_threshold)){ # 2 участка плато
            print("windows_var[[window_size_ind]][[prev_ind]]/var_val")
            print(windows_var[[window_size_ind]][[prev_ind]]/var_val)
            print("at plateuo")
            found_plateau = TRUE
            break
          }
          # ищем максимум из res
          # но так как пытаемся определить в realtime посмотрим на предыдущий и если текущий значительно меньше, то примем предыдущий за выход на плато
          if(prev_ind >= 1) {
            print("windows_var[[window_size_ind]][[prev_ind]]/var_val")
            print(windows_var[[window_size_ind]][[prev_ind]]/var_val)
          }
          last_windows[[window_size_ind]]<-exp_pop[[length(exp_pop)]][j]
          print("-")
        }
        else{
          last_windows[[window_size_ind]]<-c(last_windows[[window_size_ind]], exp_pop[[length(exp_pop)]][j])
        }
        
      }
      if(found_plateau) {
        break
      }
    }

  }
  
  if(simulator$total_population == 0){
    return(
      list('realtime_limit_reached' = simulator$realtime_limit_reached,
           'population' = data.frame(time=time,pop=pop)%>%distinct())
    )
  }
  
  if(class(simulator)[1]=='Rcpp_poisson_1d'){
    pattern = data.frame(x=simulator$get_all_x_coordinates())
  }else if(class(simulator)[1]=='Rcpp_poisson_2d'){
    pattern = data.frame(x=simulator$get_all_x_coordinates(),
                         y=simulator$get_all_y_coordinates())
  }else if(class(simulator)[1]=='Rcpp_poisson_3d'){
    pattern = data.frame(x=simulator$get_all_x_coordinates(),
                         y=simulator$get_all_y_coordinates(),
                         z=simulator$get_all_z_coordinates())
  }
  
  result <- list('realtime_limit_reached' = simulator$realtime_limit_reached,
                 'population' = data.frame(time=time,pop=pop)%>%distinct(),
                 'pattern' = pattern,
                 "found_plateau" = found_plateau,
                 "exp_pop"=data.frame(time=time[1:length(exp_pop[[1]])],exp_pop=exp_pop[[1]])%>%distinct(), 
                 "exp_pop2"=data.frame(time=time[1:length(exp_pop[[1]])],exp_pop2=exp_pop[[2]])%>%distinct(),
                 "exp_pop3"=data.frame(time=time[1:length(exp_pop[[1]])],exp_pop3=exp_pop[[3]])%>%distinct(),
                 "exp_pop4"=data.frame(time=time[1:length(exp_pop[[1]])],exp_pop4=exp_pop[[4]])%>%distinct(),
                 "exp_pop5"=data.frame(time=time[1:length(exp_pop[[1]])],exp_pop5=exp_pop[[5]])%>%distinct(),
                 "exp_pop5"=data.frame(time=time[1:length(exp_pop[[1]])],exp_pop6=exp_pop[[6]])%>%distinct(),
                 "exp_pop5"=data.frame(time=time[1:length(exp_pop[[1]])],exp_pop7=exp_pop[[7]])%>%distinct(),
                 "exp_pop5"=data.frame(time=time[1:length(exp_pop[[1]])],exp_pop8=exp_pop[[8]])%>%distinct(),
                 "exp_pop5"=data.frame(time=time[1:length(exp_pop[[1]])],exp_pop9=exp_pop[[9]])%>%distinct(),
                 "exp_pop10"=data.frame(time=time[1:length(exp_pop[[1]])],exp_pop10=exp_pop[[10]])%>%distinct())
  
  if (calculate.pcf){
    require(spatstat)
    if(class(simulator)[1]=='Rcpp_poisson_1d'){
      
      points<-ppp(simulator$get_all_x_coordinates(),
                  rep(0,simulator$total_population),
                  c(0,simulator$area_length_x),
                  c(-simulator$area_length_x/2,simulator$area_length_x/2)
      )
      
      if (missing(pcf_grid)){
        K_estimate<-Kest(points,correction="Ripley")   
        pcf_grid <- K_estimate$r
      }else{
        K_estimate<-Kest(points,r=pcf_grid,correction="Ripley") 
      }
      
      pcf_estimate <- 
        data.frame(Kest=K_estimate$iso/2,x=pcf_grid)%>%                
        mutate(pfc=(Kest-lag(Kest))/(pcf_grid-lag(pcf_grid))/simulator$area_length_x)%>% 
        pull(pfc)  
      
    }else if(class(simulator)[1]=='Rcpp_poisson_2d'){
      
      points<-ppp(simulator$get_all_x_coordinates(),
                  simulator$get_all_y_coordinates(),
                  c(0,simulator$area_length_x),
                  c(0,simulator$area_length_y)
      )
      
      if (missing(pcf_grid)){
        K_estimate<-Kest(points,correction="Ripley")   
        pcf_grid <- K_estimate$r
      }else{
        K_estimate<-Kest(points,r=pcf_grid,correction="Ripley") 
      }
      
      pcf_estimate <- pcf(K_estimate, method="b")$pcf
      
    }else if(class(simulator)[1]=='Rcpp_poisson_3d'){
      points<-pp3(simulator$get_all_x_coordinates(),
                  simulator$get_all_y_coordinates(),
                  simulator$get_all_z_coordinates(),
                  c(0,simulator$area_length_x),
                  c(0,simulator$area_length_y),
                  c(0,simulator$area_length_z)
      )
      
      if (missing(pcf_grid)){
        K_estimate<-Kest(points,correction="Ripley")   
        pcf_grid <- K_estimate$r
      }else{
        K_estimate<-Kest(points,r=pcf_grid,correction="Ripley") 
      }
      
      pcf_estimate <- pcf(K_estimate, method="b")$pcf
    }else{
      stop(glue(
        'simulator is not supported, object should be poisson_1d, poisson_2d or poisson_3d, received {cl}', 
        cl=class(simulator)[1]
      ))
    }
    
    result[['pcf']]<-data.frame(r=pcf_grid,pcf=pcf_estimate)
    result[['K']]<-data.frame(r=pcf_grid,K_iso=K_estimate$iso,K_poisson=K_estimate$theo)
}
  
  return(result)
}
