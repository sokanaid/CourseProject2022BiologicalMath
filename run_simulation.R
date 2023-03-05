
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
  function(simulator, epochs, calculate.pcf=FALSE, pcf_grid, auto_stop_at_plateau = TRUE, plateau_threshold = 0.8) {
  require(dplyr)
  time<-numeric(epochs+1)
  pop<-numeric(epochs+1)
  realtime_limit_reached<-FALSE
  
  pop[1] = simulator$total_population
  time[1] = simulator$time
  # 1. определить длительность отрезка здесь придется выбирать шаг не зависимо от симуляции (т к заранее не понятно сколько она продлится)
  # 2. убрать выбросы если нужно
  # 3. сравнивать отрезки
  
  # choose sizes of windows indenendetly on simulation len (because it is not known before the end of simulation)
  simulation_windows_sizes<-c(max(numeric(epochs+1) %/% 5, 1000),
   max(numeric(epochs + 1) %/% 10, 500))
  # TODO: проверить что здесь действительно правильно запоняетмся массив
  last_windows<-c(c(pop[1]), c(pop[1]), c(pop[1]))
  windows_mean<-c(list(), list(), list())
  windows_var<-c(list(), list(), list())
  windows_var<-c(-1,-1, -1)
  found_plateau = FALSE
  for(j in 2:(epochs + 1)){
    simulator$run_events(simulator$total_population)
    pop[j]=simulator$total_population
    if(auto_stop_at_plateau){
      # check the plateau
      for(window_size_ind in length(simulation_windows_sizes)){
        if(j %% simulation_windows_sizes[window_size_ind]==0){
          # full window
          # check plateou issue

          # TODO: здесь еще нужно выкинуть 5% максимумов и минимумов из последовательности
          mean_val<-mean(last_windows[window_size_ind])
          var_val<-var(last_windows[window_size_ind])
          windows_mean[window_size_ind]<-c(windows_mean[window_size_ind], mean_val)
          windows_var[window_size_ind]<-c(windows_var[window_size_ind], mean_val)
          if(mean_val<0.00000001) {
            # to escape zero division
            mean_val<-0.00000001
          }
          if(var_val<0.00000001) {
            # to escape zero division
            var_val<-0.00000001
          }
          prev_ind<-length(windows_mean[window_size_ind])-2
          if(prev_ind>=0 && abs(windows_mean[window_size_ind][prev_ind]/mean_val) > 0.8) {
              if(prev_ind>20) {
                prev_ind<- prev_ind-1
              }
              # ищем максимум из res
              # но так как пытаемся определить в realtime посмотрим на предыдущий и если текущий значительно меньше, то примем предыдущий за выход на плато
              res[window_size_ind]<-windows_var[window_size_ind][prev_ind]/var_val
          } else if(res[window_size_ind]!=-1 && res[window_size_ind] > windows_var[window_size_ind][prev_ind]/var_val) {
            # предыдущий отрезок был выходом на плато
            found_plateau = TRUE
            break
          }
          last_windows[window_size_ind]<-list()
        }
        else{
          last_windows[window_size_ind]<-c(last_windows[window_size_ind], pop[j])
        }
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
                 "found_plateau" = found_plateau)
  
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
