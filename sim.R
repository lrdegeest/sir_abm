library(progress)
library(igraph)
library(survival)
library(survminer)


sir_abm = function(network,
                    transmission_prob = 0.02,
                    recovery_prob = 0.01,
                    vaccine_effectiveness = 0.5,
                    sim_steps = 1e2,
                    plot = FALSE,
                    return_data = FALSE){
  #' Conducts an agent-based simulation of disease transmission and recovery.
  #'
  #' @param network The network on which the simulation will be conducted.
  #' @param transmission_prob The probability of disease transmission between connected agents (default is 0.02).
  #' @param recovery_prob The probability of an infected agent recovering in each simulation step (default is 0.01).
  #' @param vaccine_effectiveness The effectiveness of the vaccine in reducing transmission probability (default is 0.5).
  #' @param sim_steps The number of simulation steps (default is 1e2).
  #' @param plot A logical value. If TRUE, generates a survival plot based on simulation results (default is FALSE).
  #' @param return_data A logical value. If TRUE, returns raw data used in the simulation (default is FALSE).
  #'
  #' @return If `return_data` is TRUE, returns a dataframe containing simulation data.
  #'         If `plot` is TRUE, returns a survival plot. Otherwise, returns a dataframe with summary statistics.
  #'
  #' @details This function simulates disease spread using a Susceptible-Infectious-Recovered (SIR) model on the provided network.
  #'          It tracks infection periods, treatment status, and recovery status of agents.
  #'
  #' @examples
  #' # Example usage:
  #' result <- sir_abm2(my_network, transmission_prob = 0.02, recovery_prob = 0.01, vaccine_effectiveness = 0.5, sim_steps = 100, plot = TRUE)
  #'
  #' @import ggplot2
  #' @import survival
  #' @export
  #' 
  # get the number of agents
  num_agents <- vcount(network)
  
  # vaccine's effect on transmission
  vt <- vaccine_effectiveness*transmission_prob
  
  # list to store data for n-agents
  # each agent's entry is a triple
  ## 1. period in which agent gets infected (always 'NA' for `initial_infected_agents`)
  ## 2. treatment status (whether the agent did or did not receive the vaccine)
  ## 3. recovery status (all agents recover at the end; just need this for survival model)
  surv_dat <- vector(mode = 'list', length = num_agents)
  agent_degree <- degree(network)
  for (i in 1:num_agents){
    surv_dat[[i]] <- c(NA,0,0)
  }
  
  # store agent status
  agent_state <- V(network)$state
  # initial infected agents
  initial_infected_agents <- which(agent_state == 'Infectious')
  # store neighborhood for each agent
  neighborhoods <- lapply(1:num_agents, function(x) neighbors(network, x))
  
  # begin simulation
  # begin simulation
  for (step in 1:sim_steps) {
    for (agent in 1:num_agents) {
      # if the agent is infections
      if (agent_state[agent] == "Infectious") {
        # find his neighbors
        neighbors <- neighborhoods[[agent]]
        # for each neighbor
        for (neighbor in neighbors) {
          # if he is susceptible and unvaccinated
          if (agent_state[neighbor] == "Susceptible") {
            # infect him by random chance
            if (runif(1) < transmission_prob) {
              surv_dat[[neighbor]] <- c(step,0,1)
              agent_state[neighbor] <- "Infectious"
            }
          }
          # or if he is susceptible but vaccinated 
          else if (agent_state[neighbor] == "Susceptible Vaccinated") {
            # infect him by random chance -- cut by vaccine effectiveness
            if (runif(1) < vaccine_effectiveness*transmission_prob) {
              surv_dat[[neighbor]] <- c(step,1,1)
              agent_state[neighbor] <- "Infectious"
            }
          }
        }
        # now turn back to the infections agent in question
        # by chance he recovers this sweep 
        if (runif(1) < recovery_prob) {
          agent_state[neighbor] <- "Recovered"
        }
      }
    }
  }
  
  # convert to matrix
  surv_dat_tot <- do.call(rbind,surv_dat)
  # convert to dataframe
  surv_dat_df <- data.frame(surv_dat_tot)
  colnames(surv_dat_df) <- c("time","trt","status")
  # add the agent's id
  surv_dat_df$agent_id <- 1:num_agents
  # add the agent's degree
  surv_dat_df$agent_degree <- degree(network)
  # remove the initially infected agents (can't have NA in survival model -- tho these aren't always NA)
  surv_dat_df <- surv_dat_df[setdiff(1:num_agents, initial_infected_agents),]
  
  
  if(return_data) return(surv_dat_df)
  
  km_trt_fit = coxph(Surv(time,status) ~ trt, data=surv_dat_df)
  
  if(plot){
    surv_curve <- survfit(Surv(time,status) ~ trt, data=surv_dat_df)
    p <- ggsurvplot(surv_curve,
                   conf.int = TRUE,
                   legend.labs = c("Unvaccinated", "Vacinnated"),
                   legend.title = "") + 
      labs(subtitle = paste0('Treatment effect: ', 
                             round(exp(km_trt_fit$coefficients),3)))
    return(p)
  } 
  else{
    #return_list = list()
    #return(exp(km_trt_fit$coefficients))
    results_df <- data.frame(
      ate = exp(coefficients(km_trt_fit)),
      se = as.numeric(vcov(km_trt_fit)),
      network_type = network$name,
      n_agents = vcount(network),
      n_edges = ecount(network)
    )
    rownames(results_df) = NULL
    return(results_df)
  }
}




prepare_network = function(graph,
                           initial_infected = 2){
  #' Prepares a network for an agent-based simulation of disease spread.
  #'
  #' @param graph The network/graph structure for the simulation.
  #' @param initial_infected The number of initially infected agents (default is 2).
  #'
  #' @return Returns the modified network with initial states assigned to agents for simulation purposes.
  #'
  #' @details This function prepares the network by assigning initial states to agents:
  #'          - All agents are initially set to "Susceptible".
  #'          - A specified number of agents are randomly infected ("Infectious").
  #'          - Half of the susceptibles are randomly chosen to be "Susceptible Vaccinated".
  #'            Vaccinated agents cannot be infected at the simulation's start.
  #'
  #' @examples
  #' # Example usage:
  #' my_graph <- prepare_network(my_network, initial_effected = 5)
  #'
  #' @export
  #' 
  num_agents <- vcount(graph)
  
  # everyone is susceptible
  V(graph)$state <- "Susceptible"
  # but a few are randomly infected
  initial_infected_agents <- sample(1:num_agents, size = initial_infected)
  V(graph)$state[initial_infected_agents] <- "Infectious"
  
  # Set half of susceptibles to vaccinated
  vaccinated_susceptibles <- sample(1:num_agents, size = num_agents/2)
  # vaccinated agents can't be infected to start
  vaccinated_susceptibles <- setdiff(vaccinated_susceptibles,initial_infected_agents)
  V(graph)$state[vaccinated_susceptibles] <- "Susceptible Vaccinated"
  
  return(graph)
  
  
}


mc_sir = function(n_sims = 1e3){
  #' Runs multiple Monte Carlo simulations of an Agent-Based SIR model.
  #'
  #' @param n_sims The number of simulations to run (default is 1000).
  #'
  #' @return Returns a list of results from the Monte Carlo simulations.
  #'
  #' @details This function conducts a Monte Carlo simulation of an Agent-Based Susceptible-Infectious-Recovered (SIR) model.
  #'          It runs multiple iterations of the `sir_abm()` function and collects the results.
  #'
  #' @examples
  #' # Example usage:
  #' result_list <- mc_sir(n_sims = 500)
  #'
  #' @import progress
  #' @export
  #' 
  ate_list <- vector(mode = 'list', length = n_sims)
  pb <- progress_bar$new(total = n_sims,
                        clear = FALSE, 
                        width= 50,
                        format = " 🚴💨 [:bar] :percent eta: :eta")
  for(i in 1:n_sims){
    ate[[i]] <- try(sir_abm())
    pb$tick()
  }
}

run_sir = function(network_list){
  #' Runs an Agent-Based SIR simulation on a list of networks.
  #'
  #' @param network_list A list of network structures for simulation.
  #'
  #' @return Returns a list of results from the Agent-Based SIR simulations.
  #'
  #' @details This function conducts Agent-Based Susceptible-Infectious-Recovered (SIR) simulations
  #'          on a list of provided network structures.
  #'          It iterates through each network in the list, running the `sir_abm()` function,
  #'          and collects the simulation results.
  #'
  #' @examples
  #' # Example usage:
  #' my_networks <- list(network1, network2, network3)
  #' simulation_results <- run_sir(network_list = my_networks)
  #'
  #' @import progress
  #' @export
  #' 
  results_list <- vector(mode = 'list', length = length(network_list))
  pb = progress_bar$new(total = length(results_list),
                        clear = FALSE, 
                        width= 50,
                        format = " 🚴💨 [:bar] :percent eta: :eta")
  pb$tick(0)
  for(i in seq_along(network_list)){
    res <- try(sir_abm(network = network_list[[i]]))
    results_list[[i]] <- res
    pb$tick()
  }
  return(results_list)
}
