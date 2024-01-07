# sir_abm
SIR models on networks

# Demo

## One network

Load the software: 

```R
source("sim.R")
```

Create a network:

```R
n_agents = 1e2
n_connections = 1e3
g = sample_gnm(n = 1e2, m = 1e3)
```

Prepare it by adding SIR status to each node: 

```R
gm = prepare_network(g)
```

Run the simulation: 

```R
sir_abm(gm)
```

which returns a `data.frame`:

```R
> sir_abm(gm)
        ate         se            network_type n_agents n_edges
1 0.3575105 0.04941228 Erdos-Renyi (gnm) graph      100    1000
```

## Many networks

Create a list of graphs: 

```R
er_graphs  = replicate(n = 1e2,
                       expr = sample_gnm(n = 1e2, m = 1e3)
)
```

Prepare them: 

```R
er_graphs = lapply(er_graphs, prepare_network)
```

This way you can go back and collect any features of the networks after the simulation (since the simulations do not alter the network structure). 

Then run the simulations (100 graphs took about 20 seconds to run): 

```R
er_list = run_sir(er_graphs) # run_sir() shows a progress bar and returns a list
```

and combine everything into a `data.frame`:

```R
er_df = dplyr::bind_rows(er_list)
> mean(er_df$ate)
[1] 0.5381295
```







