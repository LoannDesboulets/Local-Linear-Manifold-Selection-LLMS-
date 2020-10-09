## -----------------------------------------------------------
## Clear Memory
rm(list=ls())

## -----------------------------------------------------------
## Load the function to simulate data and perform LLMS
source("Simulate_Manifold.R")
source("LLMS.R")

## -----------------------------------------------------------
## Simulate some data
X <- Simulate_Manifold()

## -----------------------------------------------------------
## Invoke the function
result <- LLMS(X)

## -----------------------------------------------------------
## Plot the Selection Paths
matplot(result$Paths,type="l",lwd=2)

## -----------------------------------------------------------
## Plot the Sparse Graphical Model
result <- LLMS(X,graph.model=TRUE) # Rerun the algorithm with graph.model
source("Plot_Graph.R")             # Customized function for plotting
G <- result$Graph[,,100]           # The graph with the maximum penalty
Plot_Graph(G)
  