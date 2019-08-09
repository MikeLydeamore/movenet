#' Builds a movement network from transfer data, incorporating length of stay at node as well as length 'away' from node.
#'
#' @param df Data frame of movement information. Must have a very specific structure, see \code{\link{Details}}
#' @param id Column containing individual identifier (unquoted)
#' @param node.column Column containing the node of the individual (unquoted)
#' @param time.column Column containing the time differences of the individual (unquoted).
#' @param nodes.list List containing all possible nodes (possibly larger than the list contained in df)
#' @param R_A Reproductive number of condition
#' @param gamma Rate of recovery of condition
#' @param parallel Compute in parallel (default \code{FALSE})
#' @param num.cores Number of cores to use in parallel computation (ignored if \code{parallel == FALSE})
#'
#' @details
#' The structure of \code{df} is very specific. For a given individual, each stay will have 2 rows (in order):
#' the time of admission, and the time of discharge. The node will be the same for  these 2 rows.
#' The \code{time.column} will then have either the stay length at a given node (even rows), or the time to
#' next admission (odd rows). The first row of this column is expected to be NA.
#'
#' The frame should be ordered by \code{id}, and then by the time (either admission or discharge). Sorting by \code{id}
#' is not required, but does make the frame nicer to look at.
#'
#' @return Adjacency matrix of the transfer network/graph.
buildNetwork <- function(df, id, node.column, time.column, nodes.list,
                         R_A, gamma,
                         parallel = F, num.cores = parallel::detectCores()-1) {
  id <- enquo(id)
  node.column <- enquo(node.column)
  time.column <- enquo(time.column)

  lbar <- df %>% slice(seq(2, n(), 2)) %>% summarise(mean(!!time.column, na.rm=T)) %>% as.numeric()

  if (parallel) {
    cl <- parallel::makeCluster(num.cores)
    parallel::clusterExport(cl, c("df", "id", "node.column", "time.column", "nodes.list",
                                  "lbar", "R_A", "gamma", "calcNetworkForIndividual"), envir=environment())

    parallel::clusterEvalQ(cl, library(dplyr))

    m_ij <- parallel::parLapply(cl, X=unique(df %>% select(!!id) %>% unique() %>% pull()),
                                fun = function(idx) {
                                  df_id <- filter(df, !!id == idx)
                                  .calcNetworkForIndividual(df_id, node.column, time.column, nodes.list, lbar, R_A, gamma)
                                })

    parallel::stopCluster(cl)
  }
  else {
    m_ij <- lapply(X = unique(df %>% select(!!id) %>% unique() %>% pull()),
                    FUN = function(idx) {
                     df_id <- filter(df, !!id == idx)
                     .calcNetworkForIndividual(df_id, node.column, time.column, nodes.list, lbar, R_A, gamma)
                   })
  }

  m <- Reduce('+', m_ij)

  return (m)
}


.calcNetworkForIndividual <- function(df, node.column, time.column, nodes.list, lbar, R_A, gamma)
{
  n <- nrow(df)
  x_ij <- matrix(0, nrow=length(nodes.list), ncol=length(nodes.list))

  i <- sapply(df %>% select(!!node.column) %>% slice(seq(2, n-2, 2)) %>% pull(), function(x) which(x == nodes.list)) %>% as.vector()
  j <- sapply(df %>% select(!!node.column) %>% slice(seq(4, n, 2)) %>% pull(), function(x) which(x == nodes.list)) %>% as.vector()

  df.time <- df %>% select(!!time.column)
  u_i <- 1-exp(-(df.time %>% slice(seq(2, n-2, 2)) %>% pull())) #P(Contracting infection) = length of stay in source ward (don't want the last visit here)
  l_j <- df.time %>% slice(seq(4, n, 2)) %>% pull() #P(transmitting infection) = length of stay in destination ward (don't want the first visit here)
  w_j <- 1-exp(-R_A*l_j/lbar)
  delta_t_ij <- df.time %>% slice(seq(3, n, 2)) %>% pull()
  v_ij <- exp(-gamma*delta_t_ij)

  u_i[is.na(u_i)] <- 0
  w_j[is.na(w_j)] <- 0
  v_ij[is.na(v_ij)] <- 0

  for (iter in 1:length(i)) {
    x_ij[i[iter],j[iter]] <- x_ij[i[iter],j[iter]] + (u_i[iter]*w_j[iter]*v_ij[iter])
  }


  return(x_ij)
}
