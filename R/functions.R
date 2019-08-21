#' Builds a movement network from transfer data, incorporating length of stay at node as well as length 'away' from node.
#'
#' @param df Data table of movement information. Must have a very specific structure, see \code{\link{Details}}
#' @param id Column containing individual identifier (unquoted)
#' @param node.column Column containing the node of the individual (unquoted)
#' @param time.column Column containing the time differences of the individual (unquoted).
#' @param exit.column Column containing entries whether line is an entry or exit to node (unquoted).
#' @param exit.code Code which indicates line is an exit (quoted).
#' @param nodes.list List containing all possible nodes (possibly larger than the list contained in df).
#' @param visit.column Column containing which visit number per individual (unquoted).
#' @param R_A Reproductive number of condition.
#' @param gamma Rate of recovery of condition.
#' @param parallel Compute in parallel (default \code{FALSE}).
#' @param num.cores Number of cores to use in parallel computation (ignored if \code{parallel == FALSE}).
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
#' Specifying \code{parallel=T} is only useful if each individual visits a large number of nodes.
#' If the number of nodes visited is small, then the overhead of running in parallel far outweighs the
#' parallelisation benefit.
#'
#' @return Adjacency matrix of the transfer network/graph.
buildNetwork <- function(df, id, node.column, time.column, exit.column, exit.code, visit.column, nodes.list,
                         R_A, gamma,
                         parallel = F, num.cores = parallel::detectCores()-1) {
  if (any(class(dt) != "data.table"))
    stop("Input df must be a data.table")

  if (parallel)
    warning("Specified parallel = T. Only use this option if each individual visits a large number of nodes.")

  id <- deparse(substitute(id))
  node.column <- deparse(substitute(node.column))
  time.column <- deparse(substitute(time.column))
  exit.column <- deparse(substitute(exit.column))
  visit.column <- deparse(substitute(visit.column))

  lbar <- df[get(exit.column) == exit.code, mean(get(time.column), na.rm=T)]

  if (parallel) {
    cl <- parallel::makeCluster(num.cores)
    parallel::clusterExport(cl, c("df", "id", "node.column", "time.column", "nodes.list",
                                  "lbar", "R_A", "gamma", ".calcNetworkForIndividual"), envir=environment())

    parallel::clusterEvalQ(cl, library(data.table))

    m_ij <- df[, list(m = parallel::parLapply(cl, X=.SD, fun = function(x) .calcNetworkForIndividual(x, node.column, time.column, exit.column, exit.code, visit.column, nodes.list, lbar, R_A, gamma))), by=id]

    parallel::stopCluster(cl)
  }
  else {
    m_ij <- df[, list(m = list(.calcNetworkForIndividual(.SD, node.column, time.column, exit.column, exit.code, visit.column, nodes.list, lbar, R_A, gamma))), by=id]
  }

  m <- Reduce('+', m_ij[, m])

  return (m)
}


.calcNetworkForIndividual <- function(df, node.column, time.column, exit.column, exit.code, visit.column, nodes.list, lbar, R_A, gamma)
{
  n <- nrow(df)
  x_ij <- matrix(0, nrow=length(nodes.list), ncol=length(nodes.list))
  max_visits <- max(df[, get(visit.column)])
  i <- df[get(exit.column) == exit.code & get(visit.column) < max_visits, vapply(get(node.column), function(x) { which(x == nodes.list)}, FUN.VALUE = 0)]
  j <- df[get(exit.column) == exit.code & get(visit.column) > 1, vapply(get(node.column), function(x) { which(x == nodes.list)}, FUN.VALUE = 0)]

  u_i <- 1-exp(-(df[get(exit.column) == exit.code & get(visit.column) < max_visits, get(time.column)]) ) #P(Contracting infection) = length of stay in source ward (don't want the last visit here)
  l_j <- df[get(exit.column) == exit.code & get(visit.column) > 1, get(time.column)]
  w_j <- 1-exp(-R_A*l_j/lbar)
  delta_t_ij <- df[get(exit.column) != exit.code & get(visit.column) > 1, get(time.column)]
  v_ij <- exp(-gamma*delta_t_ij)

  u_i[is.na(u_i)] <- 0
  w_j[is.na(w_j)] <- 0
  v_ij[is.na(v_ij)] <- 0

  for (iter in 1:length(i)) {
    x_ij[i[iter],j[iter]] <- x_ij[i[iter],j[iter]] + (u_i[iter]*w_j[iter]*v_ij[iter])
  }

  return(x_ij)
}


#' Calculates stay lengths for individuals in a given node.
#'
#' @param df Data table of movement information. Must have a very specific structure, see \code{\link{Details}}
#' @param id Column containing individual identifier (unquoted)
#' @param node.column Column containing the node of the individual (unquoted)
#' @param time.column Column containing the time differences of the individual (unquoted).
#' @param entry.column Column identifying whether a row is an entry or an exit into a node (unquoted).
#' @param entry.code Code identifying whether row is an entry into node (quoted).
#' @param order Whether the resulting \code{tibble} should be ordered by mean stay (useful for plotting)
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
#' @return data.table containing nodes and the time spent in each node. Note that there will be one row per visit (*not* per ward).
calcStayLengths <- function(df, id, node.column, time.column, entry.column, entry.code, order=T) {
  id <- deparse(substitute(id))
  node.column <- deparse(substitute(node.column))
  time.column <- deparse(substitute(time.column))
  entry.column <- deparse(substitute(entry.column))

  stay_lengths <- df[get(entry.column) == entry.code]

  if (order) {
    order <- stay_lengths[, .(mean_time=mean(get(time.column), na.rm=T)), by=node.column][order(mean_time), get(node.column)]
    print(order)
    stay_lengths <- stay_lengths[, eval(node.column):=factor(get(node.column), levels = order)]
  }

  return (stay_lengths)
}

#' Calculates number of entries into a given node
#'
#' @param df Data table of movement information. Must have a very specific structure, see \code{\link{Details}}
#' @param id Column containing individual identifier (unquoted)
#' @param node.column Column containing the node of the individual (unquoted)
#' @param time.column Column containing the time differences of the individual (unquoted).
#' @param entry.column Column containing information on whether the row is an entry or exit from the node (unquoted).
#' @param entry.code String code identifying whether row is an entry or exit from the node (quoted).
#' @param order Whether the resulting \code{tibble} should be ordered by number of entries (useful for plotting)
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
#' @return Tibble containing nodes and the number of entries into each.
calcNumEntriesIntoNode <- function(df, id, node.column, time.column, entry.column, entry.code, order=T) {
  id <- deparse(substitute(id))
  node.column <- deparse(substitute(node.column))
  time.column <- deparse(substitute(time.column))
  entry.column <- deparse(substitute(entry.column))

  next_event <- df[get(entry.column) == entry.code & (is.na(get(time.column)) | get(time.column) > 0), .N, by=node.column]

  if (order) {
    order <- next_event[order(N), get(node.column)]
    next_event <- next_event[, eval(node.column) := factor(get(node.column), levels=order)]
  }

  return (next_event)
}


#' Calculates number of entries into a given node
#'
#' @param df Data table of movement information. Must have a very specific structure, see \code{\link{Details}}
#' @param id Column containing individual identifier (unquoted)
#' @param node.column Column containing the node of the individual (unquoted)
#' @param time.column Column containing the time differences of the individual (unquoted).
#' @param entry.column Column containing information on whether the row is an entry or exit from the node (unquoted).
#' @param entry.code String code identifying whether row is an entry or exit from the node (quoted).
#' @param order Whether the resulting \code{tibble} should be ordered by mean time (useful for plotting)
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
#' @return data.table containing nodes and the time to next entry into a given node. Note that there is a row *per entry* not per node.
calcTimeBetweenEntries <- function(df, id, node.column, time.column, entry.column, entry.code, order=T) {
  id <- deparse(substitute(id))
  node.column <- deparse(substitute(node.column))
  time.column <- deparse(substitute(time.column))
  entry.column <- deparse(substitute(entry.column))

  time_between_entries <- df[get(entry.column) == entry.code & (!is.na(get(time.column)) & (get(time.column)>0))]

  if (order) {
    order <- time_between_entries[, .(mean_time = mean(get(time.column), na.rm=T)), by = node.column][order(mean_time), get(node.column)]
    time_between_entries <- time_between_entries[, eval(node.column) := factor(get(node.column), levels=order)]
  }

  return (time_between_entries)
}


#' Calculates number of direct transfers between nodes
#'
#' @param df Data table of movement information. Must have a very specific structure, see \code{\link{Details}}
#' @param id Column containing individual identifier (unquoted)
#' @param node.column Column containing the node of the individual (unquoted)
#' @param time.column Column containing the time differences of the individual (unquoted).
#' @param entry.column Column containing information on whether the row is an entry or exit from the node (unquoted).
#' @param entry.code String code identifying whether row is an entry or exit from the node (quoted).
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
#' @return data.table containing IDs and number of direct transfers between nodes (grouped by number of transfers).
calcNumberOfTransfers <- function(df, id, node.column, time.column, entry.column, entry.code) {
  id <- deparse(substitute(id))
  node.column <- deparse(substitute(node.column))
  time.column <- deparse(substitute(time.column))
  entry.column <- deparse(substitute(entry.column))

  number_of_changes <- df[ get(entry.column) == entry.code & (get(time.column) == 0 | is.na(get(time.column))), .(n=(.N)-1), by = id]
}
