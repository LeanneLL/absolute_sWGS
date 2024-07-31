library(rascal)

find_acceptable_solutions <- function(relative_copy_numbers, weights = NULL,
                                    min_ploidy = 1.1, max_ploidy = 5.5, ploidy_step = 0.01,
                                    min_cellularity = 0.01, max_cellularity = 1.0, cellularity_step = 0.01,
                                    distance_function = c("MAD", "RMSD"),
                                    distance_filter_scale_factor = 1.25,
                                    max_proportion_zero = 0.05,
                                    min_proportion_close_to_whole_number = 0.5,
                                    max_distance_from_whole_number = 0.15,
                                    solution_proximity_threshold = 5,
                                    keep_all = FALSE) {
  
  # compute distance function for grid-based range of ploidies and cellularities
  distances <- absolute_copy_number_distance_grid(
    relative_copy_numbers, weights,
    min_ploidy, max_ploidy, ploidy_step,
    min_cellularity, max_cellularity, cellularity_step,
    distance_function)
  
  # add identifier for each solution to help with marking up the best fit
  # solutions if returning all solutions
  distances <- distances %>%
    mutate(id = row_number())
  
  # find minima, i.e. those grid points which have a smaller distance than all
  # adjacent points
  distances <- distances %>%
    mutate(x = dense_rank(cellularity)) %>%
    mutate(y = dense_rank(ploidy))
  
  solutions <- select(distances, id, x, y, distance)
  
  for (xdelta in -1:1) {
    for (ydelta in -1:1) {
      if (xdelta != 0 || ydelta != 0) {
        solutions <- solutions %>%
          mutate(xc = x + xdelta, yc = y + ydelta) %>%
          left_join(select(distances, xc = x, yc = y, dc = distance), by = c("xc", "yc")) %>%
          filter(is.na(dc) | distance <= dc) %>%
          select(id, x, y, distance)
      }
    }
  }
  
  distances <- select(distances, -x, -y)
  
  solutions <- solutions %>%
    select(id) %>%
    left_join(distances, by = "id")
  
  # only retain solutions with distances no more than
  # distance_filter_scale_factor times the the minimum value
  if (is.numeric(distance_filter_scale_factor) && nrow(solutions) > 1) {
    solutions <- solutions %>%
      filter(distance < distance_filter_scale_factor * min(distance))
  }
  
  # retain solutions with acceptable ploidies and cellularities
  solutions <- solutions %>%
    rowwise() %>%
    filter(
      is_acceptable_ploidy_and_cellularity(
        ploidy, cellularity,
        relative_copy_numbers, weights,
        max_proportion_zero = max_proportion_zero,
        min_proportion_close_to_whole_number = min_proportion_close_to_whole_number,
        max_distance_from_whole_number = max_distance_from_whole_number
      )
    ) %>%
    ungroup()
  
  return(solutions)
}