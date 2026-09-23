# Internal readable aliases for the registered HDD acceleration routines.
.hdzacp_hdd_distance_cpp <- function(x) hdzacp_hdd_distance_matrix_cpp(x)
.hdzacp_hdd_stat_cpp <- function(distance, order) {
  hdzacp_hdd_stat_from_distance_cpp(distance, order)
}
