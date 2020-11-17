#' Reduce a GTFS file to a subset based on route_type
#'
#' @param gtfs_in tidytransit GTFS object
#' @param numeric route_type to reduce the GTFS data to
#'
#' @return A tidytransit GTFS object filtered by the selected route type
#' @export
#'
#' @examples
#'
reduce_gtfs_to_rail <- function(gtfs_in, route_type){


  gtfs <- gtfs_in
  gtfs$routes <- gtfs$routes[gtfs$routes$route_type==route_type, ]
  gtfs$trips <- gtfs$trips[gtfs$trips$route_id %in% gtfs$routes$route_id,]
  gtfs$stop_times <- gtfs$stop_times[gtfs$stop_times$trip_id %in% gtfs$trips$trip_id,]
  gtfs$stops <- gtfs$stops[gtfs$stops$stop_id %in% gtfs$stop_times$stop_id,]
  gtfs$transfers <- gtfs$transfers[(gtfs$transfers$from_stop_id %in% gtfs$stop_times$stop_id) & (gtfs$transfers$to_stop_id %in% gtfs$stop_times$stop_id), ]
  gtfs$shapes <- gtfs$shapes[gtfs$shapes$shape_id %in% gtfs$trips$shape_id,]

  return(gtfs)
}
