#' Produces isochrones as an sf-object
#'
#' @description gtfs_isochrones Creates one SP Object containg a isochrones calculated from a tidytransit::raptor result.
#' The function uses the minimal travel times from the start stop to all stops as a main input parameter.
#' This isochrone is a not meant to depict exact reachibility of transit in detail. Use opentripplanner for this instead.
#' The goal here is to provide a decent overview of transit access at larger scales (regional, national or continental).
#' IMPORTANT: The raptor_result table has to have only have one single stop as origin
#'
#' @param raptor_result A result from tidytransit::raptor function
#' @param stops A stops dataframe obtained from a tidytransit gtfs object: gtfs$stops
#' @param breaks A vector containing the travel time breaks for the isochrone in seconds.
#' @param hull_alpha_min The alpha values increase by increasing break interval for better graphic display of the isochrones.
#' @param buffer_value The shapes are buffered in the end for graphic display. Value in meters.
#'
#' @return An sp object with features corresponding to the isochrones of the breaks provided in 'breaks'.
#'
#'
#' @import data.table
#' @import alphahull
#' @import sf
#' @import sp
#' @import tidytransit
#' @import rgeos
#'
#'
#'
#' @export
#'
#' @examples
#'
gtfs_isochrones <- function(raptor_result, stops, breaks=c(3600, 2*3600, 3*3600), hull_alpha_min=0.2, buffer_value=5000){


  if(length(unique(raptor_result[, from_stop_id]))>1) stop(" 'raptor_result' has more than one origin stop")

  is_increasing <- function(vec) {
    return(all(diff(vec) > 0))
  }

  ah2sp <- function(x, increment=360, rnd=10, proj4string=CRS(as.character(NA))){
    if (class(x) != "ahull"){
      stop("x needs to be an ahull class object")
    }
    # Extract the edges from the ahull object as a dataframe
    xdf <- as.data.frame(x$arcs)
    # Remove all cases where the coordinates are all the same
    xdf <- subset(xdf,xdf$r > 0)
    res <- NULL
    if (nrow(xdf) > 0){
      # Convert each arc to a line segment
      linesj <- list()
      prevx<-NULL
      prevy<-NULL
      j<-1
      for(i in 1:nrow(xdf)){
        rowi <- xdf[i,]
        v <- c(rowi$v.x, rowi$v.y)
        theta <- rowi$theta
        r <- rowi$r
        cc <- c(rowi$c1, rowi$c2)
        # Arcs need to be redefined as strings of points. Work out the number of points to allocate in this arc segment.
        ipoints <- 2 + round(increment * (rowi$theta / 2),0)
        # Calculate coordinates from arc() description for ipoints along the arc.
        angles <- alphahull::anglesArc(v, theta)
        seqang <- seq(angles[1], angles[2], length = ipoints)
        x <- round(cc[1] + r * cos(seqang),rnd)
        y <- round(cc[2] + r * sin(seqang),rnd)
        # Check for line segments that should be joined up and combine their coordinates
        if (is.null(prevx)){
          prevx<-x
          prevy<-y
        } else if (x[1] == round(prevx[length(prevx)],rnd) && y[1] == round(prevy[length(prevy)],rnd)){
          if (i == nrow(xdf)){
            #We have got to the end of the dataset
            prevx<-append(prevx,x[2:ipoints])
            prevy<-append(prevy,y[2:ipoints])
            prevx[length(prevx)]<-prevx[1]
            prevy[length(prevy)]<-prevy[1]
            coordsj<-cbind(prevx,prevy)
            colnames(coordsj)<-NULL
            # Build as Line and then Lines class
            linej <- Line(coordsj)
            linesj[[j]] <- Lines(linej, ID = as.character(j))
          } else {
            prevx<-append(prevx,x[2:ipoints])
            prevy<-append(prevy,y[2:ipoints])
          }
        } else {
          # We have got to the end of a set of lines, and there are several such sets, so convert the whole of this one to a line segment and reset.
          prevx[length(prevx)]<-prevx[1]
          prevy[length(prevy)]<-prevy[1]
          coordsj<-cbind(prevx,prevy)
          colnames(coordsj)<-NULL
          # Build as Line and then Lines class
          linej <- sp::Line(coordsj)
          linesj[[j]] <- sp::Lines(linej, ID = as.character(j))
          j<-j+1
          prevx<-NULL
          prevy<-NULL
        }
      }
      # Promote to SpatialLines
      lspl <- SpatialLines(linesj)
      # Convert lines to polygons
      # Pull out Lines slot and check which lines have start and end points that are the same
      lns <- slot(lspl, "lines")
      polys <- sapply(lns, function(x) {
        crds <- slot(slot(x, "Lines")[[1]], "coords")
        identical(crds[1, ], crds[nrow(crds), ])
      })
      # Select those that do and convert to SpatialPolygons
      polyssl <- lspl[polys]
      list_of_Lines <- slot(polyssl, "lines")
      sppolys <- SpatialPolygons(list(Polygons(lapply(list_of_Lines, function(x) { Polygon(slot(slot(x, "Lines")[[1]], "coords")) }), ID = "1")), proj4string=proj4string)
      # Create a set of ids in a dataframe, then promote to SpatialPolygonsDataFrame
      hid <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "ID"))
      areas <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "area"))
      df <- data.frame(hid,areas)
      names(df) <- c("HID","Area")
      rownames(df) <- df$HID
      res <- SpatialPolygonsDataFrame(sppolys, data=df)
      res <- res[which(res$Area > 0),]
    }
    return(res)
  }



#Check input
if(!(is.data.table(raptor_result))) stop(" 'raptor_result' must be a data.table")
if(!(is_increasing(breaks))) stop(" 'breaks' must be in increasing order")

  raptor_result <- data.table(raptor_result)
  traveltimes <- data.table(raptor_result[, list(travel_time=min(travel_time)),
                             by = c("from_stop_id", "to_stop_id")])

  stops_sf <- tidytransit::stops_as_sf(stops)
  stops_sf <- merge(stops_sf, traveltimes, by.x="stop_id", by.y="to_stop_id", all.y=T)

  #Translate the cutoffs into categories in the stops_sf table
  stops_sf$cat <- length(breaks)
  for (i in length(breaks):1){
    stops_sf$cat <- ifelse(stops_sf$travel_time<breaks[i], i, stops_sf$cat)
  }

  stops <- merge(stops, stops_sf[, c("stop_id", "cat")], by="stop_id")

#Now calculate the convex hulls for each break
  hulls <- list()
  for(i in 1:(length(breaks))){
    hulls[[i]] <- alphahull::ahull(stops[stops$cat==i, ]$stop_lon, stops[stops$cat==i, ]$stop_lat, alpha = hull_alpha_min+0.5*i*i/length(breaks))
  }





  #Now transform ahull objects to sp objects:
  hull_sp <- lapply(hulls, FUN=ah2sp)

  #And then to sf objects:
  hull_sf <- lapply(hull_sp, FUN=sf::st_as_sf)
  hull_sf <- lapply(hull_sf, FUN=sf::st_set_crs, 4326)

#Now refine the display of each isochrone by unioning the alpha hull with the buffered circles around stops:

  #First generate buffers for the station points: Suggestion: define buffer depending on the scale of the map.

  buffer <- stops_sf

  #Transform to UTM for calculating buffer
  buffer <- sf::st_transform(buffer, 3045)
  hull_sf_buff <- lapply(hull_sf, FUN=sf::st_transform, crs=3045)

  #Create list of unioned isochrones
  buffererd_shapes <- list()
  for(j in 1:(length(breaks)-1)){
    points_buff <- sf::st_buffer(buffer[buffer$cat==j,] , buffer_value)
    points_union <- sf::st_union(points_buff)
    hull_buff <- sf::st_buffer(hull_sf_buff[[j]], buffer_value)

    buffererd_shapes[[j]] <- sf::st_union(points_union, hull_buff)

  }

  #Now take out overlapping
  for(i in (length(buffererd_shapes)):1){
    for(j in (length(buffererd_shapes)):1){
      if(!(i==j)){
        buffererd_shapes[[i]] <- sf::st_difference(buffererd_shapes[[i]], buffererd_shapes[[j]])
      }
    }
  }

  #Convert back to WGS84
  buffererd_shapes <- lapply(buffererd_shapes, FUN=sf::st_transform, crs=4326)


#Last step: convert it all to one single SP object, which can also be saved as one shapefile
  sp <- lapply(buffererd_shapes, as, "Spatial")

  sp_isochrone <- sp[[1]]
  sp_isochrone$break_value <- breaks[1]

  for(i in 2:length(sp)){
    sp_new <- sp[[i]]
    sp_new$break_value <- breaks[i]

    sp_isochrone <- rbind(sp_isochrone, sp_new)
  }

  return(sp_isochrone)

}
