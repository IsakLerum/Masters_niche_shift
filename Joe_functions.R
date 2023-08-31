#' @title Produce species distribution models using INLA with GMRF spatial random effects
#' 
#' @description 
#' This function performs a generalised linear regression model with a spatial SPDE random effect using INLA on
#' provided occurrence data and a set of linear covariates.
#' 
#' @param occurrenceData A \code{\link[sp::SpatialGridDataFrame]{SpatialGridDataFrame}} containing the gridded
#' occurrence data.  All values in the associated \code{data.frame} greater than zero will be treated as a one
#' for the purposes of the model.  Each column in the data frame represents a different species to be modelled.
#' @param climateData A \code{\link[sp::SpatialGridDataFrame]{SpatialGridDataFrame}} containing the climate data.
#' @param alpha A \code{numeric} scalar containing the SPDE fractional operater order.
#' @param meshParameters A \code{list} containing the parameters to use in the construction of the spatial mesh
#' to approximate the random effects.  These parameters are passed to the \code{\link[INLA::inla.mesh.2d]{inla.mesh.2d}}
#' and \code{\link[INLA::inla.nonconvex.hull]{inla.nonconvex.hull}} functions.
#' @param meshBoundary A \code{\link[sp::SpatialPolygons]{SpatialPolygons}} object containing the boundary of the
#' mesh.  If this is \code{NULL} then instead construct a border polygon using the \code{\link[INLA::inla.nonconvex.hull]{inla.nonconvex.hull}}
#' function and parameters passed to \code{meshParameters}.
#' @param progressUpdater A function created using \code{\link[createProgressUpdater]{createProgressUpdater}}
#' that can be used to update a progress indicator.  \code{NULL} means that no progress indication
#' will be provided.
#' @param responseDensity An \code{integer} scalar.  The number of sampling locations to use when building the
#' response curves.
#' @param outFolder A character scalar denoting the folder to place the SDM output.
#' @param createGeoTIFF A logical scalar.  If \code{TRUE} then a series of GeoTIFF files are created
#' in the \code{outFolder} directory with the model predictions for each species.
#' @param createRObFile A logical scalar.  If \code{TRUE} then an R object file (with extension .rds) is
#' created in the \code{outFolder} with the fitted model object for each species.
#' @param inlaVerbose A logical scalar.  If \code{TRUE}, \code{\link[INLA::inla]{inla}} is run in verbose mode.
#' @param inlaKeep A logical scalar.  If \code{TRUE}, \code{\link[INLA::inla]{inla}} retains working files.  These are
#' stored in \code{outFolder}.
#' @param inlaDebug A logical scalar.  If \code{TRUE}, \code{\link[INLA::inla]{inla}} produces debugging information.
#'
#' @return A \code{list} containing two elements:
#' \itemize{
#'  \item{modelSummaries}{A list of \code{\link[INLA::inla]{inla}} model objects for each species in the \code{occurrenceData}
#'  \code{data.frame}.}
#'  \item{spatialMesh}{A \code{\link[sp::SpatialPolygons]{SpatialPolygons}} object containing the spatial mesh used in the
#'  \code{\link[INLA::inla]{inla}} function.}
#' }
#' 
#' @examples
#' @author Joseph D. Chipperfield, \email{joseph.chipperfield@@nmbu.no}
#' @seealso \code{\link[sp::SpatialGridDataFrame]{SpatialGridDataFrame}} \code{\link[INLA::inla.mesh.2d]{inla.mesh.2d}}
#' \code{\link[INLA::inla.nonconvex.hull]{inla.nonconvex.hull}}
#' @export
#'
sdmINLASpatRandEff <- function(occurrenceData, climateData, alpha = 1.5, meshParameters = list(), meshBoundary = NULL, progressUpdater = NULL, responseDensity = 100,
                               outFolder = tempdir(), createGeoTIFF = TRUE, createRObFile = TRUE, inlaVerbose = FALSE, inlaKeep = FALSE, inlaDebug = FALSE, myName = "myModel") {
  ### 1.1 ==== Sanity check the function inputs ====
  updateProgressBar(progressUpdater, value = 0.0, details = "processing inputs")
  # Sanity check the output generation flags
  inCreateGeoTIFF <- tryCatch(as.logical(createGeoTIFF), error = function(err) {
    stop(paste("invalid entry for the GeoTIFF creation flag:", err, sep = " "))
  })
  if(length(inCreateGeoTIFF) <= 0) {
    stop("invalid entry for the GeoTIFF creation flag: zero vector length")
  } else if(length(inCreateGeoTIFF) > 1) {
    warning("GeoTIFF creation flag vector length greater than one: only the first element will be used")
    inCreateGeoTIFF <- inCreateGeoTIFF[1]
  }
  if(any(is.na(inCreateGeoTIFF))) {
    stop("invalid entry for the GeoTIFF creation flag: NA values present")
  }
  inCreateRObFile <- tryCatch(as.logical(createRObFile), error = function(err) {
    stop(paste("invalid entry for the R object file creation flag:", err, sep = " "))
  })
  if(length(inCreateRObFile) <= 0) {
    stop("invalid entry for the R object file creation flag: zero vector length")
  } else if(length(inCreateRObFile) > 1) {
    warning("R object creation flag vector length greater than one: only the first element will be used")
    inCreateRObFile <- inCreateRObFile[1]
  }
  if(any(is.na(inCreateRObFile))) {
    stop("invalid entry for the R object creation flag: NA values present")
  }
  # Sanity check the location of the output folder
  inOutFolder <- tryCatch(as.character(outFolder), error = function(err) {
    stop(paste("invalid entry for the output folder location:", err, sep = " "))
  })
  if(length(inOutFolder) <= 0) {
    stop("invalid entry for the output folder location: zero vector length")
  } else if(length(inOutFolder) > 1) {
    warning("length of vector specifying the location of the ouput folder is greater than one: only the first element will be used")
    inOutFolder <- inOutFolder[1]
  }
  # Check that the folder exists
  if(!dir.exists(inOutFolder)) {
    stop("output directory does not exist")
  }
  # Check the occurrence data
  inOccurrenceData <- tryCatch(as(occurrenceData, "SpatialGridDataFrame"), error = function(err) {
    stop(paste("invalid input for the occurrence data:", err, sep = " "))
  })
  # Check the climate data
  inClimateData <- tryCatch(as(climateData, "SpatialGridDataFrame"), error = function(err) {
    stop(paste("invalid input for the environmental covariates:", err, sep = " "))
  })
  # Ensure that they have the same grid topology and coordinate reference system
  if(!identical(as(inOccurrenceData, "SpatialGrid"), as(inClimateData, "SpatialGrid"))) {
    stop("occurrence data and environmental covariate data do not share the same grid topology and/or coordinate projection system")
  }
  # Test the value of the fractional operator order
  inAlpha <- tryCatch(as.double(alpha), error = function(err) {
    stop(paste("invalid input for the fractional operator order:", err, sep = " "))
  })
  if(is.null(inAlpha) || length(inAlpha) <= 0) {
    stop("invalid input for the fractional operator order: vector has zero length")
  } else if(length(inAlpha) > 1) {
    warning("fractional operator order specification vector length greater than one: only the first element will be used")
    inAlpha <- inAlpha[1]
  }
  if(anyNA(inAlpha)) {
    stop("invalid input for the fractional operator order: NA values detected")
  }
  # Check the mesh boundary parameter
  inMeshBoundary <- meshBoundary
  if(!is.null(inMeshBoundary)) {
    inMeshBoundary <- tryCatch(as(meshBoundary, "SpatialPolygons"), error = function(err) {
      stop(paste("invalid input for the mesh boundary:", err, sep = " "))
    })
  }
  # Check the parameters to generate the spatial mesh
  inMeshParameters <- tryCatch(as.list(meshParameters), error = function(err) {
    stop(paste("invalid input for the mesh parameters:", err, sep = " "))
  })
  # Set the mesh parameters from the input list and use default values if none supplied
  inOffset <- meshParameters$offset
  inN <- meshParameters$n
  inMaxEdge <- meshParameters$max.edge
  inMinAngle <- meshParameters$min.angle
  inCutoff <- meshParameters$cutoff
  inMaxNStrict <- meshParameters$max.n.strict
  inMaxN <- meshParameters$max.n
  inConvex <- meshParameters$convex
  if(is.null(inConvex)) {
    inConvex <- -0.15
  }
  inConcave <- meshParameters$concave
  if(is.null(inConcave)) {
    inConcave <- inConvex
  }
  inResolution <- meshParameters$resolution
  if(is.null(inResolution)) {
    inResolution <- 40
  }
  inEps <- meshParameters$eps
  # Test the response density input
  inResponseDensity <- tryCatch(as.integer(responseDensity), error = function(err) {
    stop(paste("invalid input for the response curve density value:", err, sep = " "))
  })
  if(is.null(inResponseDensity) || length(inResponseDensity) <= 0) {
    stop("invalid input for the response curve density value: specification vector length less than one")
  } else if(length(inResponseDensity) > 1) {
    warning("response curve density specification vector length greater than one: only the first element will be used")
    inResponseDensity <- inResponseDensity[1]
  }
  if(is.na(inResponseDensity) || inResponseDensity < 2) {
    stop("invalid input for the response curve density value: density values NA or less than two")
  }
  ### 1.2 ==== Generate the spatial random effects mesh ====
  updateProgressBar(progressUpdater, value = 0.2, details = "generating mesh for spatial random effect")
  # Determine which rows to use
  useRow <- !apply(X = as.matrix(inClimateData@data), FUN = anyNA, MARGIN = 1)
  boundHull <- NULL
  # Create the bounding hull of the mesh
  if(!is.null(inMeshBoundary)) {
    if(proj4string(inMeshBoundary) != proj4string(inClimateData)) {
      # Transform the mesh boundary if it is on a different projection than the climate data
      inMeshBoundary <- spTransform(inMeshBoundary, CRS(proj4string(inClimateData)))
    }
    boundHull <- inla.sp2segment(inMeshBoundary)
  } else {
    # Create a bounding hull from the climate data coordinates
    boundHull <- inla.nonconvex.hull(points = coordinates(inClimateData)[useRow, ], convex = inConvex, concave = inConcave, resolution = inResolution, eps = inEps, crs = CRS(proj4string(inClimateData)))
  }
  # Create a mesh to define the spatial random effects on
  effectsMesh <- inla.mesh.2d(boundary = boundHull, n = inN, offset = inOffset, max.edge = inMaxEdge, min.angle = inMinAngle, cutoff = inCutoff, max.n.strict = inMaxNStrict, max.n = inMaxN, crs = CRS(proj4string(inClimateData)))
  # Create a projection matrix associated with the mesh
  effectsProjMat <- inla.spde.make.A(mesh = effectsMesh, loc = coordinates(inClimateData)[useRow, ])
  ### 1.3 ==== Formulate the model specification ====
  updateProgressBar(progressUpdater, value = 0.3, details = "processing covariates")
  expandedCovarNames <- c(colnames(inClimateData@data), paste(colnames(inClimateData@data), "quadratic", sep = "_"))
  # Build the SPDE model
  spdeModel <- inla.spde2.matern(mesh = effectsMesh, alpha = inAlpha)
  # Standardise the climate data (important when doing ridge regression)
  meanClimate <- apply(X = as.matrix(inClimateData@data)[useRow, , drop = FALSE], FUN = mean, MARGIN = 2, na.rm = TRUE)
  sdClimate <- apply(X = as.matrix(inClimateData@data)[useRow, , drop = FALSE], FUN = sd, MARGIN = 2, na.rm = TRUE)
  standardClimateData <- sapply(X = 1:ncol(inClimateData@data), FUN = function(curIndex, climateData, meanClimate, sdClimate) {
    (climateData[, curIndex] - meanClimate[curIndex]) / sdClimate[curIndex]
  }, climateData = as.matrix(inClimateData@data), meanClimate = meanClimate, sdClimate = sdClimate)
  colnames(standardClimateData) <- colnames(inClimateData@data)
  # Expand the climate data by adding on the quadratic terms
  expandClimateData <- cbind(as.matrix(standardClimateData), as.matrix(standardClimateData) * as.matrix(standardClimateData))
  colnames(expandClimateData) <- expandedCovarNames
  # Create the model formula (right hand side)
  modelFormRHS <- "-1 + intercept + f(spatIndeces, model = spdeModel) + f(climIndeces, model = \"z\", Z = climCovars)"
  ### 1.4 ==== Create a set of climate data for response curve ====
  inResponseData <- do.call(rbind, lapply(X = colnames(standardClimateData), FUN = function(curColName, climateData, responseDensity) {
    # Initialise a matrix of mean values for each of the columns
    outMat <- matrix(0.0, ncol = ncol(climateData), nrow = responseDensity)
    colnames(outMat) <- colnames(climateData)
    # Calculate the range of values in the current column
    colRange <- range(climateData[, curColName], na.rm = TRUE)
    # Replace the current column with a sequence values between the known range
    outMat[, curColName] <- seq(colRange[1], colRange[2], length.out = responseDensity)
    outMat
  }, climateData = as.matrix(standardClimateData)[useRow, , drop = FALSE], responseDensity = inResponseDensity))
  expandResponseData <- cbind(inResponseData, inResponseData * inResponseData)
  colnames(expandResponseData) <- expandedCovarNames
  ### 1.5 ==== Run the model for each of the species ====
  # Create a list of model object for each species
  outResults <- lapply(X = 1:ncol(occurrenceData), FUN = function(curIndex, inOccurrenceData, expandResponseData,
                                                                  expandClimateData, useRow, spdeModel, effectsProjMat, modelFormRHS, progressUpdater, climSummaryStats, responseDensity,
                                                                  outFolder, createGeoTIFF, createRObFile, inlaVerbose, inlaKeep, inlaDebug) {
    # Retrieve the current species name
    curSpecies <- colnames(inOccurrenceData@data)[curIndex]
    updateProgressBar(progressUpdater, value = 0.4 + 0.6 * curIndex / ncol(inOccurrenceData), detail = paste("fitting distribution model for species", curSpecies, sep = " "))
    ## 1.5.1 Create a data stack ----
    # Retrieve the response data
    respData <- list(c(ifelse(inOccurrenceData@data[useRow, curSpecies] <= 0.0, 0, 1), rep(NA, nrow(expandResponseData))))
    names(respData) <- paste("species", curSpecies, sep = "_")
    # Retrieve the number of non-response rows and response rows
    numNonResponse <- nrow(expandClimateData[useRow, ])
    numResponse <- nrow(expandResponseData)
    numTotal <- numNonResponse + numResponse
    # Build the data stack for the model input
    dataStack <- inla.stack(tag = paste(curSpecies, "Data", sep = ""),
                            data = respData,
                            A = list(rbind(effectsProjMat, matrix(NA, ncol = ncol(effectsProjMat), nrow = numResponse)), 1, 1),
                            effects = list(
                              spatIndeces = 1:spdeModel$n.spde,
                              climIndeces = 1:numTotal,
                              intercept = rep(1.0, numTotal)
                            ))
    # Set the climate covariates
    climCovars <- as.matrix(rbind(as.data.frame(expandClimateData[useRow, ]), expandResponseData))
    ## 1.5.2 Create the full model specification ----
    fullModelForm <- as.formula(paste(names(respData), "~", modelFormRHS, sep = " "))
    ## 1.5.3 Run INLA to create output model ----
    outModel <- inla(fullModelForm, data = inla.stack.data(dataStack), family = "betabinomial",
                     control.predictor = list(compute = TRUE, A = inla.stack.A(dataStack)), verbose = inlaVerbose, keep = inlaKeep,
                     working.directory = paste(outFolder, "/Species_", curSpecies, "_inlaModelFiles", sep = ""), debug = inlaDebug)
    ## 1.5.4 Retrieve the predictions ----
    # Get the linear predictor (after the permutation matrix is applied)
    isPermPred <- grepl("^A[Pp]redictor\\.", rownames(outModel$summary.linear.predictor), perl = TRUE)
    linPermPred <- outModel$summary.linear.predictor[isPermPred, c("0.025quant", "mean", "0.975quant")]
    names(linPermPred) <- c("lowerEst", "mean", "upperEst")
    # Get the spatial random effect (after the permutation matrix is applied)
    spatRandPermPred <- c(as.double(effectsProjMat %*% outModel$summary.random[["spatIndeces"]][, "mean"]), rep(0, responseDensity * ncol(climSummaryStats)))
    # Get the fitted values (after inverse-link transformation)
    #		Define an inverse link function
    invLogit <- function(inValues) {
      exp(inValues) / (1.0 + exp(inValues))
    }
    #		Set the fitted values
    fittedPermPred <- data.frame(
      lowerEst = invLogit(linPermPred$lowerEst),
      mean = invLogit(linPermPred$mean),
      upperEst = invLogit(linPermPred$upperEst))
    # Gather all the predictions together
    allPreds <- data.frame(
      meanEst = fittedPermPred$mean,                                       # The mean estimate (inverse-link transformed)
      lowerEst = fittedPermPred$lowerEst,                                  # The 2.5 percentile estimate (inverse-link transformed)
      upperEst = fittedPermPred$upperEst,                                  # The 97.5 percentile estimate (inverse-link transformed)
      uncertaintyEst = fittedPermPred$upperEst - fittedPermPred$lowerEst,  # The uncertainty range (97.5 percentile - 2.5 percentile)
      meanLinearPred = linPermPred$mean,                                   # The mean estimate (linear predictor)
      meanClimatePred = linPermPred$mean - spatRandPermPred,               # The mean climate component (linear predictor)
      meanSpatialPred = spatRandPermPred)                                  # The mean spatial component (linear predictor)
    ## 1.5.5 Create the prediction geographical objects ----
    predictRaster <- SpatialGridDataFrame(as(inOccurrenceData, "SpatialGrid"), data.frame(
      # Initialise a data frame of NA values
      meanEst = rep(NA, nrow(inOccurrenceData@data)),
      lowerEst = rep(NA, nrow(inOccurrenceData@data)),
      upperEst = rep(NA, nrow(inOccurrenceData@data)),
      uncertaintyEst = rep(NA, nrow(inOccurrenceData@data)),
      meanLinearPred = rep(NA, nrow(inOccurrenceData@data)),
      meanClimatePred = rep(NA, nrow(inOccurrenceData@data)),
      meanSpatialPred = rep(NA, nrow(inOccurrenceData@data))))
    # Fill those cells with predictions
    predictRaster@data[useRow, ] <- allPreds[1:(nrow(expandClimateData[useRow, ])), ]
    ## 1.5.6 Create the response curves ----
    responseCurves <- lapply(X = 1:ncol(climSummaryStats), FUN = function(curIndex, climSummaryStats, inputPreds, responseDensity) {
      # Calculate the current response indeces
      respIndeces <- ((curIndex - 1) * responseDensity + 1):(curIndex * responseDensity)
      # Retrieve the current climate variable
      curClimVarName <- colnames(climSummaryStats)[curIndex]
      curClimVar <- inputPreds[respIndeces, curClimVarName] * climSummaryStats["sdClimate", curClimVarName] + climSummaryStats["meanClimate", curClimVarName]
      # Repackage the prediction
      data.frame(
        covarVal = curClimVar,
        meanEst = inputPreds[respIndeces, "meanEst"],
        lowerEst = inputPreds[respIndeces, "lowerEst"],
        upperEst = inputPreds[respIndeces, "upperEst"])
    }, climSummaryStats = climSummaryStats, responseDensity = responseDensity,
    inputPreds = cbind(allPreds[nrow(allPreds) - (responseDensity * ncol(climSummaryStats) - 1):0, ], expandResponseData))
    #	inputPreds = cbind(allPreds, rbind(as.data.frame(expandClimateData[useRow, ]), as.data.frame(expandResponseData)))[nrow(allPreds) - (responseDensity * ncol(climSummaryStats) - 1):0, ])
    # Set the names of the response curves object
    names(responseCurves) <- colnames(climSummaryStats)
    ## 1.5.7 Create the output files ----
    if(createGeoTIFF) {
      # If GeoTIFFs are requested then produce them
      writeGDAL(predictRaster["meanEst"], paste(outFolder, "/Species_", curSpecies, "_MeanEst.tif", sep = ""), drivername = "GTiff")
      writeGDAL(predictRaster["lowerEst"], paste(outFolder, "/Species_", curSpecies, "_LowerEst.tif", sep = ""), drivername = "GTiff")
      writeGDAL(predictRaster["upperEst"], paste(outFolder, "/Species_", curSpecies, "_UpperEst.tif", sep = ""), drivername = "GTiff")
      writeGDAL(predictRaster["uncertaintyEst"], paste(outFolder, "/Species_", curSpecies, "_UncertaintyEst.tif", sep = ""), drivername = "GTiff")
      writeGDAL(predictRaster["meanLinearPred"], paste(outFolder, "/Species_", curSpecies, "_MeanLinearPred.tif", sep = ""), drivername = "GTiff")
      writeGDAL(predictRaster["meanClimatePred"], paste(outFolder, "/Species_", curSpecies, "_MeanClimatePred.tif", sep = ""), drivername = "GTiff")
      writeGDAL(predictRaster["meanSpatialPred"], paste(outFolder, "/Species_", curSpecies, "_MeanSpatialPred.tif", sep = ""), drivername = "GTiff")
    }
    if(createRObFile) {
      # If an R object file is requested then produce that
      saveRDS(list(
        modelObject = outModel,               # The fitted model object
        spatialPredictions = predictRaster,   # The spatial predictions made by the model
        responsePredictions = responseCurves  # The response curve information
      ), file = paste(outFolder, "/", myName,".rds", sep = ""))
    }
    # Create a summary of the model output
    summary(outModel)
  }, inOccurrenceData = inOccurrenceData, expandResponseData = expandResponseData, expandClimateData = expandClimateData, useRow = useRow,
  spdeModel = spdeModel, effectsProjMat = effectsProjMat, modelFormRHS = modelFormRHS, progressUpdater = progressUpdater,
  climSummaryStats = rbind(meanClimate, sdClimate), responseDensity = inResponseDensity, outFolder = inOutFolder,
  createGeoTIFF = inCreateGeoTIFF, createRObFile = inCreateRObFile, inlaVerbose = inlaVerbose, inlaKeep = inlaKeep, inlaDebug = inlaDebug)
  names(outResults) <- colnames(inOccurrenceData@data)
  list(modelSummaries = outResults, spatialMesh = inlaMeshToSpatialPolygons(effectsMesh))
}

## 1. ------ FUNCTION TO CONVERT INLA MESH OBJECTS TO SPATIAL POLYGONS ------
#' @title Converts the Delaunay Triangulation in INLA Mesh Objects to SpatialPolygons
#' 
#' @description 
#' Converts the constructed Delaunay triangulation from a two-dimensional INLA \link[INLA::inla.mesh.create]{mesh} object
#' into a \code{\link[sp::SpatialPolygons]{SpatialPolygons}} object.  This can be useful if you want to save the mesh
#' in standard vector file formats or perform further processing on the mesh.
#' 
#' @param meshInput An \code{\link[INLA::inla.mesh.create]{inla.mesh}} object.  This mesh must be two-dimensional.
#' 
#' @return A \code{\link[sp::SpatialPolygons]{SpatialPolygons}} object containing the Delaunay triangulation from
#' the input mesh.  Note that inner and outer boundary information is not present in the output.
#' 
#' @examples
#' @author Joseph D. Chipperfield, \email{joseph.chipperfield@@nmbu.no}
#' @seealso \code{\link[sp::SpatialPolygons]{SpatialPolygons}} \code{\link[sp::SpatialPolygons]{SpatialPolygons}}
#' @export
#'
inlaMeshToSpatialPolygons <- function(meshInput) {
  # Sanity test the input
  inMeshInput <- tryCatch(as(meshInput, "inla.mesh"), error = function(err) {
    stop(paste("invalid input mesh:", err, sep = " "))
  })
  # Ensure the mesh has the correct manifold type
  if(inMeshInput$manifold != "R2") {
    stop("invalid input mesh: mesh must be two-dimensional")
  }
  # Retrieve the Delaunay triangles as polygon objects
  delaunayPolygons <- lapply(X = 1:nrow(inMeshInput$graph$tv), FUN = function(curRowIndex, pointCoords, pointIndeces) {
    # Retrieve the vertex coordinates of the current triangle
    curIndeces <- pointIndeces[curRowIndex, ]
    # Construct a Polygons object to contain the triangle
    Polygons(list(Polygon(
      pointCoords[c(curIndeces, curIndeces[1]), ], hole = FALSE
    )), ID = paste("Delaunay", curRowIndex, sep = "_"))
  }, pointCoords = inMeshInput$loc[, c(1, 2)], pointIndeces = inMeshInput$graph$tv)
  # Convert the Polygons into a SpatialPolygons object
  SpatialPolygons(delaunayPolygons, proj4string = inMeshInput$crs)
}

## 1. ------ DEFINE A FUNCTION TO CREATE A PROGRESS BAR UPDATER FUNCTION ------
#' @title Create a Function to Update a Progress Bar
#' 
#' @description 
#' This function creates another function that can be used to update a progress
#' bar.  Useful to pass to functions that may take a long time to run and needs
#' to give user feedback on the progress.
#' 
#' @param progressOb An object representing a progress indicators.  Either this
#' can be the output from the \code{\link[utils::txtProgressBar]{txtProgressBar}}
#' function or from the \code{\link[shiny::Progress]{Progress}} function.
#' 
#' @return A function that can be called to update the progress object given in
#' \code{progressOb}.  The function can take the following arguments:
#' \describe{
#'     \item{\code{message}}{A message to display to the user about the currently
#'     executing process}
#'     \item{\code{details}}{A message to display to the user about the details
#'     of the currently executing process}
#'     \item{\code{value}}{A value between the minimum and maximum values set in
#'     the progress object}
#' }
#' 
#' @examples
#' # Initialise a progress bar to use
#' examProgress <- txtProgressBar()
#' # Create a function to update the progress bar
#' examProgressUpdate <- createProgressUpdater(examProgress)
#' # Use the function to update the progress bar
#' examProgressUpdate(
#'     message = "This is a test progress update",
#'     details = "Lets say it is 50% done",
#'     value = 0.5)
#' @author Joseph D. Chipperfield, \email{joseph.chipperfield@@nmbu.no}
#' @seealso \code{\link[utils::txtProgressBar]{txtProgressBar}}
#' \code{\link[shiny::Progress]{Progress}}
#' @export
#' 
createProgressUpdater <- function(progressOb) {
  ### 1.1 ==== Sanity test the inputs ====
  # Determine the type of progress bar initialised
  if(class(progressOb) == "txtProgressBar") {
    inProgressType <- "utils"
  } else if(all(class(progressOb) == c("Progress", "R6"))) {
    inProgressType <- "shiny"
  }
  outFunc <- NULL
  ### 1.2 ==== Create an updater function for shiny progress bars ====
  if(as.character(inProgressType) == "shiny") {
    outFunc <- function(message, details, value, valProgressOb = progressOb) {
      ## 1.2.1 Sanity test the inputs ----
      # Retrieve the updater message
      inMessage <- message
      if(!is.null(inMessage)) {
        # Convert the message to a character vector
        inMessage <- tryCatch(as.character(inMessage), error = function(err) {
          c("")
        })
        # Swap out any NAs to ""
        inMessage[is.na(inMessage)] <- ""
        # Collapse the message into one element
        inMessage <- paste(inMessage, collapse = " ")
      }
      # Retrieve the updater details
      inDetails <- details
      if(!is.null(inDetails)) {
        # Convert the details to a character vector
        inDetails <- tryCatch(as.character(inDetails), error = function(err) {
          c("")
        })
        # Swap out any NAs to ""
        inDetails[is.na(inDetails)] <- ""
        # Collapse the message into one element
        inDetails <- paste(inDetails, collapse = " ")
      }
      # Retrieve the updater value
      inValue <- value
      if(!is.null(inValue)) {
        # Convert the value to a numeric vector
        inValue <- tryCatch(as.numeric(inValue), error = function(err) {
          c()
        })
        # Set the input value to NULL if the input is NA, NaN, non-finite or length 0
        if(length(inValue) <= 0 || is.na(inValue[1]) || !is.finite(inValue[1])) {
          inValue <- NULL
        } else {
          inValue <- inValue[1]
        }
      }
      ## 1.2.2 Update the progress object ----
      # Update the message if requested
      tryCatch(valProgressOb$set(message = inMessage, detail = inDetails, value = inValue), error = function(err) {
        stop(paste("unable to update progress indicator:", err, sep = " "))	
      })
    }
    ### 1.3 ==== Create an updater function for utils progress bars ====	
  } else if(as.character(inProgressType) == "utils") {
    outFunc <- function(message, details, value, valProgressOb = progressOb) {
      ## 1.3.1 Sanity test the inputs ----
      # Retrieve the updater message
      inMessage <- message
      if(!is.null(inMessage)) {
        # Convert the message to a character vector
        inMessage <- tryCatch(as.character(inMessage), error = function(err) {
          c("")
        })
        # Swap out any NAs to ""
        inMessage[is.na(inMessage)] <- ""
        # Collapse the message into one element
        inMessage <- paste(inMessage, collapse = " ")
      }
      # Retrieve the updater details
      inDetails <- details
      if(!is.null(inDetails)) {
        # Convert the details to a character vector
        inDetails <- tryCatch(as.character(inDetails), error = function(err) {
          c("")
        })
        # Swap out any NAs to ""
        inDetails[is.na(inDetails)] <- ""
        # Collapse the message into one element
        inDetails <- paste(inDetails, collapse = " ")
      }
      # Retrieve the updater value
      inValue <- value
      if(!is.null(inValue)) {
        # Convert the value to a numeric vector
        inValue <- tryCatch(as.numeric(inValue), error = function(err) {
          c()
        })
        # Set the input value to NULL if the input is NA, NaN, non-finite or length 0
        if(length(inValue) <= 0 || is.na(inValue[1]) || !is.finite(inValue[1])) {
          inValue <- NULL
        } else {
          inValue <- inValue[1]
        }
      }
      ## 1.3.2 Update the progress object ----
      if(is.null(inValue)) {
        inValue <- NA
      }
      tryCatch(setTxtProgressBar(valProgressOb, value = inValue, title = inMessage, label = inDetails), error = function(err) {
        stop(paste("unable to update progress indicator:", err, sep = " "))
      })
    }
  } else {
    stop("unknown progress indicator object type")
  }
  outFunc
}

## 2. ------ DEFINE A FUNCTION TO UPDATE THE PROGRESS BAR ------
#' @title Update Progress Bar
#' 
#' @description 
#' This function calls a progress update function and passes it the inputs
#' provided as parameters.  This functions is largely only used as a helper
#' functions inside functions that take the output from
#' \code{\link[createProgressUpdater]{createProgressUpdater}} as an argument-
#' 
#' @param updateFunc A function created using
#' \code{\link[createProgressUpdater]{createProgressUpdater}}.
#' @param value Numeric scalar.  A value between the minimum and maximum
#' values set in the progress object.
#' @param details Character scalar. A message to display to the user about
#' the details of the currently executing process.
#' @param message Character scalar. A message to display to the user about
#' the currently executing process.
#' 
#' @return For progress bars created by \code{\link[utils::txtProgressBar]{txtProgressBar}},
#' the function returns the output from \code{\link[utils::txtProgressBar]{setProgressBar}}.
#' For progress bars created by \code{\link[shiny::Progress]{Progress}}, the
#' function returns the output from the \code{set} method of the respective
#' \code{\link[shiny::Progress]{Progress}} object.
#' 
#' @examples
#' # Initialise a progress bar to use
#' examProgress <- txtProgressBar()
#' # Create a function to update the progress bar
#' examProgressUpdate <- createProgressUpdater(examProgress)
#' # Use the function to update the progress bar
#' updateProgressBar(
#'     updateFunc = examProgressUpdate,
#'     message = "This is a test progress update",
#'     details = "Lets say it is 50% done",
#'     value = 0.5)
#' @author Joseph D. Chipperfield, \email{joseph.chipperfield@@nmbu.no}
#' @seealso \code{\link[utils::txtProgressBar]{txtProgressBar}}
#' \code{\link[shiny::Progress]{Progress}}
#' @export
#' 
updateProgressBar <- function(updateFunc, value, details = NULL, message = NULL) {
  progressOut <- NULL
  if(!is.null(updateFunc)) {
    # Coerce the updateFunc input to a function
    inUpdateFunc <- tryCatch(as.function(updateFunc), error = function(err) {
      stop(paste("invalid progress bar update function:", err, sep = " "))
    })
    # Use the updateFunc to update the progress bar
    progressOut <- tryCatch(inUpdateFunc(value = value, details = details, message = message), error = function(err) {
      stop(paste("error encountered during application of progress bar update function:", err, sep = " "))
    })
  }
  progressOut
}

## 1. ------ DEFINE A FUNCTION TO PREDICT OCCURRENCE BASED ON FUNDAMENTAL NICHE ------
#' @title Predics species occurrence based on fundamental niche
#' 
#' @description 
#' This function predicts species occurrence based on fundamental niche (climate only prediction). Can be used 
#' to extrapolate predictions to new areas. [Potential issue: Might be underpredicting occurrence in the case of strong random spatial effect]
#' 
#' @param modelRDS RDS file produced by sdmINLASpatRandEff function containing modelObject, spatialPredictions, 
#' and responsePredictions.
#' @param climateData A \code{\link[sp::SpatialGridDataFrame]{SpatialGridDataFrame}} containing the climate data
#' for the target area. Is assumed to contain the same variables (and in the same order) as the climateData used 
#' for model fit.
#' @param origClimateData A \code{\link[sp::SpatialGridDataFrame]{SpatialGridDataFrame}} containing the climate 
#' data used for model (sdmINLASpatRandEff) fit. 
#'
#' @return A \code{\link[sp::SpatialGridDataFrame]{SpatialGridDataFrame}} containing probability of occurrence.
#' 
#' @examples
#' @author Adam Klimes
#' @seealso \code{\link[sp::SpatialGridDataFrame]{SpatialGridDataFrame}} 
#' @export
#'
predictSD <- function(modRDS, climateData, origClimateData = NULL){
  if (any(names(modRDS$responsePredictions) != colnames(climateData@data)))
    warning("Names of climate variables do not match between modRDS and climateData!")
  if (is.null(origClimateData)){
    origClimateData <- climateData
    warning("origClimateData not provided. Using climateData (assuming predictions are done for original area)!")
  }
  invLogit <- function(x) exp(x) / (1 + exp(x))
  climSummary <- data.frame(lapply(origClimateData@data, function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))))
  pars <- c(modRDS$modelObject$summary.fixed$mean, 
            tail(modRDS$modelObject$summary.random$climIndeces$mean, length(climSummary) * 2))
  climSt <- data.frame(Map(function(x, y) (x - y[1]) / y[2], climateData@data, climSummary))
  climStExt <- cbind(Int = 1, climSt, setNames(climSt * climSt, paste0(names(climSt), "_quadratic")))  
  pred <- data.frame(OccurrenceProb = invLogit(as.matrix(climStExt) %*% pars))
  predGrid <- climateData[1]
  predGrid@data <- pred
  predGrid
}


#get parameters of response curves from the model
getPars <- function(modRDS, climateData){
  pars <- c(modRDS$modelObject$summary.fixed$mean,
            tail(modRDS$modelObject$summary.random$climIndeces$mean, length(climateData@data) * 2))
  climNames <- names(climateData@data)
  names(pars) <- c("Intercept", climNames, paste0(climNames, "_quadratic"))
  pars
}


# Calculates the r^2 of the model, the spatial effect and environmental effects
calcVariancePartitioning <- function(inObject) {
  # Inverse logit function
  ilogit <- function(inVec) {
    1.0 / (1.0 + exp(-inVec))
  }
  # Function to calculate the R squared
  #calcRSq <- function(dataVec, meanPredVec) {
  #  dataMean <- mean(dataVec, na.rm = TRUE)
  #  SStotal <- sum((dataVec - dataMean)^2.0, na.rm = TRUE)
  #  SSres <- sum((dataVec - meanPredVec)^2.0, na.rm = TRUE)
  #  1.0 - SSres / SStotal
  #}
  # Retrieve the response variable (occurrence values)
  occurrenceValues <- inObject$modelObject$.args$data[[
    as.character(inObject$modelObject$.args$formula)[2]
  ]]
  # Retrieve the permutation matrix
  permMatrix <- inObject$modelObject$.args$control.predictor$A
  # Retrieve the spatial component of the prediction
  inSpatIndexes <- inObject$modelObject$.args$data$spatIndeces
  useSpatVal <- !is.na(inSpatIndexes)
  # Retrieve the spatial effect values at each of the data points
  spatValues <- as.double(permMatrix[, useSpatVal] %*% inObject$modelObject$summary.random$spatIndeces[inSpatIndexes[useSpatVal], "mean"]) +
    inObject$modelObject$summary.fixed["intercept", "mean"]
  # Retrieve the full prediction
  fullValues <- inObject$modelObject$summary.linear.predictor[
    grepl("^A", row.names(inObject$modelObject$summary.linear.predictor), perl = TRUE), "mean"]
  # Retrieve the climate-only prediction
  climValues <- fullValues - spatValues + inObject$modelObject$summary.fixed["intercept", "mean"]
  # Calculate the sum of squares of the data, spatial, climate and full models
  isPredData <- !is.na(occurrenceValues)
  ssData <- sum((occurrenceValues[isPredData] - mean(occurrenceValues[isPredData]))^2.0)
  ssSpat <- sum((occurrenceValues[isPredData] - ilogit(spatValues[isPredData]))^2.0)
  ssClim <- sum((occurrenceValues[isPredData] - ilogit(climValues[isPredData]))^2.0)
  ssFull <- sum((occurrenceValues[isPredData] - ilogit(fullValues[isPredData]))^2.0)
  rSqSpat <- 1.0 - ssSpat / ssData
  rSqClim <- 1.0 - ssClim / ssData
  rSqFull <- 1.0 - ssFull / ssData
  # Calculate the mean magnitude of the values
  meanMagSpat <- mean(abs(spatValues[isPredData] - inObject$modelObject$summary.fixed["intercept", "mean"]))
  meanMagClim <- mean(abs(climValues[isPredData] - inObject$modelObject$summary.fixed["intercept", "mean"]))
  meanMagSpatClim <- mean(abs(spatValues[isPredData] - inObject$modelObject$summary.fixed["intercept", "mean"]) / abs(climValues[isPredData] - inObject$modelObject$summary.fixed["intercept", "mean"]))
  list(
    rSquared = setNames(c(rSqSpat, rSqClim, rSqFull), c("spatial", "climate", "full")),
    rSquaredRatio = setNames(c(rSqSpat, rSqClim) / rSqFull, c("spatial", "climate")),
    meanMagnitude = setNames(c(meanMagSpat, meanMagClim, meanMagSpatClim), c("spatial", "climate", "spatial/climate"))
  )
}