#' fit
#'
#' @description Function to train an Artificial Hydrocarbon Network (AHN).
#'
#' @param Sigma a list with two data frames. One for the inputs X, and one for the outputs Y.
#' @param n number of particles to use.
#' @param eta learning rate of the algorithm. Default is \code{0.01}.
#' @param maxIter maximum number of iterations.
#'
#' @return an object of class "\code{ahn}" with the following components:
#' \itemize{
#'         \item network: structure of the AHN trained.
#'         \item Yo: original output variable.
#'         \item Ym: predicted output variable.
#'         \item eta: learning rate.
#'         \item minOverallError: minimum error achieved.
#'         \item variableNames: names of the input variables.
#' }
#' @export
#'
#' @examples
#' # Create data
#' x <- 2 * runif(1000) - 1;
#' x <- sort(x)
#'
#' y <- (x < 0.1) * (0.05 * runif(100) + atan(pi*x)) +
#'     (x >= 0.1 & x < 0.6) * (0.05 * runif(1000) + sin(pi*x)) +
#'     (x >= 0.6) * (0.05 * runif(1000) + cos(pi*x))
#'
#' # Create Sigma list
#' Sigma <- list(X = data.frame(x = x), Y = data.frame(y = y))
#'
#' # Train AHN
#' ahn <- fit(Sigma, 5, 0.01, 500)
#'
fit <- function(Sigma, n, eta, maxIter = 2000) {
    # Security Checking
    if (length(Sigma) != 2) {
        stop("Sigma must be a list with 2 data frames: one for the predictors and one for the outcome variables. ", call. = FALSE)
    }

    if (!is.data.frame(Sigma$X)) {
        stop("The first component of Sigma must be a data frame with the predictor variables. ", call. = FALSE)
    }

    if (!is.data.frame(Sigma$Y)) {
        stop("The second component of Sigma must be a data frame with the outcome variables. ", call. = FALSE)
    }

    if (n < 1) {
        stop("At least two particles are required.", call. = FALSE)
    }

    if ((eta <= 0) | (eta >= 1)) {
        stop("The learning rate eta must be between 0 and 1 (exclusive).", call. = FALSE)
    }

    if (maxIter < 1) {
        stop("Maximum number of iterations must be a positive integer.", call. = FALSE)
    }

    maxIter <- maxIter
    iter <- 1

    X <- Sigma$X
    Y <- Sigma$Y

    newIdx <- sort(Y[ , 1], index.return = TRUE)$ix
    X <- X[newIdx, , drop = FALSE]
    Y <- Y[newIdx, , drop = FALSE]

    numOutputs <- ncol(Y)

    overallError <- rep(Inf, maxIter)
    minOverallError <- Inf

    # Create a saturated and linear compound
    C <- CreateLinearCompound(n)

    # Initialize randomly the position of molecules
    posMolecules <- ahn_Pi <- InitialPosition(X, n)

    # Initialize AHN-structure
    H.names <- paste("molecule", 1:n, sep = "")
    H <- vector("list", length(H.names))
    names(H) <- H.names
    ahn_H <- H

    # Train compound
    while (iter < maxIter) {
        requireRelocation <- TRUE
        relocationIter <- 1
        maxRI <- 10

        while (requireRelocation & (relocationIter < maxRI)) {
            molecules <- DataInMolecules(Sigma, posMolecules)

            Yo <- Ym <- matrix(0, nrow = nrow(Y), ncol = ncol(Y))
            pointer <- 1

            errorMolecule <- matrix(0, nrow = n+1, ncol = numOutputs)

            for (i in seq_len(length(molecules$SigmaSplit$X))) {
                moleculeUsed <- paste('molecule', molecules$moleculesUsed[i], sep = '')
                Xi <- molecules$SigmaSplit$X[[moleculeUsed]]
                Yi <- molecules$SigmaSplit$Y[[moleculeUsed]]

                parameters <- ComputeMoleculeParameters(Xi, Yi, C$Omega[molecules$moleculesUsed[i]])

                H[[moleculeUsed]] <- parameters$H

                errorMolecule[molecules$moleculesUsed[i], ] <- parameters$errorMolecule

                pointerNew <- pointer + nrow(parameters$Yapprox)
                Yo[pointer:(pointerNew-1), ] <- Yi[seq_len(nrow(Yi)), ]
                Ym[pointer:(pointerNew-1), ] <- parameters$Yapprox
                pointer <- pointerNew
            }

            relocation <- RelocateMolecules(posMolecules, molecules$moleculesUsed, errorMolecule)
            requireRelocation <- relocation$requireRelocation
            posMolecules <- relocation$posMolecules
            relocationIter <- relocationIter + 1
        }

        if (relocationIter >= maxIter) {
            H <- ahn_H
            posMolecules <- ahn_Pi
        }

        # Update AHN-structure
        overallError[iter] <- sum(errorMolecule)

        if (overallError[iter] < minOverallError) {
            ahn_H <- H
            ahn_Pi <- posMolecules
            minOverallError = overallError[iter]
        }

        # Update intermolecular distances
        for (i in seq_len(n)) {
            deltaPosMolecule <- -eta * (sum(errorMolecule[i, ]) - sum(errorMolecule[i+1, ])) * rep(1, ncol(posMolecules))
            posMolecules[i, ] <- posMolecules[i, ] + deltaPosMolecule
        }

        iter <- iter + 1
    }
    network = list(H = ahn_H, Pi = ahn_Pi, n = n, C = C)
    ahn <- list(network = network,
                Yo = Yo,
                Ym = Ym,
                eta = eta,
                minOverallError = minOverallError,
                variableNames = names(Sigma$X))
    class(ahn) <- "ahn"
    ahn
}


#' Checks if argument is a \code{ahn} object
#'
#' @param x An \R object
#'
#' @export
#'
is.ahn <- function(x) inherits(x, "ahn")


# #' plotAHN
# #'
# #' @param ahn a list produced from the \link{fit} function.
# #'
# #' @return visualization of the AHN.
# #' @export
# #'
# #' @examples
# #' # Create data
# #' x <- 2 * runif(1000) - 1;
# #' x <- sort(x)
# #'
# #' y <- (x < 0.1) * (0.05 * runif(100) + atan(pi*x)) +
# #'     (x >= 0.1 & x < 0.6) * (0.05 * runif(1000) + sin(pi*x)) +
# #'     (x >= 0.6) * (0.05 * runif(1000) + cos(pi*x))
# #'
# #' # Create Sigma list
# #' Sigma <- list(X = data.frame(x = x), Y = data.frame(y = y))
# #'
# #' # Train AHN
# #' ahn <- fit(Sigma, 5, 0.01, 500)
# #'
# #' # Plot AHN
# #' plotAHN(ahn)
# #'
# plotAHN <- function(ahn) {
#     vis <- CreateNodesEdges(ahn)
#     graph <- igraph::graph_from_data_frame(vis$edges, directed = FALSE, vertices = vis$nodes)
#     ggraph(graph, layout = 'graphopt') +
#         geom_edge_link(aes(start_cap = label_rect(node1.name),
#                            end_cap = label_rect(node2.name))) +
#         geom_node_text(label = vis$nodes$label) +
#         theme_void()
# }


#' Summary Artificial Hydrocarbon Network
#'
#' @description Summary method for objects of class \code{ahn}.
#'
#' @param object an object of class "\code{ahn}" produced from the \link{fit} function.
#' @param ... further arguments passed to or from other methods.
#'
#' @return summary description of the AHN.
#' @export
#'
#' @examples
#' \dontrun{
#' # Create data
#' x <- 2 * runif(1000) - 1;
#' x <- sort(x)
#'
#' y <- (x < 0.1) * (0.05 * runif(100) + atan(pi*x)) +
#'     (x >= 0.1 & x < 0.6) * (0.05 * runif(1000) + sin(pi*x)) +
#'     (x >= 0.6) * (0.05 * runif(1000) + cos(pi*x))
#'
#' # Create Sigma list
#' Sigma <- list(X = data.frame(x = x), Y = data.frame(y = y))
#'
#' # Train AHN
#' ahn <- fit(Sigma, 5, 0.01, 500)
#'
#' # Summary AHN
#' summary(ahn)
#' }
#'
summary.ahn <- function(object, ...) {
    stopifnot(is.ahn(object))
    ahn <- object

    cat("\nArtificial Hydrocarbon Network trained:\n\n")
    cat("Number of molecules:\n", ahn$network$n, "\n\n")
    cat("Learning factor:\n", ahn$eta, "\n\n")
    cat("Overall error:\n", round(ahn$minOverallError, 4), "\n\n")

    centers <- ahn$network$Pi
    rownames(centers) <- paste('molecule', seq_len(nrow(centers)), sep = "")
    colnames(centers) <- ahn$variableNames
    cat("Centers of the molecules:\n")
    print(as.table(centers))

    cat("\nMolecules:\n")
    CreateTable(ahn)
}


#' predict
#'
#' @description Function to simulate a trained Artificial Hydrocarbon Network.
#'
#' @param object an object of class "\code{ahn}" produced from the \link{fit} function.
#' @param ... further arguments passed to or from other methods.
#'
#' @return predicted output values for inputs in \code{newdata}.
#' @export
#'
#' @examples
#' \dontrun{
#' # Create data
#' x <- 2 * runif(1000) - 1;
#' x <- sort(x)
#'
#' y <- (x < 0.1) * (0.05 * runif(100) + atan(pi*x)) +
#'     (x >= 0.1 & x < 0.6) * (0.05 * runif(1000) + sin(pi*x)) +
#'     (x >= 0.6) * (0.05 * runif(1000) + cos(pi*x))
#'
#' # Create Sigma list
#' Sigma <- list(X = data.frame(x = x), Y = data.frame(y = y))
#'
#' # Train AHN
#' ahn <- fit(Sigma, 5, 0.01, 500)
#'
#' # Test AHN
#' X <- data.frame(x = x)
#' ysim <- predict(ahn, X)
#' }
#'
predict.ahn <- function(object, ...) {
    # Security Checking
    stopifnot(is.ahn(object))
    ahn <- object

    dots <- list(...)
    newdata <- dots[[1]]

    if (!is.data.frame(newdata)) {
        stop("newdata must be a data frame with the predictor variables. ", call. = FALSE)
    }

    # Extract network components
    H <- ahn$network$H
    posMolecules <- ahn$network$Pi
    C <- ahn$network$C

    # Initial statemets
    Yapprox <- matrix(0, nrow = nrow(newdata), ncol = max(unlist(sapply(H, ncol))))
    indexes <- rep(0, nrow(newdata))

    # Distribute data over molecules
    molecules <- SimDataInMolecules(newdata, posMolecules)

    # Evaluate AHN-model
    pointer <- 1

    for (i in seq_len(length(molecules$X))) {
        Xi <- molecules$X[[i]]
        indexesi <- molecules$Index[[i]]

        ki <- C$Omega[i]
        Phi <- CH_X(Xi, ki)

        if (is.null(H[[i]])) {next}

        Ym <- Phi %*% H[[i]]

        pointerNew <- pointer + nrow(Ym)
        Yapprox[pointer:(pointerNew-1), ] <- Ym#[seq_len(nrow(Ym)), ]
        indexes[pointer:(pointerNew-1)] <- indexesi
        pointer <- pointerNew
    }

    # Sort data
    ix <- sort(indexes, index.return = TRUE)$ix;
    Y <- Yapprox[ix, ]
}


#' Visualize Artificial Hydrocarbon Network
#'
#' @description Visualize method for objects of class \code{ahn}.
#'
#' @param x an object of class "\code{ahn}" produced from the \link{fit} function.
#' @param ... further arguments passed to visNetwork functions.
#'
#' @return dynamic visualization of the AHN.
#' @export
#'
#' @examples
#' \dontrun{
#' # Create data
#' x <- 2 * runif(1000) - 1;
#' x <- sort(x)
#'
#' y <- (x < 0.1) * (0.05 * runif(100) + atan(pi*x)) +
#'     (x >= 0.1 & x < 0.6) * (0.05 * runif(1000) + sin(pi*x)) +
#'     (x >= 0.6) * (0.05 * runif(1000) + cos(pi*x))
#'
#' # Create Sigma list
#' Sigma <- list(X = data.frame(x = x), Y = data.frame(y = y))
#'
#' # Train AHN
#' ahn <- fit(Sigma, 5, 0.01, 500)
#'
#' # Visualize AHN
#' visualize(ahn)
#' }
#'
visualize <- function(x, ...) {
    stopifnot(is.ahn(x))
    vis <- CreateNodesEdges(x)
    visNetwork(vis$nodes, vis$edges, width = "100%", ...) %>%
        visGroups(groupname = "C", color = "#fbb4ae", ...) %>%
        visGroups(groupname = "H1", color = "#b3cde3", ...) %>%
        visGroups(groupname = "H2", color = "#ccebc5", ...) %>%
        visGroups(groupname = "H3", color = "#decbe4", ...) %>%
        visLegend(position = "right", main = "Legend", ...)
}
