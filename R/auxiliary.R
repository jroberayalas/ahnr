# CH_X creates the Vandermonde matrix
CH_X <- function(X, h) {
    N <- ncol(X)
    q <- nrow(X)

    Phi <- matrix(0, nrow = q, ncol = N * (h+1))

    for (i in seq_len(N)) {
        #Phi <- cbind(Phi, matrixcalc::vandermonde.matrix(X[ , i], h+1))
        Phi[ , (1:(h+1)) + (i-1)*(h+1)] <- matrixcalc::vandermonde.matrix(X[ , i], h+1)
    }
    #Phi <- Phi[ , -1]

    if (is.vector(Phi)) {
        Phi <- matrix(Phi, nrow = 1)
    }

    Phi
}

# ComputeMoleculeParameters solves the system of equations with Least Squares
ComputeMoleculeParameters <- function(Xi, Yi, Omegai) {
    Phi <- CH_X(Xi, Omegai)

    H <- pracma::mldivide(Phi, as.matrix(Yi))
    Yapprox <- Phi %*% H

    errorMolecule <- sqrt(mean(abs((Yapprox - as.matrix(Yi))^2)))

    list(H = H, Yapprox = Yapprox, errorMolecule = errorMolecule)
}

# CreateLinearCompound creates a network with n molecules
CreateLinearCompound <- function(n) {
    Omega <- 2 * rep(1, n)
    Omega[1] <- 3
    Omega[length(Omega)] <- 3

    B <- rep(1, n-1)

    list(Omega = Omega, B = B)
}

# InitialPosition defines the initial positions of the n molecules
InitialPosition <- function(X, n) {
    minValues <- purrr::map_dbl(X, min)
    maxValues <- purrr::map_dbl(X, max)

    ranges <- maxValues - minValues

    rangesMatrix <- matrix(ranges, nrow = n, ncol = length(ranges), byrow = TRUE)
    rUnifMatrix <- matrix(runif(n*length(ranges)), nrow = n, ncol = length(ranges), byrow = TRUE)
    minMatrix <- matrix(minValues, nrow = n, ncol = length(ranges), byrow = TRUE)
    posMolecules <- rangesMatrix * rUnifMatrix + minMatrix

    ###
    # rangesMatrix[1, ] <- minValues + (ranges / (n * 2))
    # for (k in seq_len(n-1)) {
    #     rangesMatrix[(k+1), ] <- minValues + (ranges * (2 * k + 1)/ (n * 2))
    # }
    # posMolecules <- rangesMatrix
    ###

    posMolecules
}

# DataInMolecules groups the data to the closer molecule
DataInMolecules <- function(Sigma, posMolecules) {
    X <- Sigma$X
    Y <- Sigma$Y

    numSamples <- nrow(X)
    n <- nrow(posMolecules)

    # dist <- matrix(0, nrow = numSamples, ncol = n)
    #
    # for (i in seq_len(n)) {
    #     for (j in seq_len(numSamples)) {
    #         dist[j, i] <- Distance(X[j, ], posMolecules[i, ])
    #     }
    # }

    dist <- as.matrix(pdist::pdist(as.matrix(X), posMolecules))

    moleculeNumber <- factor(paste("molecule", apply(dist, 1, which.min), sep = ""))

    SigmaSplitX <- split(X, moleculeNumber)
    SigmaSplitY <- split(Y, moleculeNumber)
    SigmaSplit <- list(X = SigmaSplitX, Y = SigmaSplitY)
    SigmaSplit

    list(SigmaSplit = SigmaSplit, moleculesUsed = sort(unique(apply(dist, 1, which.min))))
}

# RelocateMolecules relocates molecules with no data
RelocateMolecules <- function(posMolecules, moleculesUsed, errorMolecule) {
    n <- nrow(posMolecules)

    errorMolecules <- errorMolecule[moleculesUsed, , drop = FALSE]
    moleculeNumber <- sort(apply(errorMolecules, 1, sum), decreasing = TRUE, index.return = TRUE)$ix

    requireRelocation <- !(nrow(errorMolecules) == n)

    if (requireRelocation) {
        seekMolecule <- 1

        for (i in seq_len(n)) {
            if (!(i %in% moleculesUsed)) {
                radius <- 1
                posMolecules[i, ] <- rep(radius, ncol(posMolecules)) * runif(ncol(posMolecules)) + posMolecules[moleculeNumber[seekMolecule], ]

                seekMolecule <- seekMolecule + 1

                if (seekMolecule > length(moleculeNumber)) {
                    seekMolecule <- 1
                }
            }
        }
    }

    list(requireRelocation = requireRelocation, posMolecules = posMolecules)
}

# SimDataInMolecules groups the data to the closer molecule based on the trained network
SimDataInMolecules <- function(X, posMolecules) {
    numSamples <- nrow(X)
    n <- nrow(posMolecules)

    index <- 1:numSamples

    # dist <- matrix(0, nrow = numSamples, ncol = n)
    #
    # for (i in seq_len(n)) {
    #     for (j in seq_len(numSamples)) {
    #         dist[j, i] <- Distance(X[j, ], posMolecules[i, ])
    #     }
    # }

    dist <- as.matrix(pdist::pdist(as.matrix(X), posMolecules))

    moleculeNumber <- factor(paste("molecule", apply(dist, 1, which.min), sep = ""))

    SigmaSplitX <- split(X, moleculeNumber)
    SigmaSplitIndex <- split(index, moleculeNumber)
    SigmaSplit <- list(X = SigmaSplitX, Index = SigmaSplitIndex)
    SigmaSplit
}


# CreateLabels creates the labels for the network visualization
CreateLabels <- function(molecule, dimensions) {
    molecule_split <- split(molecule, ceiling(seq_len(length(molecule)) / (length(molecule) / dimensions)))
    do.call(paste, c(molecule_split, sep = "\n"))
}

# CreateNodes creates the data frame for the nodes of the network
CreateNodesEdges <- function(ahn) {
    label <- vector('numeric')
    id <- vector('character')
    group <- vector('character')
    from <- paste('C', seq_len(ahn$network$n - 1), sep = '')
    to <- paste('C', seq_len(ahn$network$n - 1) + 1, sep = '')
    for (i in seq_len(ahn$network$n)) {
        label <- c(label, CreateLabels(round(ahn$network$H[[paste('molecule', i, sep = "")]][ , 1], 3), ncol(ahn$network$Pi)))

        carbon <- paste('C', i, sep = '')
        hydrogens <- paste('H', i, seq_len(ahn$network$C$Omega[i]), sep = '')

        id <- c(id, carbon, hydrogens)

        group <- c(group, "C", paste('H', seq_len(ahn$network$C$Omega[i]), sep = ''))

        from <- c(from, rep(carbon, ahn$network$C$Omega[i]))
        to <- c(to, hydrogens)
    }

    nodes <- data.frame(id = id,
                        shape = 'circle',
                        shadow = FALSE,
                        label = label,
                        group = group,
                        font.color =c ("red", rep("black", length(id)-1)))

    edges <- data.frame(from = from,
                        to = to)

    list(nodes = nodes, edges = edges)
}

CreateTable <- function(ahn) {
    for (i in seq_len(ahn$network$n)) {
        cat("Molecule ", i,":\n", sep = "")

        molecule <- round(ahn$network$H[[paste('molecule', i, sep = "")]][ , 1], 3)
        dimensions <- ncol(ahn$network$Pi)
        molecule_split <- split(molecule, ceiling(seq_len(length(molecule)) / (length(molecule) / dimensions)))
        names(molecule_split) <- ahn$variableNames
        molecule_table <- as.data.frame(molecule_split)

        carbon <- paste('C', i, sep = '')
        hydrogens <- paste('H', i, seq_len(ahn$network$C$Omega[i]), sep = '')
        rownames(molecule_table) <- c(carbon, hydrogens)

        print(molecule_table)
        cat("\n")
    }
}
