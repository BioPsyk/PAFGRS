#' Convert Pedigree Data to Graph
#'
#' This function converts a pedigree data frame to a graph representation using the `igraph` package.
#' It was inspired by a similar function in the `FamAgg` package.
#'
#' It first checks if the input is of class "pedigree" or "pedigreeList" and converts it to a data frame.
#' It then ensures the pedigree data has the required columns and creates a graph representation using the `igraph` package.
#'
#' @param ped A data frame or an object of class "pedigree" or "pedigreeList" containing pedigree data.
#'
#' @return An `igraph` graph object representing the pedigree data.
#'
#' @examples
#' # Create a sample pedigree data frame
#' ped_data <- data.frame(
#'   family = c(1, 1, 1),
#'   id = c(1, 2, 3),
#'   father = c(NA, NA, 1),
#'   mother = c(NA, NA, 2),
#'   sex = c("M", "F", "M")
#' )
#'
#' # Convert the pedigree data to a graph
#' graph <- ped2graph(ped_data)
#'
#' # Plot the graph
#' plot(graph)
#'
#' @seealso \url{https://github.com/EuracBiomedicalResearch/FamAgg} for the original `FamAgg` package.
#'
#' @rdname pedigree_to_graph
#' @importFrom igraph graph
#' @export
ped2graph <- function(ped){

    if(is(ped, "pedigree") | is(ped, "pedigreeList")) ped <- .ped2df(ped)

    ped <- .checkPedCol(ped)

    ## define end/starting points
    ped[!ped[, "father"] %in% ped$id, "father"] <- NA
    ped[!ped[, "mother"] %in% ped$id, "mother"] <- NA
    eds <- apply(ped[, c("id", "father", "mother")], MARGIN=1, function(x){
        retval <- c()
        if(!is.na(x["father"]))
            retval <- x[c("father", "id")]
        if(!is.na(x["mother"]))
            retval <- c(retval, x[c("mother", "id")])
        return(retval)
    })

    Edges <- unlist(eds, use.names=FALSE)
    igr <- graph(edges=as.character(Edges))
    igr
}

## Helpers
.ped2df <- function(ped){
    ped <- do.call(cbind, ped)
    if(!any(colnames(ped)=="famid")) ped <- cbind(ped, famid=rep(1, nrow(ped)))
    ped <- ped[ , c("famid", "id", "findex", "mindex", "sex")]
    ## FIXME!!! Can I get NAs as findex too?
    ## have to rematch the fater id and mother id...
    notZ <- which(ped[, "findex"] > 0)
    ped[notZ, "findex"] <- ped[ped[notZ, "findex"], "id"]
    notZ <- which(ped[, "findex"] > 0)
    ped[notZ, "mindex"] <- ped[ped[notZ, "mindex"], "id"]
    .PEDCN <- c("family", "id", "father", "mother", "sex")
    colnames(ped) <- .PEDCN
    ped <- data.frame(ped)
    rownames(ped) <- ped$id
    ped
}

.checkPedCol <- function(ped){
    if(missing(ped))
        stop("ped is required!")
    if(!all(c("id", "father", "mother") %in% colnames(ped)))
        stop("id, father, mother are required columns in ped.")
    if(!is.data.frame(ped))
        ped <- as.data.frame(ped)
    ## add a family column if not present...
    if(!any(colnames(ped) == "family"))
        ped <- cbind(ped, family=rep(1, nrow(ped)))
    rownames(ped) <- as.character(ped$id)
    ped
}

