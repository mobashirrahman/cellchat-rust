# Inspect the pinned visium spatial object: what computeRegionDistance would actually be given.
CC <- Sys.getenv("CELLCHAT_SRC", "../CellChat")
q <- function(x) { invisible(utils::capture.output(v <- suppressWarnings(suppressMessages(x)))); v }
e <- new.env(parent = globalenv())
for (p in c("methods", "Matrix", "S4")) suppressWarnings(suppressMessages(requireNamespace(p, quietly = TRUE)))
q(sys.source(file.path(CC, "R", "CellChat_class.R"), envir = e, keep.source = FALSE))

o <- readRDS(file.path(Sys.getenv("CELLCHATRS_DATA", "data"), "visium.rds"))
say <- function(...) cat("V ", ..., "\n", sep = "")
say("datatype:", o@options$datatype)
say("data:", paste(dim(o@data), collapse = "x"))
say("data.signaling:", paste(dim(o@data.signaling), collapse = "x"))
say("nLR:", nrow(o@LR$LRsig))
say("meta cols:", paste(colnames(o@meta), collapse = ","))
say("idents levels:", paste(levels(o@idents), collapse = ","))
tb <- table(o@idents)
say("idents sizes:", paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = " "))
crd <- o@images$coordinates
say("coords:", length(crd$x_cent), "x-centroid",
    "x range", format(range(crd$x_cent), digits = 6),
    "y range", format(range(crd$y_cent), digits = 6))
say("spot.diameter:", format(o@images$scale.factors$spot.diameter, digits = 8))
say("slices:", paste(unique(o@meta$slices), collapse = ","))
say("var.features:", length(o@var.features$features))
say("ncol(coord) == ncol(data):", length(crd$x_cent) == ncol(o@data))
say("colnames(coord) == colnames(data):", identical(rownames(crd), colnames(o@data)))
say("rownames(crd)[1:3]:", paste(head(rownames(crd), 3), collapse = ","))
say("colnames(data)[1:3]:", paste(head(colnames(o@data), 3), collapse = ","))
