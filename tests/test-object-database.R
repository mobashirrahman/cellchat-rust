suppressPackageStartupMessages(library(CellChat))

# No export directory or original CellChat installation is needed. The actual
# complex/cofactor tables in the input object must determine the result.
local({
  reference <- get("cellchatrs_upstream_cached", asNamespace("CellChat"))()
  quiet <- function(expr) {
    invisible(capture.output(value <- suppressWarnings(suppressMessages(expr))))
    value
  }
  check <- function(db, lr, genes) {
    x <- matrix(rep(seq_along(genes) / length(genes), 12), nrow = length(genes),
                dimnames = list(genes, paste0("c", seq_len(12))))
    object <- methods::new("CellChat", data = x, data.signaling = x, DB = db,
      LR = list(LRsig = lr), idents = factor(rep(c("A", "B"), each = 6)),
      options = list(mode = "single", datatype = "RNA"))
    a <- quiet(reference$computeCommunProb(object, nboot = 2L))
    b <- quiet(computeCommunProb(object, nboot = 2L))
    stopifnot(identical(serialize(a@net, NULL), serialize(b@net, NULL)))
  }

  mouse <- reference$CellChatDB.mouse
  lr <- mouse$interaction["TGFB1_TGFBR1_TGFBR2", , drop = FALSE]
  genes <- unique(c(lr$ligand, unlist(mouse$complex[lr$receptor, ], use.names = FALSE)))
  mouse$interaction <- lr
  check(mouse, lr, genes[!is.na(genes) & nzchar(genes)])

  human <- reference$CellChatDB.human
  complex <- rownames(human$complex)[1L]
  old <- unlist(human$complex[complex, ], use.names = FALSE)
  human$complex[complex, 1L] <- "REVIEW_ALTERNATE"
  lr <- data.frame(ligand = complex, receptor = "REVIEW_RECEPTOR", agonist = "",
    antagonist = "", co_A_receptor = "", co_I_receptor = "",
    annotation = "Secreted Signaling", row.names = "REVIEW")
  check(human, lr, unique(c(old[nzchar(old)], "REVIEW_ALTERNATE", "REVIEW_RECEPTOR")))

  fish <- reference$CellChatDB.zebrafish
  lr <- fish$interaction[1L, , drop = FALSE]
  entities <- c(lr$ligand, lr$receptor)
  genes <- unique(unlist(lapply(entities, function(g) {
    if (g %in% rownames(fish$complex)) unlist(fish$complex[g, ], use.names = FALSE) else g
  })))
  check(fish, lr, genes[!is.na(genes) & nzchar(genes)])
})
