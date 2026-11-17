# Internal helpers deliberately avoid clinical outcomes and label learning.
.rs <- function(x) if(inherits(x,"Matrix"))Matrix::rowSums(x)else rowSums(x)
.cs <- function(x) if(inherits(x,"Matrix"))Matrix::colSums(x)else colSums(x)
.check_counts <- function(x, what = "counts", allow_zero_libraries = FALSE) {
  if (!is.matrix(x) && !inherits(x, "Matrix")) stop(what, " must be a genes x samples/cells matrix")
  if (is.null(rownames(x)) || is.null(colnames(x)) || anyNA(rownames(x)) ||
      anyNA(colnames(x)) || anyDuplicated(rownames(x)) || anyDuplicated(colnames(x)) ||
      any(!nzchar(rownames(x))) || any(!nzchar(colnames(x)))) stop(what, ": unique nonempty row and column IDs required")
  z <- if (inherits(x, "sparseMatrix")) x@x else as.vector(x)
  if (anyNA(z) || any(!is.finite(z)) || any(z < 0) || any(abs(z - round(z)) > 1e-7) || any(z > .Machine$integer.max))
    stop(what, ": finite nonnegative integer raw counts required (not TPM/log expression)")
  libraries <- if(inherits(x,"Matrix"))Matrix::colSums(x)else colSums(x)
  if (!allow_zero_libraries && any(libraries <= 0)) stop(what, ": zero library")
  invisible(TRUE)
}
.hash <- function(x) digest::digest(x, algo = "sha256", serialize = TRUE)
.canonical <- function(x, aliases = NULL) {
  x <- sub("^BCH:", "", as.character(x))
  if (!is.null(aliases)) {
    m <- match(x, aliases$sample_id); ok <- !is.na(m); x[ok] <- aliases$bulk_id[m[ok]]
  }
  x
}
.upstream_notes <- new.env(parent=emptyenv())
.upstream_notes$messages <- character()
.bayes_namespace <- function() {
  withCallingHandlers(loadNamespace("BayesPrism"),warning=function(w){
    msg<-conditionMessage(w)
    known<-grepl("replacing previous import",msg,fixed=TRUE)&&grepl("BayesPrism",msg,fixed=TRUE)&&(grepl("lowess",msg,fixed=TRUE)||grepl("register",msg,fixed=TRUE))
    if(known){.upstream_notes$messages<-unique(c(.upstream_notes$messages,msg));invokeRestart("muffleWarning")}
  })
}
.normalize_phi <- function(counts) {
  # Upstream BayesPrism normalization, not a reimplementation of preprocessing.
  fun <- get("norm.to.one", envir=.bayes_namespace(),inherits=FALSE)
  fun(counts, pseudo.min = 1e-8)
}
.reference_from_aggregates <- function(cube, nc, genes, groups, patients,
                                      reference_id, aliases = NULL, gene_map = NULL) {
  stopifnot(length(dim(cube)) == 3L, dim(cube)[1] == length(patients),
            dim(cube)[2] == nrow(groups), dim(cube)[3] == length(genes), all(nc >= 0))
  represented <- colSums(nc) > 0
  if (!all(represented)) stop("Patient split removes a reference group; explicitly redefine eligible groups before fitting")
  .check_hierarchy(groups)
  pooled <- matrix(colSums(matrix(cube,nrow=dim(cube)[1])),nrow=nrow(groups),ncol=length(genes))
  dimnames(pooled) <- list(groups$component_id, genes)
  base <- pooled * 0; sizes <- numeric(nrow(groups))
  for (k in seq_len(nrow(groups))) {
    ok <- which(nc[, k] > 0); a <- matrix(cube[ok, k, ], nrow = length(ok), ncol = length(genes))
    libs <- rowSums(a); stopifnot(all(libs > 0))
    base[k, ] <- colMeans(a / libs); sizes[k] <- mean(libs / nc[ok, k])
  }
  cb <- list()
  for (g in unique(groups$coarse_group)) {
    kk <- which(groups$coarse_group == g); a <- matrix(0,length(patients),length(genes))
    for(k in kk) a <- a + matrix(cube[,k,],nrow=length(patients),ncol=length(genes))
    ok <- rowSums(nc[, kk, drop=FALSE]) > 0; a <- a[ok, , drop=FALSE]
    cb[[g]] <- colMeans(a / rowSums(a))
  }
  names(sizes) <- groups$component_id
  groups$training_cells <- colSums(nc)
  groups$training_patients <- colSums(nc>0)
  groups$RNA_per_cell_proxy <- unname(sizes)
  out <- list(schema_version = "1.0", reference_id = reference_id,
    genes = genes, phi = .normalize_phi(pooled), pooled_counts = pooled,
    groups = groups, patients = as.character(patients), cell_sizes = sizes,
    baseline = base, coarse_baseline = cb, donor_counts = cube,
    donor_cell_counts = nc, aliases = aliases, gene_map = gene_map,
    gene_id_type = "symbol", validation_scope = "custom reference; no inherited atlas performance",
    data_license = "User-supplied; user must establish redistribution rights")
  out$reference_hash <- .hash(list(out$phi, out$groups, out$patients, out$cell_sizes))
  class(out) <- "tall_reference"; out
}
.check_hierarchy <- function(g) {
  if(!all(c('component_id','coarse_group','compartment','role')%in%names(g))||anyNA(g[,c('component_id','coarse_group','compartment','role')])||any(!nzchar(as.matrix(g[,c('component_id','coarse_group','compartment','role')]))))stop('Missing reference group identity or hierarchy')
  if(anyDuplicated(g$component_id))stop('Duplicate reference component IDs')
  for(p in intersect(g$component_id,g$coarse_group))if(!identical(as.character(g$component_id[g$coarse_group==p]),p))stop('Ambiguous group namespace: a parent ID collides with a different fine component')
  invisible(TRUE)
}

#' Build a donor-aware reference from raw single-cell counts
#' @param scrna Genes by cells raw integer counts, dense or sparse.
#' @param cell_annotations Data frame with cell_id and group. Optional coarse_group,
#'   compartment, role and L1/L2/L3 columns are retained; no hierarchy is required.
#' @param patient_ids Named vector mapping cell IDs to patient IDs, or a vector in
#'   the exact input column order.
#' @param reference_id Unique version identifier.
#' @param save_path Optional RDS path. Existing files are never overwritten.
#' @return A donor-aware tall_reference object.
#' @export
tall_build_reference <- function(scrna, cell_annotations, patient_ids,
                                 reference_id = "custom_v1", save_path = NULL) {
  .check_counts(scrna, "scrna"); a <- as.data.frame(cell_annotations)
  if (!all(c("cell_id", "group") %in% names(a)) || anyDuplicated(a$cell_id)) stop("Unique cell_id and group required")
  a$cell_id<-as.character(a$cell_id);a$group<-as.character(a$group)
  a <- a[match(colnames(scrna), a$cell_id), , drop=FALSE]
  if (anyNA(a$cell_id) || anyNA(a$group) || any(!nzchar(a$group))) stop("Unmapped cell/group")
  if (!is.null(names(patient_ids))) { if(anyDuplicated(names(patient_ids)))stop("Duplicate donor cell mapping");patient_ids <- patient_ids[colnames(scrna)] }
  if (length(patient_ids) != ncol(scrna) || anyNA(patient_ids) || any(!nzchar(patient_ids))) stop("Patient IDs must cover every cell")
  st <- unique(as.character(a$group)); pt <- sort(unique(as.character(patient_ids)))
  if(length(st)<2L)stop("Reference needs at least two represented groups")
  unsafe <- function(z) any(grepl("/",z,fixed=TRUE)|grepl("\\",z,fixed=TRUE)|grepl(":",z,fixed=TRUE)|grepl("[[:cntrl:]]",z))
  if(unsafe(st))stop("Group IDs must be safe filename identifiers; keep descriptive labels in separate columns")
  if (!"coarse_group" %in% names(a)) a$coarse_group <- a$group
  if (!"compartment" %in% names(a)) a$compartment <- "unspecified"
  if (!"role" %in% names(a)) a$role <- "primary_report"
  if(unsafe(a$coarse_group))stop("Coarse group IDs must be safe filename identifiers")
  extra <- intersect(c("coarse_group", "compartment", "role", "L1", "L2", "L3"), names(a))
  for(f in extra)a[[f]]<-as.character(a[[f]])
  # A component must have one consistent hierarchy; no labels are guessed.
  for (k in st) for (f in extra) if (length(unique(a[a$group == k, f])) != 1L) stop("Inconsistent hierarchy for ", k)
  groups <- a[match(st, a$group), c("group", extra), drop=FALSE]; names(groups)[1] <- "component_id"
  cube <- array(0, c(length(pt), length(st), nrow(scrna))); nc <- matrix(0, length(pt), length(st))
  for (d in seq_along(pt)) for (k in seq_along(st)) {
    ii <- which(patient_ids == pt[d] & a$group == st[k]); nc[d, k] <- length(ii)
    if (length(ii)) cube[d, k, ] <- .rs(scrna[, ii, drop=FALSE])
  }
  out <- .reference_from_aggregates(cube, nc, rownames(scrna), groups, pt, reference_id)
  if (!is.null(save_path)) { if (file.exists(save_path)) stop("Reference exists; use a new version path"); saveRDS(out, save_path) }
  out
}

#' Load a versioned reference
#' @param path Reference RDS path; NULL selects the installed default atlas.
#' @param version default, benchmark73 or fullAtlas81.
#' @return A validated tall_reference object.
#' @export
tall_load_reference <- function(path = NULL, version = c("default", "benchmark73", "fullAtlas81")) {
  version <- match.arg(version)
  if (is.null(path)) {
    if (version == "default") version <- "fullAtlas81"
    path <- system.file("references", paste0(version, ".rds"), package="TALLBMDeconv")
    if (!nzchar(path)) stop("Default atlas resource not installed. Install the authorised resource bundle, or provide a custom reference RDS")
  }
  r <- readRDS(path)
  if (!is.null(r$donor_resource) && !file.exists(r$donor_resource))
    r$donor_resource <- file.path(dirname(path), basename(r$donor_resource))
  .validate_reference(r)
  r
}
.validate_reference <- function(r) {
  if (!inherits(r, "tall_reference") || !all(c("phi", "groups", "patients", "cell_sizes", "reference_id") %in% names(r))) stop("Invalid reference schema")
  .check_hierarchy(r$groups)
  if (any(!is.finite(r$phi)) || any(r$phi < 0) || max(abs(rowSums(r$phi) - 1)) > 1e-8 ||
      !identical(rownames(r$phi), r$groups$component_id) || !identical(colnames(r$phi), r$genes) ||
      any(!is.finite(r$cell_sizes)) || any(r$cell_sizes <= 0)) stop("Invalid reference parameters or ordering")
  expected <- .hash(list(r$phi, r$groups, r$patients, r$cell_sizes))
  if (!identical(expected, r$reference_hash)) stop("Reference hash mismatch")
  invisible(TRUE)
}
