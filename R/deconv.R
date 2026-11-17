.config <- function(config) {
  defaults <- list(iterations=40000L, burnin=20000L, thin=4L,
                   seeds=c(123L, 20261110L), alpha=1e-8, flow_strength=1)
  if(!is.list(config)||anyDuplicated(names(config))||(length(config)&& (is.null(names(config))||any(!nzchar(names(config))))))stop("Configuration must be a uniquely named list")
  if (length(setdiff(names(config), names(defaults)))) stop("Unknown configuration field")
  x <- utils::modifyList(defaults, config)
  scalar<-setdiff(names(defaults),'seeds')
  if(any(!vapply(x,is.numeric,logical(1)))||any(lengths(x[scalar])!=1L))stop("Sampler controls must be numeric scalars except the two seeds")
  if (any(!is.finite(unlist(x))) || x$iterations <= x$burnin || x$burnin < 0 ||
      x$thin < 1 || length(x$seeds) != 2L || anyDuplicated(x$seeds) ||
      x$alpha <= 0 || x$flow_strength < 0 ||
      any(c(x$iterations, x$burnin, x$thin, x$seeds) != round(c(x$iterations, x$burnin, x$thin, x$seeds)))) stop("Invalid sampler configuration; two independent integer seeds required")
  if (length(seq.int(x$burnin + x$thin, x$iterations, x$thin)) < 2L) stop("Need at least two retained draws per chain")
  x
}
.frozen_rng <- function(seed, code) {
  kind<-RNGkind();existed<-exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE)
  if(existed)previous<-get('.Random.seed',envir=.GlobalEnv)
  on.exit({do.call(RNGkind,as.list(kind));if(existed)assign('.Random.seed',previous,envir=.GlobalEnv)else if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE))rm('.Random.seed',envir=.GlobalEnv)})
  RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(seed);force(code)
}
.fit_one <- function(x, ref, flow, cfg) {
  keep <- as.integer(seq.int(cfg$burnin + cfg$thin, cfg$iterations, cfg$thin))
  lapply(cfg$seeds, function(seed) .frozen_rng(seed, {
    if (is.null(flow)) tbm_bp_initial(as.numeric(x), ref$phi, cfg$alpha, keep)
    else tbm_flow_sample(as.numeric(x), ref$phi, cfg$alpha, keep,
                        as.numeric(ref$cell_sizes), flow$groups, flow$lo, flow$hi,
                        flow$scales, ncol(ref$phi) * cfg$flow_strength)
  }))
}
.backend_state <- new.env(parent=emptyenv())
.assert_backend <- function() {
  if(isTRUE(.backend_state$checked))return(invisible(TRUE))
  .bayes_namespace()
  existed<-exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if(existed)previous<-get(".Random.seed",envir=.GlobalEnv)
  on.exit(if(existed)assign(".Random.seed",previous,envir=.GlobalEnv)else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))rm(".Random.seed",envir=.GlobalEnv))
  phi<-matrix(c(.8,.2,.2,.8),2);keep<-seq.int(52L,300L,2L)
  upstream<-.frozen_rng(123,BayesPrism::sample.Z.theta_n(c(7,13),phi,.1,keep,FALSE))
  fast<-.frozen_rng(123,tbm_bp_initial(c(7,13),phi,.1,keep))
  if(max(abs(upstream$Z_n-fast$Z_n))>1e-10||max(abs(upstream$theta_n-fast$theta_n))>1e-10)stop("Installed BayesPrism sampler fails the frozen execution-kernel oracle")
  .backend_state$checked<-TRUE;invisible(TRUE)
}
.write_tsv <- function(x, path) {
  utils::write.table(x, file=path, sep="\t", row.names=FALSE, quote=FALSE, na="NA")
}
.write_matrix <- function(x, path) {
  con <- gzfile(path, "wt"); on.exit(close(con))
  utils::write.table(data.frame(gene=rownames(x), x, check.names=FALSE), con,
                     sep="\t", quote=FALSE, row.names=FALSE, na="NA")
}
.rhat <- function(chains) {
  n <- nrow(chains[[1]]); W <- rowMeans(sapply(chains, function(x) apply(x, 2, stats::var)))
  B <- n * apply(sapply(chains, colMeans), 1, stats::var)
  out <- sqrt(((n-1)/n*W + B/n)/W); out[!is.finite(out)] <- NA_real_; out
}

#' Run the fixed-reference initial-allocation workflow
#' @param bulk Genes by patients raw integer count matrix; unique symbols by default.
#' @param reference A tall_reference; NULL loads the installed full atlas.
#' @param flow Optional audited table. Required fields: sample_id, group, kind,
#'   denominator_scope, original_denominator, gate_identity; value for point,
#'   lower for lower, lower/upper for interval, and positive measurement_error.
#'   Fractions are 0--1. Unknown rows are retained but never used as zero.
#' @param input_type Only raw counts are supported.
#' @param gene_id_type symbol or ensembl. Ensembl conversion needs a unique
#'   gene_id/display_symbol reference dictionary. Ambiguous IDs are rejected.
#' @param config Named sampler overrides. Defaults reproduce the frozen two-chain
#'   40000/20000/thin4 initial allocation; reference updating is not performed.
#' @param out_dir New file-backed output directory.
#' @return A tall_result; expression is read on demand with tall_extract.
#' @export
tall_deconv <- function(bulk, reference=NULL, flow=NULL, input_type="counts",
                        gene_id_type="symbol", config=list(), out_dir=tempfile("TALL_run_")) {
  if (!identical(input_type, "counts")) stop("Only raw integer counts supported")
  .check_counts(bulk, "bulk"); input_libraries <- .cs(bulk)
  ref <- if (is.null(reference)) tall_load_reference() else reference
  if (is.character(ref)) ref <- tall_load_reference(ref)
  if (!inherits(ref, "tall_reference")) stop("Invalid reference")
  .validate_reference(ref)
  if (as.character(utils::packageVersion("BayesPrism")) != "2.2.3") stop("This release is audited against BayesPrism 2.2.3; install that version")
  .assert_backend()
  mapping <- data.frame(input_id=rownames(bulk), model_id=rownames(bulk), stringsAsFactors=FALSE)
  if (gene_id_type == "ensembl") {
    if (is.null(ref$gene_map) || !all(c("gene_id", "display_symbol") %in% names(ref$gene_map))) stop("Ensembl mapping dictionary unavailable")
    gid <- sub("\\.[0-9]+$", "", rownames(bulk)); gm <- ref$gene_map
    if (anyDuplicated(gid)) stop("Duplicate Ensembl IDs after version stripping")
    sym <- gm$display_symbol[match(gid, gm$gene_id)]; mapped <- !is.na(sym) & nzchar(sym)
    if (anyDuplicated(sym[mapped])) stop("Ambiguous symbol mapping; resolve explicitly")
    mapping$model_id <- sym; bulk <- bulk[mapped, , drop=FALSE]; rownames(bulk) <- sym[mapped]
  } else if (gene_id_type != "symbol") stop("Unsupported gene_id_type")
  missing <- setdiff(ref$genes, rownames(bulk))
  if (length(missing)) stop(length(missing), " model genes missing; no zero filling. Build an explicitly gene-matched reference. First: ", paste(utils::head(missing, 5), collapse=", "))
  retained <- .cs(bulk[ref$genes, , drop=FALSE]) / input_libraries
  bulk <- as.matrix(bulk[ref$genes, , drop=FALSE]); .check_counts(bulk)
  cfg <- .config(config); fl <- .prepare_flow(flow, ref, colnames(bulk))
  if(!is.null(fl)&&cfg$flow_strength==0){fl$used<-FALSE;fl$disabled_reason<-'flow_strength_zero'}
  if (dir.exists(out_dir) || file.exists(out_dir)) stop("Output exists; use a new version directory")
  dir.create(out_dir, recursive=TRUE); out_dir <- normalizePath(out_dir, winslash="/")
  for (d in c("composition", "expression/malignant", "expression/nonmalignant", "gene_source", "quality", "metadata", "fits")) dir.create(file.path(out_dir, d), recursive=TRUE, showWarnings=FALSE)
  samples <- colnames(bulk); paths <- file.path(out_dir, "fits", sprintf("patient_%04d.rds", seq_along(samples)))
  manifest <- list(schema_version="1.0", status="RUNNING", algorithm="BayesPrism2.2.3 fixed_phi_initial_read_allocation", reference_id=ref$reference_id, reference_hash=ref$reference_hash, reference_update=FALSE, config=cfg,rng_kind=c('Mersenne-Twister','Inversion','Rejection'), samples=samples, input_hash=.hash(bulk),upstream_namespace_notes=.upstream_notes$messages)
  jsonlite::write_json(manifest, file.path(out_dir,"metadata/run_manifest.json"), auto_unbox=TRUE, pretty=TRUE)
  for (i in seq_along(samples)) {
    message("Fitting ", i, "/", length(samples), ": ", samples[i])
    fits <- .fit_one(bulk[, i], ref, .flow_arrays(fl, samples[i], ref), cfg)
    z <- (fits[[1]]$Z_n + fits[[2]]$Z_n)/2
    if (any(!is.finite(z)) || any(z < 0) || max(abs(rowSums(z)-bulk[, i])) > 1e-6) stop("Allocation conservation failed for ", samples[i])
    dimnames(z) <- list(ref$genes, ref$groups$component_id)
    saveRDS(list(Z=z, chains=fits), paths[i], compress=FALSE)
  }
  result <- structure(list(path=out_dir, bulk=bulk, reference=ref, config=cfg, flow=fl,
                           samples=samples, fits=paths, retained_read_fraction=retained,
                           input_libraries=input_libraries, gene_mapping=mapping), class="tall_result")
  .export_result(result)
  saveRDS(result, file.path(out_dir,"metadata/result.rds"))
  checks<-tryCatch(tall_qc(result),error=function(e){manifest$status<-'FAILED_QA';manifest$error<-conditionMessage(e);jsonlite::write_json(manifest,file.path(out_dir,'metadata/run_manifest.json'),auto_unbox=TRUE,pretty=TRUE);stop(e)})
  .write_tsv(checks, file.path(out_dir,"quality/run_qc.tsv"))
  manifest$status <- "COMPLETE"; manifest$flow_observations <- if(is.null(fl)) 0L else sum(fl$used)
  jsonlite::write_json(manifest, file.path(out_dir,"metadata/run_manifest.json"), auto_unbox=TRUE, pretty=TRUE)
  files <- list.files(out_dir, recursive=TRUE, full.names=TRUE)
  hashes <- vapply(files, digest::digest, character(1), file=TRUE, algo="sha256")
  writeLines(paste(hashes, substring(files, nchar(out_dir)+2L), sep="  "), file.path(out_dir,"SHA256SUMS.txt"))
  result
}
.export_result <- function(r) {
  out <- r$path; ref <- r$reference; ids <- r$samples; st <- ref$groups$component_id
  coarse <- unique(ref$groups$coarse_group); groups <- unique(c(coarse, st))
  f <- matrix(0, length(ids), length(st), dimnames=list(ids, st)); quality <- list()
  for (i in seq_along(ids)) {
    fit <- readRDS(r$fits[i]); f[i, ] <- colSums(fit$Z)/sum(r$bulk[, i]); traces <- lapply(fit$chains, `[[`, "theta_trace")
    extra <- setdiff(coarse,st); qnames<-c(st,extra); qfraction<-c(f[i,],vapply(extra,function(g)sum(f[i,ref$groups$coarse_group==g]),numeric(1)))
    if(length(extra)) traces<-lapply(traces,function(t)cbind(t,sapply(extra,function(g)rowSums(t[,ref$groups$coarse_group==g,drop=FALSE]))))
    draws <- do.call(rbind, traces); rh <- .rhat(traces)
    quality[[i]] <- data.frame(sample_id=ids[i], component=qnames, RNA_fraction=qfraction, theta_mean=colMeans(draws),
      theta_lower95=apply(draws,2,stats::quantile,.025), theta_upper95=apply(draws,2,stats::quantile,.975),
      Rhat=rh, Rhat_review=!is.na(rh)&rh>=1.1, zero_mass=qfraction==0,
      low_contribution_lt_1pct=qfraction<.01, low_contribution_lt_01pct=qfraction<.001,
      absolute_chain_mean_difference_pp=100*abs(colMeans(traces[[1]])-colMeans(traces[[2]])),
      expression_uncertainty="chain_disagreement_not_credible_interval")
  }
  cfrac <- matrix(sapply(coarse, function(g) rowSums(f[, ref$groups$coarse_group==g, drop=FALSE])),nrow=length(ids),dimnames=list(ids,coarse))
  proxy <- sweep(f, 2, ref$cell_sizes, "/"); proxy <- proxy / rowSums(proxy)
  cp <- matrix(sapply(coarse, function(g) rowSums(proxy[, ref$groups$coarse_group==g, drop=FALSE])),nrow=length(ids),dimnames=list(ids,coarse))
  wt <- function(x,p) .write_tsv(data.frame(patient_id=rownames(x), x, check.names=FALSE),file.path(out,"composition",p))
  wt(f,"composition_state_rna.tsv"); wt(cfrac,"composition_rna.tsv"); wt(proxy,"composition_state_cell_proxy.tsv"); wt(cp,"composition_cell_proxy.tsv")
  mal <- ref$groups$compartment=="malignant"
  if (any(mal)) {
    wt(f[,mal,drop=FALSE],"malignant_L1_composition.tsv")
    den <- rowSums(f[,mal,drop=FALSE]); den[den==0] <- NA_real_
    wt(f[,mal,drop=FALSE]/den,"malignant_L1_within_malignant.tsv")
  }
  index <- list(); missing <- list(); reliability <- do.call(rbind,quality)
  for (g in groups) {
    kk <- which(st==g | ref$groups$coarse_group==g)
    cm <- matrix(NA_real_, length(ref$genes), length(ids), dimnames=list(ref$genes, ids)); contribution <- cm; delta <- cm
    for (i in seq_along(ids)) {
      fit <- readRDS(r$fits[i]); z <- rowSums(fit$Z[,kk,drop=FALSE]); mass <- sum(z)
      cm[,i] <- if(mass>0) z/mass*1e6 else NA_real_; contribution[,i] <- z/sum(r$bulk[,i])*1e6
      cz <- lapply(fit$chains,function(x) rowSums(x$Z_n[,kk,drop=FALSE])); den <- vapply(cz,sum,numeric(1))
      if (all(den>0)) delta[,i] <- abs(log2(1+cz[[1]]/den[1]*1e6)-log2(1+cz[[2]]/den[2]*1e6))
    }
    compartment <- if (all(ref$groups$compartment[kk]=="malignant")) "malignant" else "nonmalignant"
    ep <- paste0("expression/",compartment,"/",g,"_cpm.tsv.gz")
    .write_matrix(cm,file.path(out,ep)); .write_matrix(delta,file.path(out,"quality",paste0(g,"_chain_abs_logCPM_difference.tsv.gz")))
    index[[length(index)+1L]] <- data.frame(group=g,kind="expression",path=ep,scale="within_group_CPM",storage="tsv.gz")
    bulkcpm <- sweep(r$bulk,2,colSums(r$bulk),"/")*1e6; frac <- contribution/bulkcpm; frac[bulkcpm==0] <- NA_real_
    for (kind in c("contribution_cpm","rna_source_fraction")) {
      p <- paste0("gene_source/",g,"_",kind,".tsv.gz"); .write_matrix(if(kind=="contribution_cpm") contribution else frac,file.path(out,p))
      index[[length(index)+1L]] <- data.frame(group=g,kind=kind,path=p,scale=kind,storage="tsv.gz")
    }
    missing[[g]] <- data.frame(group=g,available=sum(colSums(is.finite(cm))==nrow(cm)),missing=sum(colSums(is.finite(cm))!=nrow(cm)),reason_if_missing="zero_posterior_allocated_mass")
  }
  idx <- do.call(rbind,index); .write_tsv(idx[idx$kind=="expression",],file.path(out,"expression/index.tsv")); .write_tsv(idx[idx$kind!="expression",],file.path(out,"gene_source/index.tsv"))
  .write_tsv(do.call(rbind,missing),file.path(out,"quality/missingness.tsv")); .write_tsv(reliability,file.path(out,"quality/reliability.tsv"))
  used <- vapply(ids,function(id) !is.null(.flow_arrays(r$flow,id,ref)),logical(1))
  .write_tsv(data.frame(sample_id=ids,flow_used=used,input_library_size=r$input_libraries,model_library_size=colSums(r$bulk),retained_read_fraction=r$retained_read_fraction),file.path(out,"metadata/sample_manifest.tsv"))
  .write_tsv(r$gene_mapping,file.path(out,"metadata/gene_mapping.tsv"))
  .write_tsv(ref$groups,file.path(out,"metadata/group_dictionary.tsv")); .write_tsv(data.frame(patient_id=ref$patients),file.path(out,"metadata/reference_members.tsv"))
  if(!is.null(r$flow)) .write_tsv(r$flow,file.path(out,"metadata/flow_observation_audit.tsv"))
  utils::capture.output(utils::sessionInfo(),file=file.path(out,"metadata/R_SESSION.txt"))
  writeLines(c("TALLBMDeconv: fixed-reference initial allocation; no reference updating.","RNA contribution is not cell percentage. Cell proxy uses gene-matched RNA-per-cell before aggregation.","Expression is within-group CPM on model genes. Zero assigned mass has NA expression.","Gene contributions reconstruct modeled bulk CPM. Fractions are NA for unexpressed genes.","Posterior theta intervals and between-chain expression differences are distinct; the latter are not expression credible intervals.","Low-contribution flags are descriptive, not validated accuracy thresholds."),file.path(out,"README.md"))
}
