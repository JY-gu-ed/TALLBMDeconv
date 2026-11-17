#' Open a file-backed result or the standardised 125-patient release
#' @param path Result directory.
#' @return A tall_result descriptor.
#' @export
tall_read_result <- function(path) {
  path <- normalizePath(path, winslash="/", mustWork=TRUE)
  manifest_path<-file.path(path,'metadata/run_manifest.json')
  if(file.exists(manifest_path)) {
    m<-jsonlite::read_json(manifest_path,simplifyVector=TRUE)
    if(!is.null(m$status)&&!m$status%in%c('COMPLETE','PASS_STANDARDIZED_RESULT_QA'))stop('Result is incomplete or failed QA: ',m$status)
  }
  rds <- file.path(path,"metadata/result.rds")
  if(file.exists(rds)) {
    r <- readRDS(rds); r$path <- path
    r$fits <- file.path(path,"fits",basename(r$fits)); return(r)
  }
  f <- file.path(path,"metadata/run_manifest.json")
  if(!file.exists(f)) stop("Not a TALL result directory")
  structure(list(path=path,manifest=jsonlite::read_json(f,simplifyVector=TRUE)),class="tall_result")
}
.read_tsv <- function(p) utils::read.delim(p,check.names=FALSE,stringsAsFactors=FALSE)
.read_matrix_subset <- function(path, genes=NULL, samples=NULL) {
  # Stream one component and selected genes, never all source matrices.
  con <- if(grepl("\\.gz$",path)) gzfile(path,"rt") else file(path,"rt")
  on.exit(close(con)); header <- strsplit(readLines(con,n=1L),"\t",fixed=TRUE)[[1]]
  header <- gsub('^"|"$', '', header)
  cols <- if(is.null(samples)) seq.int(2L,length(header)) else match(samples,header)
  if(anyNA(cols) || any(cols==1L)) stop("Requested sample missing")
  rows <- list()
  repeat {
    lines <- readLines(con,n=2000L); if(!length(lines)) break
    if(!is.null(genes)) { gi <- sub("\t.*$","",lines);gi<-gsub('^"|"$','',gi);lines<-lines[gi %in% genes] }
    if(length(lines)) rows[[length(rows)+1L]] <- lines
  }
  if(!length(rows)) stop("Requested genes missing")
  parsed <- strsplit(unlist(rows,use.names=FALSE),"\t",fixed=TRUE)
  rn <- vapply(parsed,`[`,character(1),1L);rn<-gsub('^"|"$','',rn)
  x <- do.call(rbind,lapply(parsed,function(z) suppressWarnings(as.numeric(z[cols]))))
  dimnames(x) <- list(rn,header[cols])
  if(!is.null(genes)) { if(any(!genes %in% rn)) stop("Requested genes missing");x<-x[genes,,drop=FALSE] }
  x
}

#' Extract composition or one file-backed expression/source matrix
#' @param result A tall_result or result directory.
#' @param what composition_rna, composition_state_rna, composition_cell_proxy,
#'   composition_state_cell_proxy, malignant_L1_composition,
#'   malignant_L1_within_malignant, expression, contribution_cpm or rna_source_fraction.
#' @param group Required for expression/source matrices.
#' @param genes,samples Optional subsets in the requested order.
#' @return A numeric matrix; composition rows are patients, expression rows genes.
#' @export
tall_extract <- function(result, what="composition_rna", group=NULL, genes=NULL, samples=NULL) {
  r <- if(is.character(result)) tall_read_result(result) else result
  if(startsWith(what,"composition") || startsWith(what,"malignant_L1")) {
    allowed<-c("composition_rna","composition_state_rna","composition_cell_proxy","composition_state_cell_proxy","malignant_L1_composition","malignant_L1_within_malignant")
    if(!what %in% allowed) stop("Unknown composition output")
    x<-.read_tsv(file.path(r$path,"composition",paste0(what,".tsv")));ids<-x[[1]];x<-as.matrix(x[,-1,drop=FALSE]);rownames(x)<-ids
    if(!is.null(samples)) { if(!all(samples%in%ids))stop("Requested sample missing"); x<-x[samples,,drop=FALSE] }
    return(x)
  }
  if(!what%in%c("expression","contribution_cpm","rna_source_fraction") || is.null(group))stop("Specify a matrix kind and group")
  folder<-if(what=="expression")"expression"else"gene_source"
  idx<-.read_tsv(file.path(r$path,folder,"index.tsv"));if(!"kind"%in%names(idx))idx$kind<-"expression"
  idx<-idx[idx$group==group&idx$kind==what,,drop=FALSE]
  if(nrow(idx)!=1L)stop("Group/output not uniquely available")
  .read_matrix_subset(file.path(r$path,idx$path),genes,samples)
}

#' Audit composition, missingness and same-fit RNA source conservation
#' @param result Result object or directory.
#' @param full Check every coarse source matrix; FALSE checks composition only.
#' @return Machine-readable check table. Failed checks cause an error.
#' @export
tall_qc <- function(result, full=TRUE) {
  r<-if(is.character(result))tall_read_result(result)else result
  fine<-tall_extract(r,"composition_state_rna");coarse<-tall_extract(r)
  gd<-.read_tsv(file.path(r$path,"metadata/group_dictionary.tsv"));mal<-gd$component_id[gd$compartment=="malignant"]
  checks<-data.frame(check=c("unique_sample_ids","nonnegative_RNA","RNA_closure"),
    pass=c(!anyDuplicated(rownames(fine)),all(is.finite(fine)&fine>=0),max(abs(rowSums(fine)-1))<1e-7))
  aggregate_ok<-all(vapply(colnames(coarse),function(g){kk<-gd$component_id[gd$coarse_group==g];length(kk)>0&&all(kk%in%colnames(fine))&&max(abs(rowSums(fine[,kk,drop=FALSE])-coarse[,g]))<1e-7},logical(1)))
  checks<-rbind(checks,data.frame(check=c("coarse_nonnegative_and_closed","coarse_equals_aggregated_states"),pass=c(all(is.finite(coarse)&coarse>=0)&&max(abs(rowSums(coarse)-1))<1e-7,aggregate_ok)))
  if(length(mal)&&"Malignant_total"%in%colnames(coarse)) checks<-rbind(checks,data.frame(check="L1_sum_to_malignant",pass=max(abs(rowSums(fine[,mal,drop=FALSE])-coarse[,"Malignant_total"]))<1e-7))
  if(full) {
    reconstructed<-NULL;missing_ok<-TRUE;source_expression_ok<-TRUE;malignant_source<-NULL
    for(g in colnames(coarse)) {
      z<-tall_extract(r,"contribution_cpm",g);e<-tall_extract(r,"expression",g)
      available<-colSums(is.finite(e))==nrow(e)
      missing_ok<-missing_ok && identical(unname(available),unname(coarse[,g]>0)) && all(colSums(is.na(e))%in%c(0,nrow(e)))
      expected<-sweep(e,2,coarse[,g],"*");expected[,coarse[,g]==0]<-0
      source_expression_ok<-source_expression_ok&&max(abs(z-expected))<1e-5
      if(g=="Malignant_total")malignant_source<-z
      reconstructed<-if(is.null(reconstructed)) z else reconstructed+z
    }
    x<-r$bulk
    if(is.null(x)) {x<-.read_tsv(file.path(r$path,"metadata/bulk_model_counts.tsv.gz"));gn<-x[[1]];x<-as.matrix(x[,-1,drop=FALSE]);rownames(x)<-gn}
    target<-sweep(x,2,colSums(x),"/")*1e6
    target<-target[rownames(reconstructed),colnames(reconstructed),drop=FALSE]
    checks<-rbind(checks,data.frame(check=c("undefined_expression_is_NA","expression_times_composition_equals_gene_contribution","gene_sources_reconstruct_modeled_bulk"),pass=c(missing_ok,source_expression_ok,max(abs(reconstructed-target))<1e-5)))
    if(length(mal)&&!is.null(malignant_source)) {
      l1source<-matrix(0,nrow(malignant_source),ncol(malignant_source))
      for(g in mal)l1source<-l1source+tall_extract(r,"contribution_cpm",g)
      checks<-rbind(checks,data.frame(check="L1_gene_sources_sum_to_malignant_source",pass=max(abs(l1source-malignant_source))<1e-5))
    }
  }
  if(!all(checks$pass))stop("Result QA failed: ",paste(checks$check[!checks$pass],collapse=", "))
  checks
}

#' Verify a frozen SHA256 manifest without editing the result
#' @param path Directory containing SHA256SUMS.txt.
#' @return File-level integrity table; failures cause an error.
#' @export
tall_verify <- function(path) {
  lines<-readLines(file.path(path,"SHA256SUMS.txt"));expected<-substring(lines,1,64);rel<-substring(lines,67)
  paths<-file.path(path,rel)
  if(any(!file.exists(paths)))stop("Manifest files missing")
  actual<-vapply(paths,digest::digest,character(1),algo="sha256",file=TRUE)
  out<-data.frame(file=rel,pass=unname(actual)==expected)
  if(!all(out$pass))stop("SHA256 mismatch: ",paste(out$file[!out$pass],collapse=", "))
  out
}
