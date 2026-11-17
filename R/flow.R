#' Refit a result with typed flow observations
#' @param result A tall_result with retained counts and reference.
#' @param flow Audited observation table accepted by tall_deconv.
#' @param out_dir New output directory; existing results are never overwritten.
#' @return A new tall_result from a genuine constrained fit.
#' @export
tall_add_flow <- function(result, flow, out_dir = tempfile("TALL_flow_")) {
  if (is.character(result)) result <- tall_read_result(result)
  if (is.null(result$bulk) || is.null(result$reference)) stop("Result lacks retained counts/reference needed for refitting")
  out<-tall_deconv(result$bulk, reference=result$reference, flow=flow, config=result$config, out_dir=out_dir)
  p<-file.path(out$path,"metadata/parent_run.json")
  jsonlite::write_json(list(parent_result=result$path,refitted_input="retained_modeled_raw_counts",parent_reference_id=result$reference$reference_id,
    parent_input_libraries=data.frame(sample_id=result$samples,raw_library_size=if(is.null(result$input_libraries))colSums(result$bulk)else result$input_libraries)),p,auto_unbox=TRUE,pretty=TRUE)
  cat(paste0(digest::digest(p,file=TRUE,algo="sha256"),"  metadata/parent_run.json\n"),file=file.path(out$path,"SHA256SUMS.txt"),append=TRUE)
  out
}
.prepare_flow <- function(flow, reference, samples) {
  if (is.null(flow)) return(NULL)
  f <- as.data.frame(flow, stringsAsFactors=FALSE)
  req <- c("sample_id", "group", "kind", "denominator_scope", "original_denominator", "gate_identity")
  if (!all(req %in% names(f))) stop("Flow needs: ", paste(req, collapse=", "))
  if (anyNA(f[, req]) || any(!nzchar(as.matrix(f[, req])))) stop("Flow identity/denominator fields cannot be missing")
  if (!all(f$sample_id %in% samples)) stop("Flow sample ID not present in bulk")
  if (!all(f$kind %in% c("point", "lower", "interval", "unknown"))) stop("Unsupported flow kind")
  active <- f$kind != "unknown"
  if (any(active & !f$denominator_scope %in% c("modeled_nucleated", "raw_event_upper_bound_on_modeled_nucleated_count")))
    stop("Conditional or unverified denominator unsupported; provide an explicitly audited unconditional modeled-cell observation")
  if (any(active & f$denominator_scope == "raw_event_upper_bound_on_modeled_nucleated_count" & f$kind != "lower")) stop("Raw-event denominator only supports the audited conservative lower bound")
  valid <- unique(c(reference$groups$component_id, reference$groups$coarse_group))
  if (any(active & !f$group %in% valid)) stop("Unmapped flow group")
  for (v in c("value", "lower", "upper", "measurement_error")) if (!v %in% names(f)) f[[v]] <- rep(NA_real_,nrow(f))
  if (any(active & (!is.finite(f$measurement_error) | f$measurement_error <= 0))) stop("Positive measurement_error required")
  lo <- ifelse(f$kind == "point", f$value, f$lower)
  hi <- ifelse(f$kind == "point", f$value, ifelse(f$kind == "lower", 1, f$upper))
  if (any(active & (!is.finite(lo) | !is.finite(hi) | lo < 0 | hi > 1 | lo > hi))) stop("Invalid flow fraction bounds")
  f$effective_lower <- lo; f$effective_upper <- hi; f$requested <- active;f$used <- active
  f
}
.flow_arrays <- function(flow, sample, ref) {
  if (is.null(flow)) return(NULL)
  f <- flow[flow$sample_id == sample & flow$used, , drop=FALSE]
  if (!nrow(f)) return(NULL)
  G <- matrix(0, nrow(f), nrow(ref$phi))
  for (h in seq_len(nrow(f))) {
    G[h, ] <- as.numeric(ref$groups$component_id == f$group[h] | ref$groups$coarse_group == f$group[h])
    if (!any(G[h, ] > 0)) stop("Empty flow mapping")
  }
  list(groups=G, lo=f$effective_lower, hi=f$effective_upper, scales=f$measurement_error)
}
