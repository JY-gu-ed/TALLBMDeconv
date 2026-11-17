.donors <- function(ref) {
  if(!is.null(ref$donor_counts))return(list(counts=ref$donor_counts,cell_counts=ref$donor_cell_counts))
  p<-ref$donor_resource
  if(is.null(p)||!file.exists(p))stop("Patient-split reference requires donor aggregates; install the authorised donor resource or rebuild with tall_build_reference")
  if(!is.null(ref$donor_resource_hash)&&!identical(digest::digest(p,file=TRUE,algo="sha256"),ref$donor_resource_hash))stop("Donor resource hash mismatch")
  readRDS(p)
}
.subset_reference <- function(ref, exclude) {
  keep<-!.canonical(ref$patients,ref$aliases)%in%.canonical(exclude,ref$aliases)
  if(!any(keep))stop("No reference patients remain after exclusion")
  d<-.donors(ref)
  out<-.reference_from_aggregates(d$counts[keep,,,drop=FALSE],d$cell_counts[keep,,drop=FALSE],
    ref$genes,ref$groups,ref$patients[keep],paste0(ref$reference_id,"_patient_excluded_",substr(.hash(sort(exclude)),1,8)),ref$aliases,ref$gene_map)
  out$validation_scope<-"patient-excluded reference built for this benchmark";out
}
.truth <- function(scrna, ref) {
  if(!is.list(scrna))stop("scrna must be list(counts, cell_annotations, patient_ids) or an audited pseudobulk truth list")
  if(!is.null(scrna$pseudobulk)) {
    if(is.null(scrna$cell_counts)||is.null(scrna$source_present))stop("Pseudobulk truth requires cell_counts and source_present")
    if(is.null(names(scrna$pseudobulk))||is.null(names(scrna$source_present))||anyDuplicated(names(scrna$pseudobulk)))stop("Pseudobulk truth requires unique named patient records")
    if(anyDuplicated(names(scrna$source_present))||anyNA(scrna$source_present)||!is.logical(scrna$source_present))stop('source_present must be a uniquely named TRUE/FALSE vector')
    if(anyDuplicated(rownames(scrna$cell_counts))||anyDuplicated(colnames(scrna$cell_counts)))stop('Duplicate truth cell-count identities')
    if(anyNA(scrna$cell_counts)||any(scrna$cell_counts<0)||any(scrna$cell_counts!=round(scrna$cell_counts)))stop("Invalid truth cell counts")
    for(p in names(scrna$pseudobulk)) {
      .check_counts(scrna$pseudobulk[[p]],"truth pseudobulk",allow_zero_libraries=TRUE)
      if(!p%in%rownames(scrna$cell_counts)||!p%in%names(scrna$source_present)||any(!colnames(scrna$pseudobulk[[p]])%in%colnames(scrna$cell_counts)))stop("Truth patient/group mapping incomplete")
      if(isTRUE(scrna$source_present[p])&&(!all(ref$genes%in%rownames(scrna$pseudobulk[[p]]))||!all(ref$groups$component_id%in%colnames(scrna$pseudobulk[[p]]))))stop('Truth missing model genes or group columns; confirmed absence requires explicit zero counts')
    }
    return(scrna)
  }
  if(!all(c("counts","cell_annotations","patient_ids")%in%names(scrna)))stop("Raw truth needs counts, cell_annotations, patient_ids")
  x<-scrna$counts;.check_counts(x,"single-cell truth");a<-scrna$cell_annotations
  if(!is.data.frame(a)||!all(c('cell_id','group')%in%names(a))||anyDuplicated(a$cell_id))stop('Truth needs unique cell_id and group annotations')
  a$cell_id<-as.character(a$cell_id);a$group<-as.character(a$group)
  a<-a[match(colnames(x),a$cell_id),,drop=FALSE]
  if(anyNA(a$cell_id)||anyNA(a$group)||any(!nzchar(a$group)))stop("Truth cell identity missing")
  ids<-scrna$patient_ids;if(!is.null(names(ids))){if(anyDuplicated(names(ids)))stop('Duplicate truth donor cell mapping');ids<-ids[colnames(x)]}
  if(length(ids)!=ncol(x)||anyNA(ids)||any(!nzchar(ids)))stop("Truth patient identity missing")
  if(!all(ref$genes%in%rownames(x)))stop("Truth missing evaluation genes")
  groups<-unique(c(ref$groups$component_id,as.character(a$group)));pts<-unique(as.character(ids));pb<-list();nc<-matrix(0,length(pts),length(groups),dimnames=list(pts,groups))
  for(p in pts) {
    z<-matrix(0,length(ref$genes),length(groups),dimnames=list(ref$genes,groups))
    for(g in groups){ii<-which(ids==p&a$group==g);nc[p,g]<-length(ii);if(length(ii))z[,g]<-.rs(x[ref$genes,ii,drop=FALSE])}
    pb[[p]]<-z
  }
  list(pseudobulk=pb,cell_counts=nc,source_present=stats::setNames(rep(TRUE,length(pts)),pts),
       patient_ids=pts)
}
.cor_safe <- function(x,y,method="pearson") {
  good<-is.finite(x)&is.finite(y)
  if(sum(good)<3||stats::sd(x[good])<1e-12||stats::sd(y[good])<1e-12)return(NA_real_)
  stats::cor(x[good],y[good],method=method)
}
.fit_flow <- function(r) {
  f<-r$flow
  if(is.null(f)){
    p<-file.path(r$path,'metadata/flow_observation_audit.tsv')
    if(file.exists(p))f<-.read_tsv(p)
  }
  # Version0.1.0 could label zero-strength observations as used. Correct the
  # interpretation in memory, retaining the immutable historical files.
  if(!is.null(f)&&!is.null(r$config$flow_strength)&&r$config$flow_strength==0)f$used<-FALSE
  f
}
.flow_signature <- function(f,p) {
  numeric_cols<-c('effective_lower','effective_upper','measurement_error')
  cols<-c('group','kind',numeric_cols,'denominator_scope','original_denominator','gate_identity')
  out<-as.data.frame(stats::setNames(rep(list(character()),length(cols)),cols),stringsAsFactors=FALSE)
  for(n in numeric_cols)out[[n]]<-numeric()
  if(is.null(f)||!nrow(f))return(out)
  if(!all(c('used','sample_id')%in%names(f)))stop('Imported flow usage cannot be verified')
  active<-f$sample_id==p&as.character(f$used)%in%c('TRUE','True','1')
  if(!any(active))return(out)
  if(!all(cols%in%names(f)))stop('Imported flow observation semantics are incomplete')
  f<-f[active,cols,drop=FALSE]
  for(n in setdiff(cols,numeric_cols))f[[n]]<-as.character(f[[n]])
  for(n in numeric_cols)f[[n]]<-as.numeric(f[[n]])
  f<-f[do.call(order,f),,drop=FALSE];rownames(f)<-NULL;f
}
.validate_fit_counts <- function(r,bulk,p,ref) {
  x<-r$bulk
  if(is.null(x)){
    a<-file.path(r$path,'metadata/bulk_model_counts.tsv.gz')
    b<-file.path(r$path,'metadata/frozen_model_input.rds')
    if(file.exists(a)){d<-.read_tsv(a);gn<-d[[1]];x<-as.matrix(d[,-1,drop=FALSE]);rownames(x)<-gn}
    else if(file.exists(b))x<-t(readRDS(b)$mixture)
  }
  if(is.null(x)||!p%in%colnames(x)||!all(ref$genes%in%rownames(x))||!all(ref$genes%in%rownames(bulk)))stop('Imported fit input counts cannot be verified')
  if(!isTRUE(all.equal(unname(as.numeric(x[ref$genes,p])),unname(as.numeric(bulk[ref$genes,p])),tolerance=0)))stop('Imported fit bulk counts differ from requested benchmark input')
  invisible(TRUE)
}
.validate_fit_scale <- function(r,ref) {
  size<-r$reference$cell_sizes;groups<-r$reference$groups
  if(is.null(size)) {
    p<-file.path(r$path,'metadata/gene_matched_cell_sizes.tsv')
    if(file.exists(p)){d<-.read_tsv(p);size<-stats::setNames(d$RNA_per_cell_proxy,d$component_id)}
  }
  if(is.null(groups)){
    p<-file.path(r$path,'metadata/group_dictionary.tsv');if(file.exists(p))groups<-.read_tsv(p)
  }
  ids<-ref$groups$component_id;cols<-c('component_id','coarse_group','compartment')
  if(is.null(size)||!all(ids%in%names(size))||!isTRUE(all.equal(unname(size[ids]),unname(ref$cell_sizes[ids]),tolerance=1e-12)))stop('Imported fit reference cell-size scale cannot be verified or differs')
  if(is.null(groups)||!all(cols%in%names(groups))||!all(ids%in%groups$component_id))stop('Imported fit hierarchy cannot be verified')
  actual<-groups[match(ids,groups$component_id),cols,drop=FALSE];expected<-ref$groups[,cols,drop=FALSE]
  if(!isTRUE(all.equal(actual,expected,check.attributes=FALSE)))stop('Imported fit reference hierarchy differs')
  invisible(TRUE)
}

#' Patient-level paired validation without reference leakage
#' @param bulk Paired genes by patients raw counts.
#' @param scrna List with counts, cell_annotations and patient_ids; alternatively
#'   named pseudobulk matrices, a patient-by-group cell_counts matrix, and named
#'   source_present. Pseudobulk matrices are genes by reference state.
#' @param reference A versioned reference.
#' @param split independent, holdout (exclude all validation patients together),
#'   or lopo (exclude the evaluated patient). Splitting requires donor aggregates.
#' @param results Optional existing tall_result/directory, or named per-patient
#'   list for LOPO. Reuses fits, never silently reruns them.
#' @param flow Optional typed observations, kept distinct from independent tests.
#' @param min_cells Minimum true cells for expression evaluation; composition
#'   retains valid zeros. Missing predictions remain NA and coverage is reported.
#' @param config Sampler overrides for new fits only.
#' @param groups Optional predeclared evaluation groups. NULL evaluates all.
#' @param out_dir New benchmark directory.
#' @return Tables for composition, expression, reference-mean comparison,
#'   patient-specific variation, coverage and leakage audit.
#' @export
tall_benchmark <- function(bulk,scrna,reference,split=c("independent","holdout","lopo"),
                           results=NULL,flow=NULL,min_cells=20L,config=list(),groups=NULL,out_dir=tempfile("TALL_benchmark_")) {
  split<-match.arg(split);.check_counts(bulk);ref<-reference
  if(is.character(ref))ref<-tall_load_reference(ref)
  .validate_reference(ref)
  tr<-.truth(scrna,ref);pt<-colnames(bulk);truthpt<-names(tr$pseudobulk)
  if(anyDuplicated(.canonical(pt,ref$aliases))||anyDuplicated(.canonical(truthpt,ref$aliases)))stop('Duplicate biological patients after alias resolution')
  requested_flow<-.prepare_flow(flow,ref,pt)
  ti<-match(.canonical(pt,ref$aliases),.canonical(truthpt,ref$aliases))
  if(anyNA(ti))stop("Paired truth patient missing")
  truthpt<-truthpt[ti]
  overlap<-.canonical(pt,ref$aliases)%in%.canonical(ref$patients,ref$aliases)
  if(split=="independent"&&any(overlap))stop("Reference leakage: validation patients occur in training reference")
  if(dir.exists(out_dir))stop("Benchmark directory exists; use new version")
  dir.create(out_dir,recursive=TRUE)
  refs<-if(split=="holdout")rep(list(.subset_reference(ref,pt)),length(pt))else if(split=="lopo")lapply(pt,function(p).subset_reference(ref,p))else rep(list(ref),length(pt))
  names(refs)<-pt;pred<-list()
  for(p in pt) {
    r<-if(is.null(results))NULL else if(inherits(results,"tall_result")||is.character(results))results else results[[p]]
    imported<-!is.null(r)
    if(!is.null(results)&&is.null(r))stop('Supplied results missing patient: ',p,'; no automatic refitting')
    if(is.null(r)) r<-tall_deconv(bulk[,p,drop=FALSE],refs[[p]],flow=if(is.null(flow))NULL else flow[flow$sample_id==p,,drop=FALSE],config=config,out_dir=file.path(out_dir,sprintf('fit_%04d',match(p,pt))))
    if(is.character(r))r<-tall_read_result(r)
    .validate_fit_counts(r,bulk,p,refs[[p]])
    .validate_fit_scale(r,refs[[p]])
    r$flow<-.fit_flow(r)
    if(imported&&!is.null(flow)&&!isTRUE(all.equal(.flow_signature(requested_flow,p),.flow_signature(r$flow,p),tolerance=1e-12,check.attributes=FALSE)))stop('Imported fit did not apply the requested flow observations; refit explicitly with tall_add_flow')
    members<-.read_tsv(file.path(r$path,"metadata/reference_members.tsv"))$patient_id
    if(any(.canonical(if(split=="holdout")pt else p,ref$aliases)%in%.canonical(members,ref$aliases)))stop("Imported fit has reference leakage")
    expected<-.canonical(refs[[p]]$patients,ref$aliases)
    if(!setequal(.canonical(members,ref$aliases),expected))stop("Imported fit reference members do not match specified split")
    input_path<-file.path(r$path,"metadata/frozen_model_input.rds")
    actual_phi<-if(!is.null(r$reference))r$reference$phi else if(file.exists(input_path))readRDS(input_path)$phi else NULL
    if(is.null(actual_phi)||!identical(actual_phi,refs[[p]]$phi))stop("Imported fit reference parameters cannot be verified or differ from the specified reference")
    pred[[p]]<-r
  }
  composition<-list();expression<-list();patient_matrices<-list();coverage<-list();qc<-list();n<-0L;m<-0L
  allgroups<-unique(c(ref$groups$coarse_group,ref$groups$component_id))
  if(!is.null(groups)){if(!all(groups%in%allgroups))stop("Unknown evaluation group");allgroups<-groups}
  for(i in seq_along(pt)) {
    p<-pt[i];tp<-truthpt[i];rr<-refs[[p]];r<-pred[[p]]
    if(!isTRUE(as.logical(tr$source_present[tp]))) {qc[[p]]<-data.frame(patient=p,truth_status="source_absent",included=FALSE,modeled_RNA_fraction_of_provided_truth=NA_real_);next}
    z<-tr$pseudobulk[[tp]][rr$genes,,drop=FALSE];nc<-tr$cell_counts[tp,]
    model<-intersect(colnames(z),rr$groups$component_id);total<-sum(z[,model,drop=FALSE]);ntotal<-sum(nc[model]);if(total<=0||ntotal<=0)stop("Truth has no modeled mass")
    qc[[p]]<-data.frame(patient=p,truth_status="present",included=TRUE,modeled_RNA_fraction_of_provided_truth=total/sum(z))
    f<-tall_extract(r,"composition_state_rna",samples=p);cp<-tall_extract(r,"composition_state_cell_proxy",samples=p)
    for(g in allgroups) {
      kk<-rr$groups$component_id[rr$groups$component_id==g|rr$groups$coarse_group==g];kk<-intersect(kk,colnames(z))
      if(!length(kk))next
      truth_rna<-sum(z[,kk,drop=FALSE])/total;truth_cell<-sum(nc[kk])/ntotal
      predicted_rna<-if(all(kk%in%colnames(f)))sum(f[p,kk])else NA_real_
      predicted_cell<-if(all(kk%in%colnames(cp)))sum(cp[p,kk])else NA_real_
      fl<-r$flow;inflow<-!is.null(fl)&&any(fl$sample_id==p&fl$used&fl$group==g)
      n<-n+1L;composition[[n]]<-data.frame(patient=p,group=g,truth_cells=sum(nc[kk]),truth_RNA=truth_rna,predicted_RNA=predicted_rna,RNA_error_pp=100*(predicted_rna-truth_rna),truth_cell=truth_cell,predicted_cell_proxy=predicted_cell,cell_proxy_error_pp=100*(predicted_cell-truth_cell),evaluation_role=if(inflow)"input_gate_calibration"else"not_directly_input_group")
      eligible<-sum(nc[kk])>=min_cells
      if(!eligible)next
      tx<-rowSums(z[,kk,drop=FALSE]);if(sum(tx)<=0)next
      truthlog<-log2(1+tx/sum(tx)*1e6)
      px<-tall_extract(r,"expression",g,genes=rr$genes,samples=p)[,1]
      avail<-all(is.finite(px));plog<-log2(1+px)
      base<-if(g%in%names(rr$coarse_baseline))rr$coarse_baseline[[g]]else rr$baseline[g,]
      baseline<-log2(1+base*1e6)
      m<-m+1L;expression[[m]]<-data.frame(patient=p,group=g,truth_cells=sum(nc[kk]),truth_eligible=TRUE,prediction_available=avail,
        RMSE=if(avail)sqrt(mean((plog-truthlog)^2))else NA_real_,Pearson=.cor_safe(plog,truthlog),Spearman=.cor_safe(plog,truthlog,"spearman"),
        reference_mean_RMSE=sqrt(mean((baseline-truthlog)^2)),improvement_over_mean=if(avail)sqrt(mean((baseline-truthlog)^2))-sqrt(mean((plog-truthlog)^2))else NA_real_,composition_gate_input=inflow,evaluation_role="single_cell_expression_truth_not_input")
      if(is.null(patient_matrices[[g]]))patient_matrices[[g]]<-list(truth=list(),predicted=list())
      patient_matrices[[g]]$truth[[p]]<-truthlog;patient_matrices[[g]]$predicted[[p]]<-plog
    }
  }
  comp<-if(length(composition))do.call(rbind,composition)else data.frame()
  expr<-if(length(expression))do.call(rbind,expression)else data.frame()
  diffs<-list()
  for(g in names(patient_matrices)) {
    t<-do.call(cbind,patient_matrices[[g]]$truth);x<-do.call(cbind,patient_matrices[[g]]$predicted)
    good<-colSums(is.finite(x))==nrow(x);t<-t[,good,drop=FALSE];x<-x[,good,drop=FALSE];nn<-ncol(x)
    if(nn>=3L){rho<-vapply(seq_len(nrow(t)),function(j).cor_safe(t[j,],x[j,]),numeric(1));st<-apply(t,1,stats::sd);sx<-apply(x,1,stats::sd);ratio<-sx/st;ratio[st<1e-8]<-NA_real_}else{rho<-ratio<-NA_real_}
    med<-function(z)if(any(is.finite(z)))stats::median(z[is.finite(z)])else NA_real_
    diffs[[g]]<-data.frame(group=g,evaluable_patients=nn,genes_with_patient_correlation=sum(is.finite(rho)),median_patient_Pearson=med(rho),median_patient_SD_ratio=med(ratio))
  }
  diff<-if(length(diffs))do.call(rbind,diffs)else data.frame()
  if(nrow(expr))coverage<-do.call(rbind,lapply(split(expr,expr$group),function(z)data.frame(group=z$group[1],truth_eligible_N=nrow(z),prediction_available_N=sum(z$prediction_available),RMSE_N=sum(is.finite(z$RMSE)),mean_RMSE=if(any(is.finite(z$RMSE)))mean(z$RMSE,na.rm=TRUE)else NA_real_)))else coverage<-data.frame()
  audit<-data.frame(patient=pt,split=split,reference_patients=vapply(refs,function(x)length(x$patients),integer(1)),overlap_after_split=vapply(seq_along(pt),function(i)any(.canonical(if(split=="holdout")pt else pt[i],ref$aliases)%in%.canonical(refs[[i]]$patients,ref$aliases)),logical(1)))
  composition_summary<-if(nrow(comp))do.call(rbind,lapply(split(comp,list(comp$group,comp$evaluation_role),drop=TRUE),function(z){
    good<-is.finite(z$RNA_error_pp);ce<-is.finite(z$cell_proxy_error_pp)
    data.frame(group=z$group[1],evaluation_role=z$evaluation_role[1],RNA_N=sum(good),RNA_MAE_pp=if(any(good))mean(abs(z$RNA_error_pp[good]))else NA_real_,RNA_RMSE_pp=if(any(good))sqrt(mean(z$RNA_error_pp[good]^2))else NA_real_,cell_proxy_N=sum(ce),cell_proxy_MAE_pp=if(any(ce))mean(abs(z$cell_proxy_error_pp[ce]))else NA_real_)
  }))else data.frame()
  truth_coverage<-if(length(qc))do.call(rbind,qc)else data.frame()
  out<-list(composition=comp,composition_summary=composition_summary,expression=expr,coverage=coverage,patient_specific_difference=diff,split_audit=audit,truth_coverage=truth_coverage)
  for(nm in names(out)) .write_tsv(out[[nm]],file.path(out_dir,paste0(nm,".tsv")))
  saveRDS(out,file.path(out_dir,"benchmark.rds"));out
}
