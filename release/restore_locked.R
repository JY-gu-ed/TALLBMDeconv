# Supported entry: fresh official digest-pinned Linux/amd64 container.
setwd('/work');Sys.setenv(TZ='UTC',BIOCONDUCTOR_USE_CONTAINER_REPOSITORY='FALSE')
stopifnot(getRversion()=='4.5.2')
verify_only<-identical(Sys.getenv('TALL_RESTORE_VERIFY_ONLY'),'TRUE')
if(dir.exists('R_library')&&!verify_only)stop('Use a NEW empty workdir; existing library is preserved')
dir.create('source_cache',showWarnings=FALSE);dir.create('tmp',showWarnings=FALSE)
if(!verify_only)dir.create('R_library')
.libPaths(c('/work/R_library',.libPaths()));Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
lock<-read.delim('release/DEPENDENCY_SOURCE_LOCK.tsv',stringsAsFactors=FALSE)
lock<-lock[order(lock$install_order),,drop=FALSE];stopifnot(!anyDuplicated(lock$Package))
hash<-function(p)strsplit(system2('sha256sum',shQuote(p),stdout=TRUE),'[[:space:]]+')[[1]][1]
get_source<-function(url,dest,expected,pkg=NULL){
 if(!file.exists(dest)){
  ok<-try(download.file(url,dest,mode='wb',quiet=FALSE),silent=TRUE)
  if((inherits(ok,'try-error')||ok!=0)&&!is.null(pkg)){
   fallback<-paste0('https://cloud.r-project.org/src/contrib/Archive/',pkg,'/',basename(dest))
   download.file(fallback,dest,mode='wb')
  }
 }
 stopifnot(file.exists(dest),identical(hash(dest),expected))
}
for(i in seq_len(nrow(lock))){
 r<-lock[i,];stopifnot(basename(r$archive)==r$archive)
 p<-file.path('source_cache',r$archive);get_source(r$source_url,p,r$sha256,r$Package)
 if(!verify_only)stopifnot(system2(file.path(R.home('bin'),'R'),c('CMD','INSTALL','-l','/work/R_library',shQuote(p)))==0)
}
bp<-jsonlite::fromJSON('release/BAYESPRISM_SOURCE_LOCK.json')
get_source(bp$url,'BayesPrism_source.tar.gz',bp$source_archive_sha256)
code<-jsonlite::fromJSON('release/CODE_ARCHIVE_LOCK.json');stopifnot(hash('TALLBMDeconv_0.1.1.tar.gz')==code$sha256)
if(!verify_only){
 dir.create('bp_source');untar('BayesPrism_source.tar.gz',exdir='bp_source');p<-file.path(list.files('bp_source',full.names=TRUE),'BayesPrism');stopifnot(length(p)==1)
 stopifnot(system2(file.path(R.home('bin'),'R'),c('CMD','INSTALL','-l','/work/R_library',shQuote(p)))==0)
 stopifnot(system2(file.path(R.home('bin'),'R'),c('CMD','INSTALL','-l','/work/R_library','TALLBMDeconv_0.1.1.tar.gz'))==0)
 actual<-vapply(lock$Package,function(p)as.character(packageVersion(p)),character(1))
 expected<-as.character(package_version(lock$Version))
 matched<-unname(actual)==expected
 write.table(data.frame(Package=lock$Package,source_version=lock$Version,canonical_expected=expected,installed=unname(actual),matched=matched),
             'INSTALLED_VERSION_AUDIT.tsv',sep='\t',quote=FALSE,row.names=FALSE)
 if(!all(matched))stop('Installed dependency version differs from lock: ',paste(lock$Package[!matched],collapse=', '))
}
writeLines(if(verify_only)'PASS_EXISTING_SOURCE_CACHE_HASHES_AND_RESTORE_PARSE'else'PASS_LOCKED_SOURCE_RESTORE','RESTORE_ACCEPTANCE.txt')
