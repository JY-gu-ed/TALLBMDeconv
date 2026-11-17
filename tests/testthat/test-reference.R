test_that("custom reference preserves donor identity, names and cell size", {
  f<-fixture();expect_s3_class(f$ref,"tall_reference");expect_equal(rowSums(f$ref$phi),c(M1=1,B=1));expect_true(all(f$ref$cell_sizes>0))
  p<-tempfile(fileext=".rds");saveRDS(f$ref,p);expect_identical(tall_load_reference(p)$reference_hash,f$ref$reference_hash)
  bad<-f$ref;bad$phi[1,1]<-.5;saveRDS(bad,p);expect_error(tall_load_reference(p),"Invalid reference|hash")
  x<-f$x;colnames(x)[2]<-colnames(x)[1];expect_error(tall_build_reference(x,f$a,f$donors),"unique")
})
test_that("installed resources have distinct reference and validation scope", {
  skip_if(!nzchar(system.file("references","benchmark73.rds",package="TALLBMDeconv")),"Authorised atlas resources absent in public-code distribution")
  a<-tall_load_reference(version="benchmark73");b<-tall_load_reference()
  expect_length(a$patients,73);expect_length(b$patients,81)
  expect_false(identical(a$reference_id,b$reference_id));expect_equal(ncol(a$phi),15970)
  expect_match(b$validation_scope,"not.*inherit|not.*independent|no.*independent")
  expect_false(any(TALLBMDeconv:::.canonical(a$patients,a$aliases)%in%a$aliases$bulk_id))
})
test_that("counts and IDs fail loudly rather than being silently corrected", {
  f<-fixture();x<-f$bulk;x[1,1]<-.25
  expect_error(tall_deconv(x,f$ref),"integer")
  expect_error(tall_deconv(f$bulk,f$ref,input_type="TPM"),"raw integer")
  expect_error(tall_deconv(f$bulk[-1,,drop=FALSE],f$ref),"genes missing")
  expect_error(tall_deconv(f$bulk,f$ref,config=list(bogus=1)),"Unknown")
})
