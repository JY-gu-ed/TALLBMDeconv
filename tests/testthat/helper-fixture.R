fixture <- function() {
  x<-matrix(c(20,4,8,3,5,1, 18,5,9,3,4,2, 2,18,1,16,2,8,
              3,17,2,15,3,9, 19,4,8,2,6,1, 2,19,1,15,2,9),nrow=6)
  # Three donors with both cell identities; simulation is a software fixture only.
  colnames(x)<-paste0("cell",seq_len(ncol(x)));rownames(x)<-paste0("G",1:6)
  a<-data.frame(cell_id=colnames(x),group=c("M1","M1","B","B","M1","B"),
    coarse_group=c("Malignant_total","Malignant_total","B","B","Malignant_total","B"),
    compartment=c("malignant","malignant","nonmalignant","nonmalignant","malignant","nonmalignant"))
  donors<-c("D1","D2","D1","D2","D3","D3")
  ref<-tall_build_reference(x,a,donors,"fixture_v1")
  bulk<-cbind(V1=rowSums(x[,c(1,3),drop=FALSE]),V2=rowSums(x[,c(2,4),drop=FALSE]))
  list(x=x,a=a,donors=donors,ref=ref,bulk=bulk,
       cfg=list(iterations=600L,burnin=200L,thin=4L))
}
flow_fixture <- function(sample="V1",kind="point",value=.6) data.frame(sample_id=sample,group="Malignant_total",kind=kind,
  denominator_scope="modeled_nucleated",original_denominator="all modeled nucleated cells",gate_identity="audited synthetic gate",
  value=value,lower=.5,upper=.8,measurement_error=.05)
