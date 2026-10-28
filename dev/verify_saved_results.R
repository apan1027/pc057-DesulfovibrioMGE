# Read-only comparison of saved RDS components with archived CSV exports.
suppressPackageStartupMessages(library(data.table))
args<-commandArgs(FALSE); filearg<-sub('^--file=','',args[grepl('^--file=',args)])
cli<-commandArgs(TRUE); if(length(cli)!=2 || cli[1]!='--bundle')stop('Use --bundle PREPARED_INPUT_DIRECTORY'); B<-normalizePath(cli[2]); D<-file.path(B,'inputs/data')
r<-readRDS(file.path(B,'inputs/models/revised_results.rds'))
map<-c(meta='metadata_corrected',per_genome='per_genome_corrected',enrichment='enrichment_corrected',cluster_intervals='cluster_intervals_corrected',spatial='pdcs31_gene_position_null',ecology_models='ecology_models_corrected',outgroup='outgroup_corrected',outgroup_tests='outgroup_tests_corrected',outgroup_models='outgroup_models_corrected',shared_species='shared_species_corrected')
normal<-function(x){s<-as.character(x);s[is.na(s)|s=='NA']<-'';s[tolower(s)=='true']<-'TRUE';s[tolower(s)=='false']<-'FALSE';s}
for(name in names(map)){
  p<-file.path(D,paste0(map[[name]],'.csv'))
  if(name=='cluster_intervals')p<-file.path(B,'inputs/previous-export/cluster_intervals_corrected.csv')
  a<-as.data.frame(r[[name]]); b<-as.data.frame(fread(p,na.strings=c('','NA')))
  stopifnot(nrow(a)==nrow(b),all(names(b)%in%names(a)))
  for(col in names(b)){
    x<-a[[col]];y<-b[[col]]
    if(is.numeric(x)&&is.numeric(y)){
      stopifnot(identical(is.na(x),is.na(y)))
      ix<-!is.na(x);stopifnot(all(abs(x[ix]-y[ix])<=1e-10*pmax(1,abs(x[ix]))))
    }else stopifnot(identical(normal(x),normal(y)))
  }
  message(name,': ',nrow(b),' rows / ',ncol(b),' CSV columns agree')
}
message('PASS: all ten saved-result components match their CSV exports; no inferential computation executed.')
