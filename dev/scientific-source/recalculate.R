suppressPackageStartupMessages({library(data.table);library(dplyr);library(tidyr);library(lme4);library(lmerTest);library(MASS);library(jsonlite)})
cli<-commandArgs(TRUE)
if(length(cli)!=5 || cli[1]!='--input' || cli[3]!='--output' || cli[5]!='--allow-refit')stop('Not part of the default workflow. Explicit use: --input INPUT_DIR --output NEW_DIR --allow-refit')
if(file.exists(cli[4]))stop('Refuse to overwrite an existing output directory')
O<-normalizePath(cli[4],mustWork=FALSE);D<-file.path(O,'working/data');dir.create(D,recursive=TRUE)
stopifnot(all(file.copy(list.files(cli[2],full.names=TRUE),D)))
read <- function(n) fread(file.path(D,paste0(n,'.csv')),na.strings=c('','NA'))
save <- function(x,n) fwrite(as.data.table(x),file.path(D,paste0(n,'.csv')),na='')
truth <- function(x) tolower(as.character(x)) %in% c('true','1')
m <- read('metadata_corrected'); m[, species_assigned:=truth(species_assigned)]
t <- read('call_tiers'); g <- read('genome_mge_summary')
meta <- as_tibble(m) %>% left_join(as_tibble(g),by=c('genome_id','accession'='assembly_accession'))
load <- as_tibble(t) %>% count(cohort,category,tier,genome_id,name='n')
pg <- expand_grid(genome_id=m$genome_id,tier=c('T2','T3','T4'),category=c('Defense','Antidefense')) %>%
 left_join(load %>% filter(cohort=='Desulfovibrio') %>% dplyr::select(-cohort),by=c('genome_id','tier','category')) %>% mutate(n=replace_na(n,0L)) %>%
 left_join(meta,by='genome_id') %>% mutate(log_size=log(genome_size_bp/1e6),eco=factor(ecology_category,levels=c('animal_gut','human_gut')))
save(pg,'per_genome_corrected')
save(m[,.(n=.N),by=.(ecology_category)],'cohort_ecology_counts')
shared <- m[ecology_category %in% c('human_gut','animal_gut') & species_assigned,.(n=.N),by=.(gtdb_species,ecology_category)] %>%
 pivot_wider(names_from=ecology_category,values_from=n,values_fill=0) %>% mutate(shared=human_gut>0 & animal_gut>0,balanced=human_gut>=2 & animal_gut>=2)
save(shared,'shared_species_corrected')
nest <- m[ecology_category %in% c('human_gut','animal_gut') & species_assigned,.(n=.N,n_ecology=uniqueN(ecology_category),labels=paste(sort(unique(ecology_category)),collapse=';')),by=bioproject]
save(nest,'bioproject_nesting_corrected')
message('Fitting corrected ecology models')
fits <- list(); fitlog <- list(); counter<-0
record <- function(fit,label,tier,category,n,subset='full_assigned_gut',scale='difference') {
 co <- coef(summary(fit)); cn <- if('ecohuman_gut'%in%rownames(co))'ecohuman_gut' else 'cohortOutgroup'
 if(!cn%in%rownames(co))stop('Missing contrast in ',label)
 est<-co[cn,1];se<-co[cn,2];pcol<-grep('Pr\\(',colnames(co),value=TRUE)
 pv<-if(length(pcol))co[cn,pcol[1]] else 2*pnorm(-abs(est/se))
 data.frame(model=label,tier=tier,category=category,subset=subset,n=n,estimate=est,se=se,ci_lo=est-1.96*se,ci_hi=est+1.96*se,p=pv,scale=scale,singular=if(inherits(fit,'merMod'))isSingular(fit) else FALSE)
}
safe <- function(fun,label) tryCatch(fun(),error=function(e){message(label,': ',conditionMessage(e));NULL})
for(cat in c('Defense','Antidefense'))for(tier0 in c('T2','T3','T4')){
 dd <- pg %>% filter(tier==tier0,category==cat,!is.na(eco),species_assigned) %>% as.data.frame()
 dd$gtdb_species<-factor(dd$gtdb_species);dd$bioproject<-factor(dd$bioproject)
 # Baseline model is retained as an association model, with study adjustment reported alongside.
 settings<-list(species_random=dd,clone_representatives=dd %>% group_by(clone_cluster) %>% arrange(desc(completeness),contamination,accession) %>% slice_head(n=1) %>% ungroup(),one_per_project_species=dd %>% group_by(bioproject,gtdb_species) %>% arrange(desc(completeness),contamination,accession) %>% slice_head(n=1) %>% ungroup())
 for(nm in names(settings)){
  ds<-as.data.frame(settings[[nm]])
  ff<-safe(function()lmerTest::lmer(n~eco+log_size+(1|gtdb_species),data=ds,REML=TRUE,control=lmerControl(check.conv.singular='ignore')),nm)
  if(!is.null(ff)){counter<-counter+1;fits[[counter]]<-record(ff,nm,tier0,cat,nrow(ds));fitlog[[paste(cat,tier0,nm)]]<-ff}
 }
 ff<-safe(function()lmerTest::lmer(n~eco+log_size+(1|gtdb_species)+(1|bioproject),data=dd,REML=TRUE,control=lmerControl(check.conv.singular='ignore')),'species_project_random')
 if(!is.null(ff)){counter<-counter+1;fits[[counter]]<-record(ff,'species_project_random',tier0,cat,nrow(dd));fitlog[[paste(cat,tier0,'study')]]<-ff}
 # NB is fit with lme4; glmmTMB is deliberately not used because of its local TMB mismatch.
 ff<-safe(function()glmer(n~eco+log_size+(1|gtdb_species),data=dd,family=poisson,control=glmerControl(optimizer='bobyqa',optCtrl=list(maxfun=50000))),'poisson')
 if(!is.null(ff)){counter<-counter+1;fits[[counter]]<-record(ff,'poisson_species',tier0,cat,nrow(dd),scale='log_IRR');fitlog[[paste(cat,tier0,'poisson')]]<-ff}
 if(cat=='Defense'){
  ff<-safe(function()glmer.nb(n~eco+log_size+(1|gtdb_species),data=dd,control=glmerControl(optimizer='bobyqa',optCtrl=list(maxfun=50000))),'negative_binomial')
  if(!is.null(ff)){counter<-counter+1;fits[[counter]]<-record(ff,'NB_species',tier0,cat,nrow(dd),scale='log_IRR');fitlog[[paste(cat,tier0,'NB')]]<-ff}
 }
}
ef<-bind_rows(fits);save(ef,'ecology_models_corrected');saveRDS(fitlog,file.path(O,'working/model_fits.rds'))
diagnostics<-bind_rows(lapply(names(fitlog),function(key){
 f<-fitlog[[key]];parts<-strsplit(key,' ',fixed=TRUE)[[1]]
 theta<-if(parts[3]=='NB')getME(f,'glmer.nb.theta') else NA_real_
 data.frame(fit_key=key,category=parts[1],tier=parts[2],model_key=parts[3],optimizer_message=paste(f@optinfo$conv$lme4$messages,collapse=';'),singular=isSingular(f),nb_theta=theta,poisson_boundary=!is.na(theta) && theta>1e5)
}));save(diagnostics,'model_diagnostics')

message('Recomputing enrichment under separated categories and common name normalization')
tests<-list();boots<-list();bt<-0;bb<-0
set.seed(20261020);B<-4000L
# Exact sum of independent finite-population hypergeometric draws, conditional on each stratum's margins.
conditional <- function(N,K,n,x) {
 ok<-n>0 & N>0;N<-N[ok];K<-K[ok];n<-n[ok]
 dist<-1
 for(i in seq_along(N)){
  vals<-0:min(n[i],K[i]);q<-dhyper(vals,K[i],N[i]-K[i],n[i]);out<-numeric(length(dist)+length(q)-1)
  for(j in seq_along(q))out[seq_along(dist)+j-1]<-out[seq_along(dist)+j-1]+dist*q[j]
  dist<-out
 }
 dist<-pmax(dist,0);dist<-dist/sum(dist);v<-seq_along(dist)-1;mu<-sum(n*K/N)
 list(expected=mu,p=sum(dist[abs(v-mu)>=abs(x-mu)-1e-10]),p_enrich=sum(dist[v>=x]),p_deplete=sum(dist[v<=x]))
}
for(cat in c('Defense','Antidefense'))for(tier0 in c('T2','T3')){
 z<-t[cohort=='Desulfovibrio' & category==cat & tier==tier0]
 base<-mean(z$near);sizes<-z[,.(total=.N,n_near=sum(near),n_genomes=uniqueN(genome_id)),by=cname][total>=10]
 for(ii in seq_len(nrow(sizes))){
  typ<-sizes$cname[ii];a<-z[,.(N=.N,K=sum(near),n=sum(cname==typ)),by=genome_id]
  v<-conditional(a$N,a$K,a$n,sizes$n_near[ii])
  # Clone representative sensitivity uses representatives chosen without conditioning on subtype presence.
  reps<-m[order(-completeness,contamination,accession),.SD[1],by=clone_cluster]$genome_id
  zz<-z[genome_id%in%reps];ac<-zz[,.(N=.N,K=sum(near),n=sum(cname==typ)),by=genome_id];xc<-sum(zz$cname==typ & zz$near==1);vc<-conditional(ac$N,ac$K,ac$n,xc)
  bt<-bt+1;tests[[bt]]<-data.frame(category=cat,tier=tier0,cname=typ,total=sizes$total[ii],n_near=sizes$n_near[ii],n_genomes=sizes$n_genomes[ii],baseline=base,fraction=sizes$n_near[ii]/sizes$total[ii],p_binomial=binom.test(sizes$n_near[ii],sizes$total[ii],base)$p.value,conditional_expected=v$expected,p_conditional=v$p,clone_total=sum(zz$cname==typ),clone_near=xc,clone_expected=vc$expected,p_conditional_clone=vc$p)
 }
 # Cluster intervals are paired differences on the SAME counting convention; no sign-tail bootstrap P is manufactured.
 for(cluster0 in c('genome','species')){
  zz<-merge(z,m[,.(genome_id,gtdb_species)],by='genome_id',all.x=TRUE)
  if(cluster0=='species')zz<-zz[!is.na(gtdb_species)]
  zz[,cluster:=if(cluster0=='species')gtdb_species else genome_id]
  cls<-sort(unique(zz$cluster));types<-sizes$cname
  totals<-zz[,.(N=.N,K=sum(near)),by=cluster][match(cls,cluster)]
  subN<-matrix(0,length(cls),length(types),dimnames=list(cls,types));subK<-subN
  a<-zz[cname%in%types,.(N=.N,K=sum(near)),by=.(cluster,cname)]
  subN[cbind(match(a$cluster,cls),match(a$cname,types))]<-a$N;subK[cbind(match(a$cluster,cls),match(a$cname,types))]<-a$K
  w<-t(rmultinom(B,length(cls),rep(1/length(cls),length(cls))))
  baseb<-as.vector(w%*%totals$K)/as.vector(w%*%totals$N)
  ns<-w%*%subN;ks<-w%*%subK;delta<-ks/ns-baseb
  for(j in seq_along(types)){
   obsN<-sum(subN[,j]);obsK<-sum(subK[,j]);v<-delta[,j];v<-v[is.finite(v)]
   degenerate<-obsK==0 || obsK==obsN
   bb<-bb+1;boots[[bb]]<-data.frame(category=cat,tier=tier0,cname=types[j],cluster_unit=cluster0,scope=if(cluster0=='species')'assigned_species_activity_carriers' else 'activity_carrier_genomes',n_clusters=length(cls),n_units=nrow(zz),subtype_n=obsN,subtype_near=obsK,baseline=sum(totals$K)/sum(totals$N),delta=obsK/obsN-sum(totals$K)/sum(totals$N),ci_lo=if(degenerate)NA_real_ else unname(quantile(v,.025)),ci_hi=if(degenerate)NA_real_ else unname(quantile(v,.975)),n_boot=length(v),zero_or_all_event=degenerate)
  }
 }
}
en<-bind_rows(tests)%>%group_by(category,tier)%>%mutate(q_binomial=p.adjust(p_binomial,'BH'),q_conditional=p.adjust(p_conditional,'BH'),q_conditional_clone=p.adjust(p_conditional_clone,'BH'),direction=ifelse(fraction>baseline,'enriched','depleted'))%>%ungroup()
save(en,'enrichment_corrected');save(bind_rows(boots),'cluster_intervals_corrected')

message('PDC-S31 gene-position control')
feat<-read('pdcs31_contig_features');feat<-feat[tolower(feature_type)=='cds'];feat[,len:=stop-start+1]
pdc<-t[cohort=='Desulfovibrio' & category=='Defense' & tier=='T3' & cname=='pdcs31'];ph<-read('prophage_consensus')
relation<-function(chrom,s,e,dist){ctg0<-chrom;p<-ph[ph$chrom==ctg0];if(!nrow(p))return(FALSE);any(s<=p$end+dist & e>=p$start-dist)}
stopifnot(sum(mapply(relation,pdc$chrom,pdc$start,pdc$stop,MoreArgs=list(dist=0)))==sum(pdc$overlap),sum(mapply(relation,pdc$chrom,pdc$start,pdc$stop,MoreArgs=list(dist=10000)))==sum(pdc$near))
spatial<-list();sn<-0;pools<-list()
for(scope0 in c('all_loci','clone_representatives')){
 reps<-m[order(-completeness,contamination,accession),.SD[1],by=clone_cluster]$genome_id
 pp<-if(scope0=='all_loci')pdc else pdc[genome_id%in%reps]
 for(distance0 in c(0,5000,10000)){
  strat<-list();ii<-0
  for(chrom0 in unique(pp$chrom)){
   p<-pp[chrom==chrom0];lens<-p$stop-p$start+1
   f<-feat[chrom==chrom0 & len>=min(lens)*.5 & len<=max(lens)*2]
   k<-sum(mapply(relation,f$chrom,f$start,f$stop,MoreArgs=list(dist=distance0)))
   x<-sum(mapply(relation,p$chrom,p$start,p$stop,MoreArgs=list(dist=distance0)))
   stopifnot(nrow(f)>=nrow(p));ii<-ii+1;strat[[ii]]<-data.frame(chrom=chrom0,N=nrow(f),K=k,n=nrow(p),x=x)
  }
  ss<-bind_rows(strat);v<-conditional(ss$N,ss$K,ss$n,sum(ss$x));sn<-sn+1
  spatial[[sn]]<-data.frame(scope=scope0,distance_bp=distance0,n_loci=sum(ss$n),observed=sum(ss$x),expected=v$expected,p=v$p,p_enrichment=v$p_enrich)
  pools[[sn]]<-ss %>% mutate(scope=scope0,distance_bp=distance0) %>% dplyr::select(scope,distance_bp,everything())
 }
}
sp<-bind_rows(spatial)%>%group_by(scope)%>%mutate(q=p.adjust(p,'BH'))%>%ungroup();save(sp,'pdcs31_gene_position_null');save(bind_rows(pools),'pdcs31_gene_position_strata')

message('Matched outgroup recalculation')
og<-read('outgroup_per_genome');og<-as_tibble(og)
og<-og%>%dplyr::select(-defense,-defense_per_mb)%>%left_join(meta%>%dplyr::select(genome_id,ecology_category,contig_n50,gtdb_species,bioproject),by='genome_id')
outload<-as_tibble(t)%>%filter(category=='Defense')%>%count(genome_id,tier,name='defense')
ogall<-expand_grid(genome_id=og$genome_id,tier=c('T2','T3','T4'))%>%left_join(og,by='genome_id')%>%left_join(outload,by=c('genome_id','tier'))%>%mutate(defense=replace_na(defense,0),defense_per_mb=defense/genome_mb,cohort=factor(cohort,levels=c('Desulfovibrio','Outgroup')),qc_pass=truth(qc_pass))
save(ogall,'outgroup_corrected')
ot<-list();om<-list();on<-0;mn<-0
for(scope0 in c('all_core_vs_qc','confirmed_gut_vs_qc'))for(tier0 in c('T2','T3','T4')){
 ds<-ogall%>%filter(tier==tier0,qc_pass,cohort=='Outgroup' | scope0=='all_core_vs_qc' | ecology_category%in%c('human_gut','animal_gut'))
 a<-ds%>%filter(cohort=='Desulfovibrio');b<-ds%>%filter(cohort=='Outgroup')
 for(y in c('provirus_per_mb','defense_per_mb','provirus','defense')){
  test<-wilcox.test(a[[y]],b[[y]],exact=FALSE);on<-on+1
  ot[[on]]<-data.frame(scope=scope0,tier=tier0,outcome=y,n_core=nrow(a),n_outgroup=nrow(b),core_median=median(a[[y]]),outgroup_median=median(b[[y]]),p=test$p.value)
 }
 for(k in c(1,2,3,5)){
  tab<-matrix(c(sum(a$provirus>=k),sum(a$provirus<k),sum(b$provirus>=k),sum(b$provirus<k)),nrow=2,byrow=TRUE)
  ft<-fisher.test(tab);on<-on+1;ot[[on]]<-data.frame(scope=scope0,tier=tier0,outcome=paste0('carriage_ge',k),n_core=nrow(a),n_outgroup=nrow(b),core_median=mean(a$provirus>=k),outgroup_median=mean(b$provirus>=k),p=ft$p.value)
 }
 ds$log_contigs<-log(ds$n_contigs_3kb);ds$log_mb<-log(ds$genome_mb)
 ff<-glm(I(provirus>=1)~cohort+log_mb+log_contigs,data=ds,family=binomial)
 mn<-mn+1;om[[mn]]<-record(ff,'carriage_logistic_contiguity',tier0,'Provirus',nrow(ds),scope0,'log_OR')
 for(outcome in c('provirus','defense')){
  ff<-safe(function()MASS::glm.nb(as.formula(paste(outcome,'~cohort+log_contigs+offset(log_mb)')),data=ds),'outgroup_NB')
  if(!is.null(ff)){mn<-mn+1;om[[mn]]<-record(ff,paste0(outcome,'_NB_contiguity'),tier0,outcome,nrow(ds),scope0,'log_IRR')}
 }
}
save(bind_rows(ot),'outgroup_tests_corrected');save(bind_rows(om),'outgroup_models_corrected')
saveRDS(list(meta=meta,per_genome=pg,enrichment=en,cluster_intervals=bind_rows(boots),spatial=sp,ecology_models=ef,outgroup=ogall,outgroup_tests=bind_rows(ot),outgroup_models=bind_rows(om),shared_species=shared),file.path(O,'working/revised_results.rds'))
write_json(list(assigned_species=length(unique(na.omit(m$gtdb_species))),unassigned_genomes=sum(is.na(m$gtdb_species)),ecology=as.list(table(m$ecology_category)),gut_assigned=sum(m$ecology_category%in%c('human_gut','animal_gut') & m$species_assigned),shared_species=sum(shared$shared),balanced_species=sum(shared$balanced),bioprojects_gut=nrow(nest),projects_with_both_labels=sum(nest$n_ecology>1),host_defense_T3=sum(t$cohort=='Desulfovibrio' & t$category=='Defense' & t$tier=='T3'),bootstrap_draws=B),file.path(D,'analysis_summary.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(O,'working/R-session-info.txt'))
message('Recalculation complete')
