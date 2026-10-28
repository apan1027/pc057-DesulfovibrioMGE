# Plot expressions inherited from the validated frozen redraw implementation.
# Only selection, shared setup and output plumbing were refactored.

pc057_plot_fig1_cohort <- function(bundle, output) {
 B<-normalizePath(bundle);D<-file.path(B,"inputs/data");F<-output
 set.seed(20261022)
r<-readRDS(file.path(B,'inputs/models/revised_results.rds'));m<-r$meta;pg<-r$per_genome
theme_set(theme_bw(base_size=11.5,base_family='Helvetica')+theme(panel.grid.minor=element_blank(),plot.tag=element_text(face='bold'),legend.position='bottom',strip.background=element_rect(fill='grey96'),plot.margin=margin(5,9,5,5)))
pal<-c(animal_gut='#0072B2',human_gut='#D55E00');gl<-c(animal_gut='Animal gut',human_gut='Human gut')
savefig<-function(p,name,w,h){ggsave(file.path(F,paste0(name,'.png')),p,width=w,height=h,dpi=400,bg='white');ggsave(file.path(F,paste0(name,'.pdf')),p,width=w,height=h,bg='white')}
sp<-m%>%count(gtdb_species,sort=TRUE)%>%mutate(gtdb_species=ifelse(is.na(gtdb_species),'Species unassigned',sub('Desulfovibrio ','D. ',gtdb_species)))%>%slice_head(n=11)
e<-m%>%count(ecology_category,sort=TRUE)%>%mutate(label=gsub('_',' ',ecology_category))
p1<-ggplot(sp,aes(n,reorder(gtdb_species,n)))+geom_col(fill='#0072B2',width=.7)+geom_text(aes(label=n),hjust=-.2,size=2.7)+scale_x_continuous(expand=expansion(mult=c(0,.12)))+labs(x='Genomes',y=NULL)
p2<-ggplot(e,aes(n,reorder(label,n)))+geom_col(fill='#009E73',width=.7)+geom_text(aes(label=n),hjust=-.2,size=2.7)+scale_x_continuous(expand=expansion(mult=c(0,.14)))+labs(x='Genomes',y=NULL)+theme(axis.text.y=element_text(size=9.5))
p3<-ggplot(m,aes(completeness,contamination,colour=log10(num_contigs)))+geom_point(size=1.3,alpha=.75)+geom_vline(xintercept=90,lty=2,colour='grey55')+geom_hline(yintercept=5,lty=2,colour='grey55')+scale_colour_viridis_c(name='log10 contigs')+labs(x='CheckM2 completeness (%)',y='Contamination (%)')
ct<-m%>%filter(!is.na(country),nzchar(country))%>%count(country,sort=TRUE)%>%slice_head(n=8)
p4<-ggplot(ct,aes(n,reorder(country,n)))+geom_col(fill='#56B4E9',width=.7)+geom_text(aes(label=n),hjust=-.2,size=2.7)+scale_x_continuous(expand=expansion(mult=c(0,.14)))+labs(x='Genomes',y=NULL)
savefig((p1|p2)/(p3|p4)+plot_annotation(tag_levels='A'),'fig1_cohort',7.6,6.4)


 invisible(list.files(output,full.names=TRUE))
}

pc057_plot_fig2_catalog <- function(bundle, output) {
 B<-normalizePath(bundle);D<-file.path(B,"inputs/data");F<-output
 set.seed(20261022)
r<-readRDS(file.path(B,'inputs/models/revised_results.rds'));m<-r$meta;pg<-r$per_genome
theme_set(theme_bw(base_size=11.5,base_family='Helvetica')+theme(panel.grid.minor=element_blank(),plot.tag=element_text(face='bold'),legend.position='bottom',strip.background=element_rect(fill='grey96'),plot.margin=margin(5,9,5,5)))
pal<-c(animal_gut='#0072B2',human_gut='#D55E00');gl<-c(animal_gut='Animal gut',human_gut='Human gut')
savefig<-function(p,name,w,h){ggsave(file.path(F,paste0(name,'.png')),p,width=w,height=h,dpi=400,bg='white');ggsave(file.path(F,paste0(name,'.pdf')),p,width=w,height=h,bg='white')}
z<-pg%>%filter(tier=='T3',ecology_category%in%names(pal))%>%mutate(group=ifelse(category=='Defense','Defense candidates','Anti-defense candidates'))
p1<-ggplot(z,aes(ecology_category,n,fill=ecology_category))+geom_violin(alpha=.5,trim=TRUE,colour=NA)+geom_boxplot(width=.18,fill='white',outlier.size=.4)+facet_wrap(~group,scales='free_y')+scale_fill_manual(values=pal,guide='none')+scale_x_discrete(labels=gl)+labs(x=NULL,y='Merged candidates per genome')
z2<-pg%>%filter(tier=='T3',category=='Defense')
p2<-ggplot(z2,aes(n_consensus_prophage,n))+geom_jitter(width=.08,height=0,colour='grey35',alpha=.35,size=1)+labs(x='Consensus prophages per genome',y='Merged defense candidates')
savefig(p1/p2+plot_layout(heights=c(1,1))+plot_annotation(tag_levels='A'),'fig2_catalog',7.2,5.6)


 invisible(list.files(output,full.names=TRUE))
}

pc057_plot_fig3_enrichment <- function(bundle, output) {
 B<-normalizePath(bundle);D<-file.path(B,"inputs/data");F<-output
 set.seed(20261022)
r<-readRDS(file.path(B,'inputs/models/revised_results.rds'));m<-r$meta;pg<-r$per_genome
theme_set(theme_bw(base_size=11.5,base_family='Helvetica')+theme(panel.grid.minor=element_blank(),plot.tag=element_text(face='bold'),legend.position='bottom',strip.background=element_rect(fill='grey96'),plot.margin=margin(5,9,5,5)))
pal<-c(animal_gut='#0072B2',human_gut='#D55E00');gl<-c(animal_gut='Animal gut',human_gut='Human gut')
savefig<-function(p,name,w,h){ggsave(file.path(F,paste0(name,'.png')),p,width=w,height=h,dpi=400,bg='white');ggsave(file.path(F,paste0(name,'.pdf')),p,width=w,height=h,bg='white')}
labels<-c(pdcs31='PDC-S31',pdlambda1='PD-Lambda-1',pdcm11='PDC-M11',pdcs07='PDC-S07',rm='RM',cas='Cas')
en<-r$enrichment%>%filter(category=='Defense',tier=='T3',cname%in%names(labels))%>%mutate(label=factor(labels[cname],levels=rev(labels)),significant=q_conditional<.05,direction=ifelse(fraction>baseline,'Enriched','Depleted'))
p1<-ggplot(en,aes(fraction*100,label,colour=direction,shape=significant))+geom_vline(xintercept=100*en$baseline[1],lty=2,colour='grey45')+geom_point(size=2.8)+geom_text(aes(label=paste0(n_near,'/',total)),nudge_x=3,hjust=0,size=3,show.legend=FALSE)+scale_shape_manual(values=c('TRUE'=16,'FALSE'=1),labels=c('TRUE'='q < 0.05','FALSE'='q >= 0.05'),name='Within-genome test')+scale_colour_manual(values=c(Enriched='#D55E00',Depleted='#0072B2'),guide='none')+scale_x_continuous(limits=c(-1,82))+labs(x='Candidates in prophage context (%)',y=NULL)
ci<-r$cluster_intervals%>%filter(category=='Defense',tier=='T3',cname%in%c('pdcs31','pdlambda1','pdcm11'))%>%mutate(label=factor(labels[cname],levels=rev(labels[c('pdcs31','pdlambda1','pdcm11')])),unit=ifelse(cluster_unit=='genome','Genome clusters','Species clusters'))
p2<-ggplot(ci,aes(delta,label,colour=unit))+geom_vline(xintercept=0,lty=2,colour='grey45')+geom_errorbar(aes(xmin=ci_lo,xmax=ci_hi),position=position_dodge(width=.4),width=.15,orientation='y')+geom_point(position=position_dodge(width=.4),size=2)+scale_colour_manual(values=c('Genome clusters'='#0072B2','Species clusters'='#D55E00'),name=NULL)+labs(x='Subtype fraction minus background',y=NULL)
savefig((p1|p2)+plot_annotation(tag_levels='A'),'fig3_enrichment',7.7,3.8)


 invisible(list.files(output,full.names=TRUE))
}

pc057_plot_fig4_neighborhood <- function(bundle, output) {
 B<-normalizePath(bundle);D<-file.path(B,"inputs/data");F<-output
 set.seed(20261022)
r<-readRDS(file.path(B,'inputs/models/revised_results.rds'));m<-r$meta;pg<-r$per_genome
theme_set(theme_bw(base_size=11.5,base_family='Helvetica')+theme(panel.grid.minor=element_blank(),plot.tag=element_text(face='bold'),legend.position='bottom',strip.background=element_rect(fill='grey96'),plot.margin=margin(5,9,5,5)))
pal<-c(animal_gut='#0072B2',human_gut='#D55E00');gl<-c(animal_gut='Animal gut',human_gut='Human gut')
savefig<-function(p,name,w,h){ggsave(file.path(F,paste0(name,'.png')),p,width=w,height=h,dpi=400,bg='white');ggsave(file.path(F,paste0(name,'.pdf')),p,width=w,height=h,bg='white')}
t<-fread(file.path(D,'call_tiers.csv'));f<-fread(file.path(D,'pdcs31_contig_features.csv'));ph<-fread(file.path(D,'prophage_consensus.csv'));cm<-fread(file.path(D,'contig_map.csv'))
pdc<-t[cohort=='Desulfovibrio' & category=='Defense' & tier=='T3' & cname=='pdcs31']%>%as_tibble()%>%left_join(m%>%dplyr::select(genome_id,accession,gtdb_species,strain),by='genome_id')
pick<-function(species,context){v<-pdc%>%filter(gtdb_species==species,context_class==context)%>%arrange(genome_id);stopifnot(nrow(v)>0);v[1,]}
rep<-bind_rows(pick('Desulfovibrio sp910584445','contained'),pick('Desulfovibrio sp910584445','flank_5kb'),pick('Desulfovibrio desulfuricans_E','contained'),pick('Desulfovibrio sp910584445','no_context'))
save_rows<-list();panels<-list();gp<-c('PDC-S31'='#D55E00','Integrase/recombinase'='#E69F00','Phage hallmark'='#0072B2','Other CDS'='grey65')
for(i in seq_len(nrow(rep))){
 x<-rep[i,];mid<-(x$start+x$stop)/2;ctg<-x$chrom
 ff<-as_tibble(f[chrom==ctg & feature_type=='cds'])%>%filter(stop>=mid-15000,start<=mid+15000)%>%mutate(function_group=case_when(start<=x$stop & stop>=x$start~'PDC-S31',grepl('integrase|recombinase',product,ignore.case=TRUE)~'Integrase/recombinase',grepl('phage|terminase|capsid|tail protein|portal',product,ignore.case=TRUE)~'Phage hallmark',TRUE~'Other CDS'),x0=(ifelse(strand=='+',start,stop)-mid)/1000,x1=(ifelse(strand=='+',stop,start)-mid)/1000)
 pp<-as_tibble(ph[chrom==ctg])%>%filter(end>=mid-15000,start<=mid+15000)%>%mutate(xmin=pmax((start-mid)/1000,-15),xmax=pmin((end-mid)/1000,15))
 title<-paste0(LETTERS[i],'  ',sub('Desulfovibrio ','D. ',x$gtdb_species),'  ',x$accession,'  [',gsub('_',' ',x$context_class),']')
 p<-ggplot()+geom_rect(data=pp,aes(xmin=xmin,xmax=xmax,ymin=-.28,ymax=.28),fill='#56B4E9',alpha=.2)+geom_segment(data=ff,aes(x=x0,xend=x1,y=0,yend=0,colour=function_group),linewidth=2.1,arrow=arrow(type='closed',length=unit(1.2,'mm')))+geom_vline(xintercept=0,lty=3,colour='grey35')+scale_colour_manual(values=gp,name=NULL)+coord_cartesian(xlim=c(-15,15),ylim=c(-.35,.4))+labs(x=if(i==4)'Position relative to PDC-S31 midpoint (kb)' else NULL,y=NULL,title=title)+theme(panel.grid=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(size=9),plot.title=element_text(size=10.5,face='bold'),legend.position='none')
 if(i!=4)p<-p+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
 panels[[i]]<-p;save_rows[[i]]<-ff%>%mutate(panel=LETTERS[i],candidate_accession=x$accession,candidate_chrom=ctg)
}
sp<-r$spatial%>%filter(scope=='all_loci')%>%dplyr::select(distance_bp,observed,expected,q)%>%pivot_longer(c(observed,expected),names_to='type',values_to='n')%>%mutate(distance=factor(distance_bp,levels=c(0,5000,10000),labels=c('Overlap','Within 5 kb','Within 10 kb')))
pp<-ggplot(sp,aes(distance,n,fill=type))+geom_col(position='dodge',width=.65)+geom_text(aes(label=sprintf('%.1f',n)),position=position_dodge(width=.65),vjust=-.3,size=3)+scale_fill_manual(values=c(observed='#D55E00',expected='grey65'),labels=c(expected='Expected from same-contig CDS positions',observed='Observed PDC-S31'),name=NULL)+scale_y_continuous(limits=c(0,34))+labs(x=NULL,y='Loci out of 42',title='E  Spatial comparison conditioned on the candidate-bearing contigs')+theme(plot.title=element_text(size=9,face='bold'),legend.text=element_text(size=8))
savefig(wrap_plots(c(panels,list(pp)),ncol=1,heights=c(1,1,1,1,2)),'fig4_neighborhood',7.7,7.2)
# Input tables are read-only: no neighborhood tables are rewritten.


 invisible(list.files(output,full.names=TRUE))
}

pc057_plot_fig5_outgroup <- function(bundle, output) {
 B<-normalizePath(bundle);D<-file.path(B,"inputs/data");F<-output
 set.seed(20261022)
r<-readRDS(file.path(B,'inputs/models/revised_results.rds'));m<-r$meta;pg<-r$per_genome
theme_set(theme_bw(base_size=11.5,base_family='Helvetica')+theme(panel.grid.minor=element_blank(),plot.tag=element_text(face='bold'),legend.position='bottom',strip.background=element_rect(fill='grey96'),plot.margin=margin(5,9,5,5)))
pal<-c(animal_gut='#0072B2',human_gut='#D55E00');gl<-c(animal_gut='Animal gut',human_gut='Human gut')
savefig<-function(p,name,w,h){ggsave(file.path(F,paste0(name,'.png')),p,width=w,height=h,dpi=400,bg='white');ggsave(file.path(F,paste0(name,'.pdf')),p,width=w,height=h,bg='white')}
og<-r$outgroup%>%filter(tier=='T3',qc_pass,cohort=='Outgroup' | ecology_category%in%c('human_gut','animal_gut'))
op<-c(Desulfovibrio='#0072B2',Outgroup='#D55E00');ol<-c(Desulfovibrio='Desulfovibrio gut subset',Outgroup='Selected gut outgroup')
cv<-bind_rows(lapply(1:5,function(k)og%>%group_by(cohort)%>%summarise(k=k,fraction=mean(provirus>=k),.groups='drop')))
p1<-ggplot(cv,aes(k,fraction,colour=cohort))+geom_line()+geom_point(size=2)+scale_colour_manual(values=op,labels=ol,name=NULL)+scale_x_continuous(breaks=1:5,labels=paste0('>=',1:5))+scale_y_continuous(labels=percent)+labs(x='Proviruses per genome',y='Genomes carrying')
box<-function(y,lab)ggplot(og,aes(cohort,.data[[y]],colour=cohort,fill=cohort))+geom_boxplot(width=.45,alpha=.15,outlier.shape=NA)+geom_jitter(width=.15,size=.5,alpha=.35)+scale_colour_manual(values=op,labels=ol,name=NULL)+scale_fill_manual(values=op,guide='none')+labs(x=NULL,y=lab)+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position='none')
savefig((p1|box('provirus_per_mb','Proviruses per Mb')|box('defense_per_mb','Merged defense candidates per Mb'))+plot_layout(guides='collect')+plot_annotation(tag_levels='A'),'fig5_outgroup',7.7,3.6)


 invisible(list.files(output,full.names=TRUE))
}

pc057_plot_figS1_ecology <- function(bundle, output) {
 B<-normalizePath(bundle);D<-file.path(B,"inputs/data");F<-output
 set.seed(20261022)
r<-readRDS(file.path(B,'inputs/models/revised_results.rds'));m<-r$meta;pg<-r$per_genome
theme_set(theme_bw(base_size=11.5,base_family='Helvetica')+theme(panel.grid.minor=element_blank(),plot.tag=element_text(face='bold'),legend.position='bottom',strip.background=element_rect(fill='grey96'),plot.margin=margin(5,9,5,5)))
pal<-c(animal_gut='#0072B2',human_gut='#D55E00');gl<-c(animal_gut='Animal gut',human_gut='Human gut')
savefig<-function(p,name,w,h){ggsave(file.path(F,paste0(name,'.png')),p,width=w,height=h,dpi=400,bg='white');ggsave(file.path(F,paste0(name,'.pdf')),p,width=w,height=h,bg='white')}
em<-r$ecology_models%>%filter(category=='Defense',scale=='difference')%>%mutate(model=factor(model,levels=rev(c('species_random','clone_representatives','one_per_project_species','species_project_random')),labels=rev(c('Species adjustment','Clone representatives','One per project and species','Species and study adjustment'))),significance=p<.05)
p1<-ggplot(em,aes(estimate,model,colour=significance))+geom_vline(xintercept=0,lty=2,colour='grey40')+geom_errorbar(aes(xmin=ci_lo,xmax=ci_hi),orientation='y',width=.15)+geom_point(size=2)+facet_wrap(~tier,nrow=1)+scale_colour_manual(values=c('TRUE'='#D55E00','FALSE'='grey50'),labels=c('TRUE'='P < 0.05','FALSE'='P >= 0.05'),name=NULL)+labs(x='Human minus animal defense candidates',y=NULL)
savefig(p1,'figS1_ecology',7.7,3.9)

 invisible(list.files(output,full.names=TRUE))
}

pc057_plot_figS2_terl <- function(bundle, output) {
 B<-normalizePath(bundle);D<-file.path(B,"inputs/data");F<-output
 set.seed(20261022)
r<-readRDS(file.path(B,'inputs/models/revised_results.rds'));m<-r$meta;pg<-r$per_genome
theme_set(theme_bw(base_size=11.5,base_family='Helvetica')+theme(panel.grid.minor=element_blank(),plot.tag=element_text(face='bold'),legend.position='bottom',strip.background=element_rect(fill='grey96'),plot.margin=margin(5,9,5,5)))
pal<-c(animal_gut='#0072B2',human_gut='#D55E00');gl<-c(animal_gut='Animal gut',human_gut='Human gut')
savefig<-function(p,name,w,h){ggsave(file.path(F,paste0(name,'.png')),p,width=w,height=h,dpi=400,bg='white');ggsave(file.path(F,paste0(name,'.pdf')),p,width=w,height=h,bg='white')}
file.copy(file.path(B,'frozen-figures/figS2_terl.png'),file.path(F,'figS2_terl.png'),overwrite=TRUE)

 invisible(list.files(output,full.names=TRUE))
}
