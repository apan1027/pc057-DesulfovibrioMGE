from pathlib import Path
import csv,re,json,collections
import argparse
ap=argparse.ArgumentParser(description='Export current supplementary CSV tables from frozen records; no inference.')
ap.add_argument('--bundle',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
args=ap.parse_args();D=args.bundle.resolve()/'inputs/data';OUT=args.output.resolve()
assert not OUT.exists(), 'Output directory already exists';OUT.mkdir(parents=True)
rd=lambda n:list(csv.DictReader((D/(n+'.csv')).open()))
meta=rd('metadata_corrected');md={x['genome_id']:x for x in meta};calls=rd('calls_classified');tiers=rd('call_tiers')
pg={x['genome_id']:dict(genome_id=x['genome_id'],accession=x['accession'],gtdb_species=x['gtdb_species'],ecology_category=x['ecology_category'],species_assigned=x['species_assigned'],bioproject=x['bioproject'],clone_cluster=x['clone_cluster']) for x in meta}
for x in rd('per_genome_corrected'):
    pg[x['genome_id']][x['category']+'_'+x['tier']]=x['n']
    for k in ['n_consensus_prophage','n_viral_total','n_plasmid_candidate','genome_size_bp','completeness','contamination','num_contigs']:
        pg[x['genome_id']][k]=x[k]
pdc=rd('pdcs31_call_evidence');assert len(pdc)==42
quality={x['contig_id']:x for x in rd('checkv_quality')};vc=[]
for x in rd('viral_catalog'):
    q=quality.get(x['virus_id'],{})
    vc.append(dict(**x,checkv_quality=q.get('checkv_quality','Not-assessed'),checkv_completeness=q.get('completeness',''),checkv_contamination=q.get('contamination','')))
diag={(x['category'],x['tier'],x['model_key']):x for x in rd('model_diagnostics')}
ecology=rd('ecology_models_corrected')
for x in ecology:
    key={'species_project_random':'study','poisson_species':'poisson','NB_species':'NB'}.get(x['model'],x['model'])
    a=diag[(x['category'],x['tier'],key)]
    x.update({k:a[k] for k in ['optimizer_message','nb_theta','poisson_boundary']})
specs=[
 ('S1_Metadata',meta,'metadata_corrected.csv'),
 ('S2_Load',list(pg.values()),'per_genome_corrected.csv, pivoted to one row per genome'),
 ('S3_DefenseCalls',[x for x in calls if x['activity']!='Antidefense'],'calls_classified.csv, Defense + Defense_candidate'),
 ('S4_AntiDefense',[x for x in calls if x['activity']=='Antidefense'],'calls_classified.csv, Antidefense only'),
 ('S5_MergedDefense',[x for x in tiers if x['cohort']=='Desulfovibrio' and x['category']=='Defense' and x['tier']=='T3'],'call_tiers.csv, T3 defense'),
 ('S6_Enrichment',rd('enrichment_corrected'),'enrichment_corrected.csv'),
 ('S7_ClusterIntervals',rd('cluster_intervals_corrected'),'cluster_intervals_corrected.csv'),
 ('S8_PDC_S31',pdc,'pdcs31_call_evidence.csv; raw PADLOC matches and specimen records'),
 ('S9_SpatialNull',rd('pdcs31_gene_position_null'),'pdcs31_gene_position_null.csv'),
 ('S10_EcologyModels',ecology,'ecology_models_corrected.csv + model_diagnostics.csv; Wald intervals, nominal P; NB boundary flagged'),
 ('S11_Outgroup',rd('outgroup_corrected'),'outgroup_corrected.csv; 593 genomes x3 conventions; use qc_pass and ecology scope'),
 ('S12_OutgroupTests',rd('outgroup_tests_corrected'),'outgroup_tests_corrected.csv'),
 ('S13_OutgroupModels',rd('outgroup_models_corrected'),'outgroup_models_corrected.csv; genome-level exploratory adjustment'),
 ('S14_ARG',rd('arg_consensus'),'archived arg_consensus, unchanged'),
 ('S15_VF_High',[x for x in rd('vf_consensus') if x['tier'].upper()=='HIGH'],'archived vf_consensus, HIGH tier, unchanged'),
 ('S16_ViralCatalog',vc,'archived viral_catalog + checkv_quality, unchanged'),
 ('S17_ViralLabels',rd('vcontact3_assignments'),'archived vcontact3_assignments, unchanged'),
 ('S18_MetadataAudit',rd('metadata_audit'),'metadata_audit.csv; original/revised labels and source evidence'),
 ('S19_SpatialStrata',rd('pdcs31_gene_position_strata'),'pdcs31_gene_position_strata.csv; full/clone sets at 0, 5 and 10 kb; 177 strata'),
 ('S20_SharedSpecies',rd('shared_species_corrected'),'shared_species_corrected.csv; assigned species only'),
 ('S21_NeighborhoodGenes',rd('neighborhood_figure_genes'),'neighborhood_figure_genes.csv; Figure4 representative windows')]
exported=[]
for name,rows,source in specs:
    target=OUT/(name+'.csv');fields=list(rows[0])
    with target.open('w',newline='')as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
    exported.append(dict(table=name,rows=len(rows),columns=len(fields),source=source))
with(OUT/'table-export-manifest.json').open('w')as f:json.dump(exported,f,indent=2)
assert len(exported)==21
print(json.dumps(exported,indent=2))
