from pathlib import Path
import sqlite3,csv,json,re,hashlib,shutil,collections,datetime
import argparse
ap=argparse.ArgumentParser(description='Prepare revised records from the original annotation database, read-only.')
ap.add_argument('--database',type=Path,required=True);ap.add_argument('--bundle',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
args=ap.parse_args();D=args.output.resolve();S=args.bundle.resolve()
assert args.database.is_file();assert not D.exists(), 'Refuse to overwrite output';D.mkdir(parents=True)
c=sqlite3.connect('file:'+str(args.database.resolve())+'?mode=ro',uri=True);c.row_factory=sqlite3.Row
def write(name,rows):
    rows=list(rows)
    if not rows:return
    with (D/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def table(name):return [dict(x) for x in c.execute('select * from '+name)]
raw={r['accession']:r for r in map(json.loads,(S/'inputs/provenance/ncbi_summary.jsonl').open())}
meta=[x for x in table('genome_metadata') if x['analysis_set']=='core']
prov={x['assembly_accession']:x for x in csv.DictReader((S/'inputs/provenance/table_1907_genome_provenance.csv').open())}
gut=re.compile(r'\bgut\b|gastrointestinal|gastro.intestinal|\bintestin(?:e|es|al)\b|caec|cecum|\bcolon\b|jejun|ileum|feces|faeces|fecal|faecal|stool|rumen|digesta|\bmanure\b|\bdung\b|digestive tract|digestive system',re.I)
human=re.compile(r'\bhuman\b|homo sapiens',re.I)
animal=re.compile(r'\bpigs?\b|piglet|\bswine\b|sus scrofa|gallus|chicken|mouse|murine|mus musculus|\brat\b|rattus|marmot|marmota|rumen|cattle|bovine|bos taurus|sheep|ovine|ovis|goat|capra|deer|rangifer|camel|buffalo|vole|myodes|gorilla|termite|trichonympha|hodotermopsis|chinchilla|\bcat\b|felis|\bdog\b|canis|equine|horse|rabbit|kangaroo|koala|wombat',re.I)
audit=[]
for m in meta:
    old=m['ecology_category']; j=raw.get(m['accession'],{});a=j.get('assembly_info',{});bs=a.get('biosample',{})
    attrs=bs.get('attributes',[])
    kv={x['name']:x.get('value','') for x in attrs if re.search(r'source|site|body|host|env_|tissue|gut|habitat|isolation',x.get('name',''),re.I)}
    local=' | '.join(f'{k}: {v}' for k,v in kv.items())
    primary=' | '.join(str(m.get(k) or '') for k in ['isolation_source','host','bs_title'])+' | '+local
    txt=primary.lower().replace('_',' ')
    iso=(m['isolation_source'] or '').lower()
    # Explicit specimen/environment evidence takes precedence over host names.
    if re.search(r'mouth|\boral\b|subgingival',iso+' '+str(kv.get('env_medium','')).lower()):
        new='host_non_gut';reason='explicit oral sampling site'
    elif re.search(r'ab[cs]cess|abcess|\bblood\b|wound',iso):
        new='clinical_non_gut';reason='explicit extra-intestinal clinical site'
    elif re.search(r'ant-built patch|fungus garden',txt):
        new='ant_associated_material';reason='ant-built material, not an animal intestinal sample'
    elif re.search(r'epidermal|skin surface|subgingival',txt):
        new='host_non_gut';reason='explicit non-gut body site'
    elif re.search(r'in vitro culture derivated|grow together with monocercomonoides',txt):
        new='gut_derived_culture';reason='mixed in-vitro community derived from gut; excluded from direct gut contrast'
    elif re.search(r'wastewater|sewage|sludge|bioreactor|biodigester|anaerobic digester|landfill|leachate|effluent',txt):
        new='wastewater';reason='explicit engineered waste/reactor/landfill source'
    elif re.search(r'\bsoil\b|\bclay\b|rhizosphere',txt):
        new='soil';reason='explicit terrestrial substrate'
    elif re.search(r'sediment',txt):
        new='sediment';reason='explicit sediment source'
    elif re.search(r'freshwater|\briver\b|\blake\b|groundwater|\bpond\b',iso):
        new='freshwater';reason='explicit freshwater source'
    elif gut.search(txt):
        if human.search(txt):new='human_gut';reason='gut-specific field and human provenance'
        elif animal.search(txt):new='animal_gut';reason='gut-specific field and nonhuman provenance'
        else:new='gut_host_unknown';reason='gut site supported, human/nonhuman host unresolved'
    elif re.match(r'^MGBC\d+',m['bs_title'] or ''):
        new='animal_gut';reason='Mouse Gastrointestinal Bacterial Catalogue specimen identifier'
    elif human.search(txt) or animal.search(str(m.get('host') or '')):
        new='host_site_unknown';reason='host identity alone does not establish intestinal sampling site'
    else:new='unknown';reason='insufficient site/host evidence'
    m['ecology_original']=old;m['ecology_category']=new;m['ecology_evidence']=reason
    m['gtdb_species_original']=m['gtdb_species'];m['gtdb_species']=(m['gtdb_species'] or '').strip() or None
    m.update({k:prov.get(m['accession'],{}).get(k,'') for k in ['genome_id','bioproject','submitter','clone_cluster']})
    m['species_assigned']=bool(m['gtdb_species']);m['gut_comparison_eligible']=bool(new in ['human_gut','animal_gut'] and m['gtdb_species'])
    audit.append(dict(accession=m['accession'],genome_id=m['genome_id'],original=old,revised=new,changed=old!=new,reason=reason,isolation_source=m['isolation_source'],host=m['host'],biosample_title=m['bs_title'],biosample_evidence=local,biosample_record_available=bool(j),species_original=m['gtdb_species_original'],species_assigned=bool(m['gtdb_species'])))
write('metadata_corrected.csv',meta);write('metadata_audit.csv',audit)
activity={x['sys_uid']:x['activity'] for x in table('defense_finder_systems')}
calls=table('defense_catalog'); unique={x['defense_id']:x for x in calls}
assert len(calls)==7440 and len(unique)==7430
canon={x['defense_id']:x for x in csv.DictReader((S/'inputs/provenance/table_1705_defense_canonical_map.csv').open())}
ov=table('prophage_defense_overlap');byov=collections.defaultdict(list)
for x in ov:byov[x['defense_id']].append(x)
plasmids={x['seq_name'] if 'seq_name' in x else x.get('chrom') for x in table('plasmid_candidates')}
for x in unique.values():
    x['activity']=activity[x['source_record_id']] if x['tool']=='defensefinder' else 'Defense_candidate'
    x['host_defense_included']=x['activity']!='Antidefense'
    x.update({k:canon[x['defense_id']][k] for k in ['cname','sys_canonical','sys_locus']})
    rel=byov[x['defense_id']]
    x['near']=bool(rel);x['overlap']=any(z['overlap_bp']>0 for z in rel)
    x['contained']=any(x['start']>=z['prophage_start'] and x['stop']<=z['prophage_end'] for z in rel)
    x['context_class']='contained' if x['contained'] else 'boundary_overlap' if x['overlap'] else 'flank_5kb' if any(z['relation']=='flank_5kb' for z in rel) else 'flank_10kb' if rel else 'no_context'
    x['predicted_plasmid']=x['chrom'] in plasmids
write('calls_classified.csv',unique.values());write('calls_raw.csv',calls);write('overlaps.csv',ov)
for t in ['genome_mge_summary','prophage_consensus','contig_map','plasmid_candidates','arg_consensus','vf_consensus','viral_catalog','checkv_quality','vcontact3_assignments','padloc_hits']:
    write(t+'.csv',table(t))
# Export only CDS and genes around PDC-S31-bearing contigs, sufficient for a same-contig gene-position null.
pdcs31=[x for x in unique.values() if x['system_type'].upper()=='PDC-S31'];ctgs={x['chrom'] for x in pdcs31}
features=[dict(x) for x in c.execute('select * from bakta_features where chrom in ('+','.join('?'*len(ctgs))+')',list(ctgs))]
write('pdcs31_contig_features.csv',features)
write('outgroup_calls.csv',table('outgroup_defense_catalog'));write('outgroup_per_genome.csv',table('outgroup_comparison_per_genome'))

print('Prepared revised metadata and call records; no model fitting or database write performed.')
