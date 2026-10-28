from pathlib import Path
import csv,collections,re,json
import argparse,shutil
ap=argparse.ArgumentParser(description='Rebuild separated call tiers from exported annotation records.')
ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
args=ap.parse_args();D=args.output.resolve();assert not D.exists(), 'Refuse to overwrite output'
D.mkdir(parents=True)
for p in args.input.resolve().glob('*.csv'):shutil.copy2(p,D/p.name)
def read(n):return list(csv.DictReader((D/n).open()))
def write(n,rows):
    rows=list(rows)
    with (D/n).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
syn=dict(eleos='dynamins',casadaptation='cas',dgtpase='dxtpase',avast='avs',dsr1='dsr',qatabcd='qat',ptsspabcd='sspbcde',argonaute='pago',juk='jukab',pdt71or5='pdt75',tirnlr='tirnlr')
def norm(x):
    y=re.sub(r'^(gao|rst)[-_]','',x.lower());y=re.sub('[-_ ]','',y);y=re.sub('(other|fused)$','',y);y=re.sub(r'class(i+|1|2)$','',y)
    return syn.get(y,y)
allcalls=read('calls_classified.csv');og=read('outgroup_calls.csv')
atype={x['system_type']:x['activity'] for x in allcalls if x['tool']=='defensefinder'}
anti_types={x for x,y in atype.items() if y=='Antidefense'}
for x in og:
    x['activity']='Antidefense' if x['tool']=='defensefinder' and (x['system_type'] in anti_types or x['system_type'].startswith('Anti_')) else 'Defense' if x['tool']=='defensefinder' else 'Defense_candidate'
    x['host_defense_included']=x['activity']!='Antidefense'
    x['source_record_id']=x['defense_id'];x['n_genes_resolved']=x['n_genes'];x['coords_resolved']=1
# Exact duplicate rows in exported outgroup table are removed before graph construction.
og=list({tuple(sorted((k,str(v)) for k,v in x.items())):x for x in og}.values())
write('outgroup_calls_classified.csv',og)
phages=collections.defaultdict(list)
for p in read('prophage_consensus.csv'):phages[p['chrom']].append((int(p['start']),int(p['end'])))
plasmids={x['seq_name'] for x in read('plasmid_candidates.csv')}
def location(chrom,s,e):
    inside=overlap=0;dist=10**12
    for a,b in phages[chrom]:
        inside|=s>=a and e<=b;overlap|=s<=b and e>=a
        dist=min(dist,max(a-e,s-b,0))
    label='contained' if inside else 'boundary_overlap' if overlap else 'flank_5kb' if dist<=5000 else 'flank_10kb' if dist<=10000 else 'no_context'
    return dict(context_class=label,contained=int(bool(inside)),overlap=int(bool(overlap)),near_5kb=int(dist<=5000),near=int(dist<=10000),distance_bp='' if dist==10**12 else dist,predicted_plasmid=int(chrom in plasmids))
def make_tiers(rows,cohort,category):
    rows=[dict(x,cname=norm(x['system_type']),start=int(x['start']),stop=int(x['stop'])) for x in rows if (x['activity']=='Antidefense')==(category=='Antidefense')]
    assert all(x['start']>0 and x['stop']>=x['start'] for x in rows)
    bychrom=collections.defaultdict(list)
    for i,x in enumerate(rows):bychrom[x['chrom']].append(i)
    parents={t:list(range(len(rows))) for t in ['T3','T4']}
    def find(p,i):
        while p[i]!=i:p[i]=p[p[i]];i=p[i]
        return i
    def union(p,i,j):p[find(p,j)]=find(p,i)
    for ids in bychrom.values():
        ids.sort(key=lambda i:rows[i]['start'])
        for ii,i in enumerate(ids):
            for j in ids[ii+1:]:
                if rows[j]['start']>rows[i]['stop']:break
                union(parents['T4'],i,j)
                if rows[i]['tool']!=rows[j]['tool'] and rows[i]['cname']==rows[j]['cname']:union(parents['T3'],i,j)
    out=[]
    for tier in ['T2','T3','T4']:
        groups=collections.defaultdict(list)
        for i,x in enumerate(rows):groups[i if tier=='T2' else find(parents[tier],i)].append(x)
        for count,group in enumerate(groups.values(),1):
            names={x['cname'] for x in group};s=min(x['start'] for x in group);e=max(x['stop'] for x in group)
            item=dict(unit_id=f'{cohort}_{category}_{tier}_{count}',cohort=cohort,category=category,tier=tier,genome_id=group[0]['genome_id'],chrom=group[0]['chrom'],start=s,stop=e,cname=';'.join(sorted(names)),display_name=' / '.join(sorted({x['system_type'].upper() for x in group})),n_constituent_calls=len(group),tools=';'.join(sorted({x['tool'] for x in group})),member_defense_ids=';'.join(x['defense_id'] for x in group))
            item.update(location(item['chrom'],s,e) if cohort=='Desulfovibrio' else {k:'' for k in ['context_class','contained','overlap','near_5kb','near','distance_bp','predicted_plasmid']})
            out.append(item)
    return out
out=[]
for rows,cohort in [(allcalls,'Desulfovibrio'),(og,'Outgroup')]:
    for category in ['Defense','Antidefense']:out.extend(make_tiers(rows,cohort,category))
write('call_tiers.csv',out)
summary=[dict(cohort=cohort,category=category,tier=tier,n_units=len(g),n_near=sum(int(x['near']) for x in g) if cohort=='Desulfovibrio' else '') for (cohort,category,tier),g in __import__('itertools').groupby(sorted(out,key=lambda x:(x['cohort'],x['category'],x['tier'])),key=lambda x:(x['cohort'],x['category'],x['tier'])) for g in [list(g)]]
write('tier_totals.csv',summary)
print(json.dumps(summary,indent=2))
# Compare reconstructed coordinate relations against the archived row-level relation table.
old={x['defense_id']:x for x in allcalls};changed=[]
for x in out:
    if x['cohort']=='Desulfovibrio' and x['tier']=='T2':
        a=old[x['member_defense_ids']]
        if int(x['near'])!=int(a['near']=='True'):changed.append(dict(defense_id=a['defense_id'],old_near=a['near'],new_near=x['near']))
write('coordinate_relation_differences.csv',changed or [dict(defense_id='',old_near='',new_near='')])
print('Archived/recomputed T2 near relation differences:',len(changed))
