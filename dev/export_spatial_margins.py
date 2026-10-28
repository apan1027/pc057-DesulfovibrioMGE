from pathlib import Path
import csv, json, hashlib, math
from collections import defaultdict
def rows(p):
    with p.open(newline='') as f:
        return list(csv.DictReader(f))

def write_csv(p, data, fields=None):
    with p.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields or list(data[0]))
        w.writeheader(); w.writerows(data)

def digest(p):
    return hashlib.file_digest(p.open('rb'), 'sha256').hexdigest()

def export_margins(data):
    meta = rows(data/'metadata_corrected.csv')
    for r in meta:
        assert r['completeness'] and r['contamination'] and r['accession']
    reps = {}
    for r in sorted(meta, key=lambda r: (-float(r['completeness']), float(r['contamination']), r['accession'])):
        reps.setdefault(r['clone_cluster'], r['genome_id'])
    chosen = set(reps.values())
    calls = [r for r in rows(data/'call_tiers.csv') if (r['cohort'],r['category'],r['tier'],r['cname']) == ('Desulfovibrio','Defense','T3','pdcs31')]
    features, prophages = defaultdict(list), defaultdict(list)
    for r in rows(data/'pdcs31_contig_features.csv'):
        if r['feature_type'].lower()=='cds': features[r['chrom']].append((int(r['start']),int(r['stop'])))
    for r in rows(data/'prophage_consensus.csv'):
        prophages[r['chrom']].append((int(r['start']),int(r['end'])))
    def context(chrom, s, e, d):
        return any(s<=pe+d and e>=ps-d for ps,pe in prophages[chrom])
    result, checks = [], []
    saved = {(r['scope'],int(r['distance_bp'])):r for r in rows(data/'pdcs31_gene_position_null.csv')}
    for scope in ['all_loci','clone_representatives']:
        groups = defaultdict(list)
        for r in calls:
            if scope=='all_loci' or r['genome_id'] in chosen: groups[r['chrom']].append(r)
        for distance in [0,5000,10000]:
            block=[]
            for chrom, cc in groups.items():
                lengths=[int(r['stop'])-int(r['start'])+1 for r in cc]
                eligible=[(s,e) for s,e in features[chrom] if min(lengths)*.5<=e-s+1<=max(lengths)*2]
                block.append(dict(scope=scope,distance_bp=distance,chrom=chrom,N=len(eligible),K=sum(context(chrom,s,e,distance) for s,e in eligible),n=len(cc),x=sum(context(chrom,int(r['start']),int(r['stop']),distance) for r in cc)))
            reference=saved[scope,distance]
            n=sum(r['n'] for r in block); x=sum(r['x'] for r in block)
            expected=sum(r['n']*r['K']/r['N'] for r in block)
            assert n==int(reference['n_loci']) and x==int(reference['observed'])
            assert math.isclose(expected,float(reference['expected']),rel_tol=1e-13,abs_tol=1e-12)
            checks.append(dict(scope=scope,distance_bp=distance,n_strata=len(block),n_loci=n,observed=x,expected_from_margins=expected,expected_saved=float(reference['expected']),match=True))
            result.extend(block)
    old=rows(data.parent/'previous-export/pdcs31_gene_position_strata.csv')
    new=[r for r in result if r['scope']=='all_loci' and r['distance_bp']==10000]
    normalize=lambda rr: sorted((r['chrom'],*[int(r[k]) for k in ['N','K','n','x']]) for r in rr)
    assert normalize(old)==normalize(new)
    return result, checks


if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(description='Export existing spatial margins only; no inference.')
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--bundle',type=Path,required=True)
    args=ap.parse_args();assert not args.output.exists(), 'Refuse to overwrite an existing output'
    data=args.bundle.resolve()/'inputs/data'
    margins,checks=export_margins(data);write_csv(args.output,margins)
    print(json.dumps(checks,indent=2))
