"""Prepare or verify the registered frozen inputs; no downloads or inference."""
from pathlib import Path
import argparse, collections, csv, hashlib, json, math, zipfile
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
    with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def rows(p):
    with p.open(newline='')as f:return list(csv.DictReader(f))
def verify(bundle):
    manifest=json.loads((ROOT/'metadata/input-manifest.json').read_text())['files']
    for r in manifest:
        p=bundle/r['path']
        assert p.is_file()and p.stat().st_size==r['bytes']and sha(p)==r['sha256'],r['path']
    d=bundle/'inputs/data';m=rows(d/'metadata_corrected.csv');calls=rows(d/'calls_classified.csv');t=rows(d/'call_tiers.csv')
    defense=sum(r['activity']!='Antidefense'for r in calls);anti=len(calls)-defense
    merged=sum((r['cohort'],r['category'],r['tier'])==('Desulfovibrio','Defense','T3')for r in t)
    assert(len(m),defense,anti,merged)==(436,6937,493,5217)
    assert sum(not r['gtdb_species']for r in m)==11
    assert len({r['gtdb_species']for r in m if r['gtdb_species']})==62
    assert sum(r['ecology_category']in('human_gut','animal_gut')for r in m)==397
    assert len(rows(d/'pdcs31_call_evidence.csv'))==42
    strata=rows(d/'pdcs31_gene_position_strata.csv');results=rows(d/'pdcs31_gene_position_null.csv');checks=[]
    assert len(strata)==177
    for r in results:
        group=[b for b in strata if(b['scope'],b['distance_bp'])==(r['scope'],r['distance_bp'])]
        n=sum(int(b['n'])for b in group);x=sum(int(b['x'])for b in group);expected=sum(int(b['n'])*int(b['K'])/int(b['N'])for b in group)
        assert n==int(r['n_loci'])and x==int(r['observed'])and math.isclose(expected,float(r['expected']),abs_tol=1e-12,rel_tol=1e-13)
        checks.append(dict(scope=r['scope'],distance_bp=int(r['distance_bp']),n=n,observed=x,expected=expected))
    register=rows(bundle/'supporting-secondary/outgroup_assembly_register.csv')
    assert len(register)==157 and sum(r['quality_eligible']=='TRUE'for r in register)==95
    return dict(status='pass',input_files=len(manifest),core_genomes=436,defense_calls=defense,antidefense_calls=anti,T3_defense_candidates=merged,spatial_strata=177,spatial_margin_checks=checks,no_inference=True)
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--verify',type=Path);args=ap.parse_args()
    if args.verify:
        print(json.dumps(verify(args.verify.resolve()),indent=2));return
    info=json.loads((ROOT/'metadata/archive.json').read_text());archive=ROOT/info['path']
    assert sha(archive)==info['sha256']and archive.stat().st_size==info['bytes'],'Archive fingerprint mismatch'
    dest=ROOT/'analyses/data/00-raw/d00-resource/pc057'
    if dest.exists():
        print(json.dumps(verify(dest),indent=2));return
    dest.parent.mkdir(parents=True,exist_ok=True)
    expected={r['path']for r in json.loads((ROOT/'metadata/input-manifest.json').read_text())['files']}
    with zipfile.ZipFile(archive)as z:
        assert set(z.namelist())==expected and z.testzip()is None
        for n in z.namelist():
            p=Path(n);assert not p.is_absolute()and '..'not in p.parts
        dest.mkdir()
        for n in z.namelist():z.extract(n,dest)
    # An interrupted extraction fails verification on the next invocation;
    # existing input payloads are never silently replaced or repaired.
    print(json.dumps(verify(dest),indent=2))
if __name__=='__main__':main()
