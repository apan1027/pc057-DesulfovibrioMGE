#!/usr/bin/env python3
"""Pick one TerL marker protein per viral contig and per reference phage genus.

Reads the hmmsearch domain tables produced by build_terl_tree.sh, writes one
domain-envelope FASTA per TerL Pfam profile plus a tip-metadata TSV. Reference
phages come from INPHARED and are thinned to one genome per ICTV genus so that
the reference backbone samples described diversity evenly.
"""
import argparse
import collections
import gzip
import sqlite3


def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if name:
                    yield name, ''.join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        yield name, ''.join(buf)


def best_domain(domtbl, profiles, min_score):
    """Highest-scoring envelope per (target protein, profile)."""
    best = {}
    with open(domtbl) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.split()
            prof, target, score = f[3], f[0], float(f[13])
            if prof not in profiles or score < min_score:
                continue
            key = (target, prof)
            if key not in best or score > best[key][0]:
                best[key] = (score, int(f[19]), int(f[20]))
    return best


def inphared_taxonomy(tsv_gz):
    tax = {}
    with gzip.open(tsv_gz, 'rt') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        i_acc, i_cls = header.index('Accession'), header.index('Classification')
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) <= i_cls:
                continue
            lineage = f[i_cls].split()
            family = next((x for x in lineage if x.endswith('viridae')), 'unclassified')
            genus = next((x for x in lineage if x.endswith('virus') and x != 'Caudoviricetes'),
                         'unclassified')
            tax[f[i_acc]] = (family, genus)
    return tax


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--query-faa', required=True)
    ap.add_argument('--query-domtbl', required=True)
    ap.add_argument('--ref-faa', required=True)
    ap.add_argument('--ref-domtbl', required=True)
    ap.add_argument('--inphared-tsv', required=True)
    ap.add_argument('--db', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--profiles', nargs='+', default=['Terminase_6N', 'GpA_ATPase'])
    ap.add_argument('--min-score', type=float, default=25.0)
    ap.add_argument('--pad', type=int, default=30,
                    help='residues kept on each side of the Pfam envelope')
    ap.add_argument('--min-len', type=int, default=60)
    ap.add_argument('--max-per-family', type=int, default=6,
                    help='reference genera kept per ICTV family, highest-scoring first')
    args = ap.parse_args()

    profiles = set(args.profiles)
    q_dom = best_domain(args.query_domtbl, profiles, args.min_score)
    r_dom = best_domain(args.ref_domtbl, profiles, args.min_score)

    con = sqlite3.connect(args.db)
    family_of = dict(con.execute(
        'select virus_id, family_prediction from vcontact3_assignments'))
    top5 = [r[0] for r in con.execute(
        'select family_prediction from vcontact3_assignments '
        'group by 1 order by count(*) desc limit 5')]
    con.close()
    short = {f: 'nf{}o{}'.format(f.split('_')[2], f.split('_')[6]) for f in top5}

    seqs = dict(read_fasta(args.query_faa))
    seqs.update(dict(read_fasta(args.ref_faa)))
    tax = inphared_taxonomy(args.inphared_tsv)

    rows = []
    for prof in args.profiles:
        per_contig, per_genus = {}, {}
        for (target, p), (score, start, end) in q_dom.items():
            if p != prof:
                continue
            contig = target.rsplit('_', 1)[0]
            if contig not in per_contig or score > per_contig[contig][0]:
                per_contig[contig] = (score, target, start, end)
        per_accession = {}
        for (target, p), (score, start, end) in r_dom.items():
            if p != prof:
                continue
            accession = target.rsplit('_', 1)[0]
            if accession not in per_accession or score > per_accession[accession][0]:
                per_accession[accession] = (score, target, start, end)
        for accession, hit in per_accession.items():
            genus = tax.get(accession, ('unclassified', 'unclassified'))[1]
            if genus not in per_genus or hit[0] > per_genus[genus][0]:
                per_genus[genus] = hit + (accession,)
        # Thin the backbone to a fixed number of genera per ICTV family. Every
        # family stays represented, and the families large enough to serve as the
        # positive control keep enough genera to be scored.
        by_family = collections.defaultdict(list)
        for genus, hit in per_genus.items():
            by_family[tax.get(hit[4], ('unclassified', 'unclassified'))[0]].append((genus, hit))
        per_genus = {}
        for family, entries in by_family.items():
            entries.sort(key=lambda kv: kv[1][0], reverse=True)
            for genus, hit in entries[:args.max_per_family]:
                per_genus[genus] = hit

        def envelope(target, start, end):
            seq = seqs[target]
            a = max(1, start - args.pad)
            b = min(len(seq), end + args.pad)
            return seq[a - 1:b]

        path = '{}/{}.dom.faa'.format(args.outdir, prof)
        n_query = n_ref = 0
        with open(path, 'w') as out:
            for contig, (score, target, start, end) in sorted(per_contig.items()):
                sub = envelope(target, start, end)
                if len(sub) < args.min_len:
                    continue
                n_query += 1
                tip = 'Q{:04d}'.format(n_query)
                fam = family_of.get(contig, 'unassigned')
                rows.append([prof, tip, 'query', contig, target, fam,
                             short.get(fam, 'other'), '', '{:.1f}'.format(score), len(sub)])
                out.write('>{}\n{}\n'.format(tip, sub))
            for genus, (score, target, start, end, accession) in sorted(per_genus.items()):
                sub = envelope(target, start, end)
                if len(sub) < args.min_len:
                    continue
                n_ref += 1
                tip = 'R{:04d}'.format(n_ref)
                fam, gen = tax.get(accession, ('unclassified', 'unclassified'))
                rows.append([prof, tip, 'reference', accession, target, fam,
                             'reference', gen, '{:.1f}'.format(score), len(sub)])
                out.write('>{}\n{}\n'.format(tip, sub))
        print('{}: {} query + {} reference genera -> {}'.format(prof, n_query, n_ref, path))

    meta_path = '{}/terl_taxa_metadata.tsv'.format(args.outdir)
    with open(meta_path, 'w') as fh:
        fh.write('profile\ttip_id\ttip_type\tsource_id\tprotein_id\tfamily\t'
                 'family_short\tref_genus\tdomain_score\tdomain_len\n')
        for r in rows:
            fh.write('\t'.join(map(str, r)) + '\n')
    print('metadata -> ' + meta_path)

    # per-family TerL recovery across the whole placed query set
    hit_contigs = {t.rsplit('_', 1)[0] for (t, p) in q_dom}
    rec = collections.Counter()
    tot = collections.Counter()
    for virus, fam in family_of.items():
        tot[fam] += 1
        if virus in hit_contigs:
            rec[fam] += 1
    rec_path = '{}/terl_recovery_by_family.tsv'.format(args.outdir)
    with open(rec_path, 'w') as fh:
        fh.write('family\tn_query\tn_terl\tpct_terl\n')
        for fam, n in tot.most_common():
            fh.write('{}\t{}\t{}\t{:.1f}\n'.format(fam, n, rec[fam], 100.0 * rec[fam] / n))
    print('recovery -> ' + rec_path)


if __name__ == '__main__':
    main()
