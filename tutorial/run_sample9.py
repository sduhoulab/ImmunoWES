#!/usr/bin/env python3
"""Single-patient adapter for the supplied gastric_cancer code snapshot."""
import csv, gzip, hashlib, json, shlex, shutil
import subprocess, sys, tarfile, traceback
from datetime import datetime
from pathlib import Path

cfgfile = Path(sys.argv[1] if len(sys.argv) > 1 else 'sample9.json').resolve()
C = json.loads(cfgfile.read_text())
W = (cfgfile.parent / C['work']).resolve()
S = W / 'source/gastric_cancer'
R = W / 'results'
for d in [W / 'logs', W / 'evidence', W / 'adapted', R]:
    d.mkdir(parents=True, exist_ok=True)
mode = sys.argv[2] if len(sys.argv) > 2 else 'all'
logpath = W / 'logs' / (datetime.now().strftime('%Y%m%d_%H%M%S_%f') + '_' + mode + '.log')
LOG = logpath.open('w', buffering=1)
GENES = ['DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
# ENA filereports, retrieved 2026-10-09. Do not infer mapping from directory order.
RUNS = [('YSU_N9', 'SRR1501968', 21138874,
         ['f35070d88c622e9645d15e3265e3d81f', '9695788ef4c3ed3cf31796c4781e97e7']),
        ('YSU_T9', 'SRR1501970', 12243645,
         ['733ba51269fe1477e40c5c340e12e793', 'b454b4f81ab09ad91660cbe7a836f28c'])]

def say(x):
    print(x, flush=True); print(x, file=LOG, flush=True)

def need(paths=(), tools=()):
    bad = [str(p) for p in paths if not Path(p).is_file()]
    bad += [t for t in tools if not shutil.which(t)]
    if bad:
        raise RuntimeError('BLOCKED: missing/unreadable: ' + ', '.join(bad))

def run(args, cwd=None):
    say('$ ' + shlex.join(map(str, args)))
    subprocess.run(list(map(str, args)), cwd=cwd, stdout=LOG, stderr=subprocess.STDOUT, check=True)

def digest(p, algorithm='sha256'):
    with Path(p).open('rb') as f:
        return hashlib.file_digest(f, algorithm).hexdigest()

def table(p, rows):
    with Path(p).open('w', newline='') as f:
        csv.writer(f, delimiter='\t').writerows(rows)

def prepare():
    need([C['archive']])
    # Extract only regular code files and final sample9 HLA reports; never run archive scripts.
    with tarfile.open(C['archive']) as t:
        for m in t:
            p = Path(m.name)
            selected = m.name.endswith(('.sh', '.py', '.pl', '.R')) or (
                any('/' + x + '/' in m.name for x in ['YSU_N9', 'YSU_T9']) and
                (m.name.endswith('_final.result.txt') or p.name in ['R1_bestguess.txt', 'R1_bestguess_G.txt']))
            if not m.isfile() or not selected:
                continue
            if p.is_absolute() or '..' in p.parts:
                raise RuntimeError('Unsafe archive path')
            dest = W / 'source' / p
            data = t.extractfile(m).read()
            if dest.exists() and dest.read_bytes() != data:
                raise RuntimeError('Source snapshot differs: ' + str(dest))
            dest.parent.mkdir(parents=True, exist_ok=True)
            if not dest.exists(): dest.write_bytes(data)
    say('Prepared original source snapshot (no source edits).')

def audit():
    rows = [['sample', 'file', 'MD5', 'expected_MD5', 'records', 'status']]
    bad = []
    for sample, accession, expected_n, hashes in RUNS:
        for mate, expected in enumerate(hashes, 1):
            p = Path(C['raw']) / f'{accession}_{mate}.fastq.gz'
            n, md5, status = 0, '', 'OK'
            if not p.is_file():
                status = 'MISSING'
            else:
                before = p.stat()
                md5 = digest(p, 'md5')
                try:
                    with gzip.open(p, 'rb') as f:
                        while h := f.readline():
                            seq, plus, qual = f.readline().rstrip(), f.readline(), f.readline().rstrip()
                            if not h.startswith(b'@') or not plus.startswith(b'+') or len(seq) != len(qual):
                                raise ValueError(f'bad FASTQ record {n+1}')
                            if h.split()[0] != f'@{accession}.{n+1}'.encode():
                                raise ValueError(f'unexpected read ID at {n+1}')
                            n += 1
                except (OSError, EOFError, ValueError) as e:
                    status = str(e)
                if md5 != expected: status += '; MD5 mismatch'
                if n != expected_n: status += '; record count mismatch'
                after = p.stat()
                if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                    status += '; input changed during audit'
            rows.append([sample, p.name, md5, expected, n, status])
            say(str(rows[-1]))
            if status != 'OK': bad.append(p.name)
    table(W / 'evidence/input_audit.tsv', rows)
    if bad: raise RuntimeError('BLOCKED: invalid/missing FASTQ: ' + ', '.join(bad))

def adapt(source, replacements):
    text = Path(source).read_text()
    for old, new in replacements.items(): text = text.replace(old, str(new))
    if Path(source).name == '2_0.extract_p.py':
        # Operational repair only: propagate the original caught processing error.
        text = text.replace('except Exception as e:', 'except Exception as e:\n            raise')
    dest = W / 'adapted' / Path(source).name
    dest.write_text(text)
    return dest

def hla():
    d = R / 'hla'; d.mkdir(exist_ok=True)
    reports = []
    for sample, *_ in RUNS:
        p = S / f'hla_estimation/{sample}/result/{sample}_final.result.txt'
        need([p])
        out = d / sample / 'result'; out.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, out / p.name)
        with p.open() as f:
            reports.append({r[0]: ';'.join(r[1:]) for r in csv.reader(f, delimiter='\t')})
    table(R/'hla_comparison.tsv', [['gene','normal_archived_all_columns','tumor_archived_all_columns']] +
          [[g,reports[0].get(g,''),reports[1].get(g,'')] for g in ['A','B','C']+GENES])
    # Exact original extraction/splitting logic; only replace the hard-coded root.
    for name in ['extract_hla.py', 'split_hla.py']:
        p = adapt(S / 'hla_estimation' / name,
                  {'/data/yuan/gastric_cancer/hla_estimation': d})
        run([sys.executable, p])
    rows = list(csv.reader((d / 'all_normal_HLA.txt').open(), delimiter='\t'))
    if len(rows) != 2 or rows[1][0] != 'YSU_N9' or len(rows[1]) != 11:
        raise RuntimeError('Unexpected HLA extraction output')
    say('Reused archived HLA-HD reports; not a new HLA typing run. ' + str(rows[1]))

def hed():
    hla()
    script = S / 'HED_UKB/UKB-HED/CalculateIndividualDivergence.pl'
    status = [['locus', 'status']]
    for g, prefix in zip(GENES, ['DRB', 'DQA', 'DQB', 'DPA', 'DPB']):
        p = R / 'hla' / (g + '.txt')
        row = list(csv.reader(p.open(), delimiter='\t'))[1]
        if len(row) != 3 or not all(row[1:]):
            status.append([g, 'NA: incomplete genotype; no zero imputation'])
            continue
        fasta = Path(C['hed_fasta_dir']) / (prefix + '.fas')
        need([script, fasta, C['grantham']], ['perl'])
        # Prevent the legacy script from silently producing NA for missing allele sequences.
        ids = {x[1:] for x in fasta.read_text().splitlines() if x.startswith('>')}
        if not set(row[1:]) <= ids: raise RuntimeError('HED allele missing from alignment')
        run(['perl', script, '-d', C['grantham'], '-f', fasta, '-i', p])
        out = Path(str(p) + '.IndividualDivergence.txt')
        need([out])
        values = list(csv.DictReader(out.open(), delimiter='\t'))
        if len(values) != 1: raise RuntimeError('Invalid HED output')
        float(values[0]['Divergence_Average'])
        status.append([g, 'computed'])
    table(R / 'hed_status.tsv', status)

def preflight():
    names = ['fastqc','multiqc','trim_galore','cutadapt','bwa','samtools','gatk',
             'java','hlahd.sh','bowtie2','table_annovar.pl','Rscript','MixMHC2pred_unix']
    rows = [['tool', 'resolved_path']] + [[n, shutil.which(n) or 'MISSING on PATH'] for n in names]
    table(W / 'evidence/tools.tsv', rows)
    for row in rows[1:]: say('\t'.join(row))
    for k in ['ref','dbsnp','indels','bed','grantham','proteome','humandb','mixmhc']:
        say(k + '\t' + C[k] + '\t' + ('exists' if Path(C[k]).exists() else 'MISSING/unreadable'))

def step(name, command, inputs, outputs):
    """Resume only successful steps with matching code/config, input stats and output hashes."""
    need(inputs)
    versions = W/'evidence/wes_versions.txt'
    signature = [digest(__file__), C, command, versions.read_text() if versions.exists() else '',
                 [(str(p), Path(p).stat().st_size, Path(p).stat().st_mtime_ns) for p in inputs]]
    token = hashlib.sha256(json.dumps(signature, sort_keys=True).encode()).hexdigest()
    stamp = R / (name + '.done.json')
    if stamp.exists():
        saved = json.loads(stamp.read_text())
        if saved['token'] == token and all(Path(p).is_file() and digest(p) == saved['outputs'].get(str(p)) for p in outputs):
            say('REUSE ' + name); return
        raise RuntimeError('Changed inputs/config/output: use a new work directory: ' + name)
    run(['bash', '-euo', 'pipefail', '-c', command])
    need(outputs)
    for p in outputs:
        if str(p).endswith('.bam'): run(['samtools', 'quickcheck', '-v', p])
        if str(p).endswith('.vcf') and '#CHROM\t' not in Path(p).read_text():
            raise RuntimeError('Invalid VCF header: ' + str(p))
    stamp.write_text(json.dumps({'token': token, 'outputs': {str(p): digest(p) for p in outputs}}, indent=2))

def wes():
    audit()  # Mandatory before WES, even when resuming.
    need(tools=['fastqc','multiqc','trim_galore','cutadapt','bwa','samtools','gatk'])
    ref, snp, indel, bed = [C[k] for k in ['ref','dbsnp','indels','bed']]
    need([ref, ref+'.fai', str(Path(ref).with_suffix('.dict')), snp, snp+'.tbi',
          indel, indel+'.tbi', bed] + [C['bwa_index']+x for x in ['.amb','.ann','.bwt','.pac','.sa']])
    with (W/'evidence/wes_versions.txt').open('w') as f:
        for args in [['fastqc','--version'],['multiqc','--version'],['trim_galore','--version'],
                     ['cutadapt','--version'],['bwa'],['samtools','--version'],['gatk','--version'],['java','-version']]:
            p = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            f.write(shlex.join(args)+'\n'+p.stdout+'\n')
    q = lambda x: shlex.quote(str(x))
    t = int(C['threads'])
    for d in ['input','raw_qc','clean','clean_qc','bam','vcf']: (R/d).mkdir(exist_ok=True)
    for sample, accession, *_ in RUNS:
        a,b = [Path(C['raw'])/f'{accession}_{m}.fastq.gz' for m in [1,2]]
        x,y = [R/'clean'/f'{sample}_{m}_val_{m}.fq.gz' for m in [1,2]]
        step(sample+'_rawqc', f'fastqc -t {t} -o {q(R/"raw_qc")} {q(a)} {q(b)}', [a,b],
             [R/'raw_qc'/f'{accession}_{m}_fastqc.zip' for m in [1,2]])
        links = [R/'input'/f'{sample}_{m}.fastq.gz' for m in [1,2]]
        for link, original in zip(links, [a,b]):
            if link.is_symlink() or link.exists():
                if link.resolve() != original.resolve(): raise RuntimeError('Input alias conflict')
            else: link.symlink_to(original.resolve())
        step(sample+'_trim', f'trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 3 --paired -o {q(R/"clean")} {q(links[0])} {q(links[1])}', [a,b], [x,y])
        step(sample+'_cleanqc', f'fastqc -t {t} -o {q(R/"clean_qc")} {q(x)} {q(y)}', [x,y],
             [R/'clean_qc'/f'{sample}_{m}_val_{m}_fastqc.zip' for m in [1,2]])
        bam, marked, recal, bqsr = [R/'bam'/(sample+s) for s in ['.bam','_marked.bam','_recal.table','_bqsr.bam']]
        rg = '@RG\\tID:'+sample+'\\tSM:'+sample+'\\tLB:WXS\\tPL:Illumina'
        step(sample+'_align', f'bwa mem -M -t {t} -R {q(rg)} {q(C["bwa_index"])} {q(x)} {q(y)} | samtools sort -@ {t} -m 1G -o {q(bam)} -', [x,y,ref]+[C['bwa_index']+z for z in ['.amb','.ann','.bwt','.pac','.sa']], [bam])
        metrics = R/'bam'/(sample+'.metrics')
        step(sample+'_dedup', f'gatk --java-options "-Xmx20G" MarkDuplicates -I {q(bam)} --REMOVE_DUPLICATES true -O {q(marked)} -M {q(metrics)} && samtools index -@ {t} {q(marked)}', [bam], [marked,Path(str(marked)+'.bai'),metrics])
        step(sample+'_bqsr', f'gatk --java-options "-Xmx50G" BaseRecalibrator -R {q(ref)} -I {q(marked)} --known-sites {q(snp)} --known-sites {q(indel)} -O {q(recal)} && gatk --java-options "-Xmx50G" ApplyBQSR -R {q(ref)} -I {q(marked)} -bqsr {q(recal)} -O {q(bqsr)} && samtools index {q(bqsr)}', [marked,ref,snp,indel], [recal,bqsr,Path(str(bqsr)+'.bai')])
    for d in ['raw_qc','clean_qc']:
        run(['multiqc', R/d, '-o', R/d/'multiqc', '--force'])
    tumor, normal = [R/'bam'/f'YSU_{s}9_bqsr.bam' for s in ['T','N']]
    v, f, passed = [R/'vcf'/('YSU_T9_'+s+'.vcf') for s in ['mutect2','somatic','filter']]
    step('mutect2', f'gatk --java-options "-Xmx20G" Mutect2 -R {q(ref)} -I {q(tumor)} -tumor YSU_T9 -I {q(normal)} -normal YSU_N9 -L {q(bed)} -O {q(v)}', [tumor,normal,ref,bed], [v,Path(str(v)+'.stats')])
    step('filtermutect', f'gatk FilterMutectCalls -R {q(ref)} -V {q(v)} -O {q(f)}', [v,Path(str(v)+'.stats'),ref], [f])
    # Exact filtering expression from 9.vcf.sh.
    expr = '{if(/^#/){print}else{next unless $F[6] eq "PASS";next if $F[0] =~/_/;print } }'
    step('pass', f'perl -alne {q(expr)} {q(f)} > {q(passed)}', [f], [passed])

def annotate():
    v = R/'vcf/YSU_T9_filter.vcf'
    need([v], ['table_annovar.pl'])
    protocols = ['refGene','avsnp144','gnomad41_genome','dbnsfp42c']
    need([Path(C['humandb'])/('hg38_'+p+'.txt') for p in protocols] + [Path(C['humandb'])/'hg38_refGeneMrna.fa'])
    out = R/'annotation'; out.mkdir(exist_ok=True)
    # Documented interface repair: -vcfinput receives VCF, never .avinput.
    run(['table_annovar.pl', v, C['humandb'], '-buildver','hg38','-out',out/'YSU_T9',
         '-remove','-protocol',','.join(protocols),'-operation','g,f,f,f','-nastring','.','-vcfinput'])
    need([out/'YSU_T9.hg38_multianno.vcf',out/'YSU_T9.hg38_multianno.txt'])

def peptides():
    need([R/'annotation/YSU_T9.hg38_multianno.vcf', C['proteome']])
    run([sys.executable, '-c', 'import pandas; import Bio'])
    base = R/'downstream'; base.mkdir(exist_ok=True)
    shutil.copy2(R/'annotation/YSU_T9.hg38_multianno.vcf', base)
    replacements = {'/data/yuan/gastric_cancer/downstream_ukb/gastric_vcf': base,
                    '/data/yuan/IMM/Ref-seq/UP000005640_9606.fasta': C['proteome'],
                    '/data/LiuRT/mhc2/MixMHCIIpred/gastric_cancer/15mer_mutation': R/'15mer'}
    for name in ['1.filter_vcf.sh','2.extract_vcf.sh','2_0.extract_p.py','3.match-seq.py','4.map_mutation.py','5.map-15mer.py']:
        p = adapt(S/'downstream_ukb/code'/name, replacements)
        run(['bash','-euo','pipefail',p] if name.endswith('.sh') else [sys.executable,p])
    need([R/'15mer/YSU_T9_filtered_filtered.txt'])
    say('Legacy peptide preparation completed; inspect mapping losses in log.')

def predict():
    # No automatic conversion of ambiguous/incomplete HLA to predictor alleles.
    if not C['reviewed_hla']: raise RuntimeError('BLOCKED: reviewed predictor-format HLA file required; archived normal HLA is incomplete/discordant')
    p = R/'15mer/YSU_T9_filtered_filtered.txt'
    binary = Path(C['mixmhc'])/'MixMHC2pred_unix'
    need([p, C['reviewed_hla'], binary])
    alleles = Path(C['reviewed_hla']).read_text().split()
    if not alleles or not p.read_text().strip(): raise RuntimeError('Empty predictor input')
    need([Path(C['mixmhc'])/'PWMdef'/(a+'.txt') for a in alleles])
    out = R/'YSU_T9_prediction.txt'
    run([binary,'-i',p,'-o',out,'-a',*alleles,'--no_context'], cwd=C['mixmhc'])
    rows = list(csv.DictReader((s for s in out.open() if not s.startswith('#')), delimiter='\t'))
    count = sum(float(r['%Rank_best']) < 10 for r in rows)
    table(R/'rank_best_counts.tsv', [['sample','Count<10% Rank'],['YSU_T9',count]])

actions = dict(prepare=prepare,audit=audit,preflight=preflight,hla=hla,hed=hed,
               wes=wes,annotate=annotate,peptides=peptides,predict=predict)
try:
    say('Config: '+str(cfgfile)+'; mode: '+mode+'; log: '+str(logpath))
    if mode == 'all':
        for action in [prepare,preflight,wes,annotate,hla,hed,peptides,predict]: action()
    else: actions[mode]()
    say('SUCCESS: '+mode)
except Exception as e:
    say(str(e)); traceback.print_exc(file=LOG); sys.exit(1)
finally:
    LOG.close()
