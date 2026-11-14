"""Reporting-only correctness: preserve published values, missingness and private logs."""
import gzip
import json
from pathlib import Path

from pfvarx.reporting.build import _distance_matrix, _log_content, _num, _optional_json, _redact, _vcf_values


def test_missing_values_are_not_zero():
    assert [_num(v) for v in (None, '', 'NA', '-', 'nan', 'inf', 'not available')] == [None] * 7
    assert _num('0') == 0
    assert _num('-0.25') == -0.25


def test_vcf_plain_gzip_and_missing_fields(tmp_path: Path):
    text = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n' \
           'chr1\t2\t.\tA\tG\t37\tPASS\tANN=G|missense_variant|MODERATE\tGT:DP\t0/1:9\n' \
           'chr1\t3\t.\tT\tC\t.\tPASS\t.\tGT:DP\t./.:0\n'
    for name in ('test.vcf', 'test.vcf.gz'):
        path = tmp_path / name
        if name.endswith('gz'):
            with gzip.open(path, 'wt') as handle:
                handle.write(text)
        else:
            path.write_text(text)
        warnings = []
        values = _vcf_values(path, warnings)
        assert not warnings
        assert values[('chr1', 2, 'A', 'G')]['quality'] == 37
        assert values[('chr1', 2, 'A', 'G')]['genotype'] == '0/1'
        assert values[('chr1', 2, 'A', 'G')]['consequence'] == 'missense_variant'
        assert values[('chr1', 3, 'T', 'C')]['quality'] is None
        assert values[('chr1', 3, 'T', 'C')]['genotype'] is None


def test_log_redaction_preserves_provenance_and_rejects_external_symlinks(tmp_path: Path):
    run = tmp_path / 'run'
    logs = run / 'logs'
    logs.mkdir(parents=True)
    original_hash = 'a' * 64
    payload = {'token': 'SECRET_VALUE', 'hash': original_hash,
               'text': 'api_key=PRIVATE123 Bearer BEARER123 https://user:PASS@example.org\n'
                       '-----BEGIN PRIVATE KEY-----\nPRIVATE_BODY\n-----END PRIVATE KEY-----'}
    redacted = _redact(payload)
    assert redacted['hash'] == original_hash
    encoded = json.dumps(redacted)
    for secret in ('SECRET_VALUE', 'PRIVATE123', 'BEARER123', 'user:PASS', 'PRIVATE_BODY'):
        assert secret not in encoded
    (logs / 'test.log').write_text(payload['text'])
    external = tmp_path / 'private.log'
    external.write_text('EXTERNAL_SECRET')
    (logs / 'external.log').symlink_to(external)
    entries = _log_content(run)
    assert len(entries) == 1
    assert 'PRIVATE123' not in entries[0]['text']
    assert 'EXTERNAL_SECRET' not in json.dumps(entries)


def test_optional_corruption_and_measured_distances(tmp_path: Path):
    path = tmp_path / 'optional.json'
    path.write_text('{bad json')
    warnings = []
    assert _optional_json(path, warnings) is None
    assert warnings
    matrix = tmp_path / 'tree.mldist'
    matrix.write_text('2\nS1 0 0.125\nS2 0.125 0\n')
    values = _distance_matrix(matrix, [])
    assert values[1] == {'sample_1': 'S1', 'sample_2': 'S2', 'distance': 0.125}
    matrix.write_text('2\nS1 0\n')
    assert _distance_matrix(matrix, warnings) == []


def _minimal_run(tmp_path, mode='FO002', state='RUNNING'):
    run = tmp_path / 'runs' / 'test-run'
    run.mkdir(parents=True)
    effective = {'hash': 'frozen-effective-hash', 'config': {'project_id':'test', 'samples_tsv':'samples.tsv',
        'purpose':'analysis', 'final_output':mode, 'input_type':'FASTQ', 'assay':{'design':'WGS'},
        'reference':{'fasta':'missing-reference.fasta'}, 'analyses':{}}, 'parameters':{'layer6.max_embedded_records':{'value':1}},
        'quality':{'policy_id':'test-policy'}, 'compute':{}}
    records = {'frozen_config.json':effective, 'frozen_plan.json':{'samples':[{'sample_id':'S1'}],
        'final_output':mode,'reference':{},'layers':[]}, 'status.json':{'status':state},
        'snapshot.json':{'files':{}}}
    for name, value in records.items():
        (run / name).write_text(json.dumps(value))
    gate = run / 'results/qc/samples/S1/qc_decision.json'
    gate.parent.mkdir(parents=True)
    gate.write_text(json.dumps({'verdict':'PASS', 'metrics':{}}))
    return run


def test_completion_requires_final_artifacts_and_no_fake_variant_zero(tmp_path):
    from pfvarx.reporting.build import harvest
    run = _minimal_run(tmp_path, mode='FO001')
    data = harvest(tmp_path, run)
    assert data['samples'][0]['processing'] == 'RUNNING'
    assert data['variants_total'] is None
    bam = run / 'results/alignments/S1.bam'
    bam.parent.mkdir(parents=True)
    bam.touch()
    assert harvest(tmp_path, run)['samples'][0]['processing'] == 'COMPLETED'


def test_full_variant_summary_precedes_embedding_cap(tmp_path):
    from pfvarx.reporting.build import harvest
    run = _minimal_run(tmp_path)
    evidence = run / 'results/variants/samples/S1.evidence.tsv'
    evidence.parent.mkdir(parents=True)
    evidence.write_text('chrom\tpos\tref\talt\tkind\tstatus\twsaf\twsaf_lo\twsaf_hi\tn\ty\tscreen\n'
        'chr1\t1\tA\tG\tSNV\tPASS\t0.5\t0.2\t0.8\t10\t5\tinf\n'
        'chr1\t2\tT\tTA\tINS\tLOW_CONFIDENCE\tNA\tNA\tNA\t2\t1\t0\n')
    data = harvest(tmp_path, run)
    assert data['variants_total'] == 2
    assert len(data['variants']) == 1
    assert data['variant_summary']['types'] == {'SNV':1, 'INS':1}
    assert data['variants'][0]['screen'] == 'inf'
    assert data['variants'][0]['quality'] is None


def test_report_keeps_the_pipeline_version_that_produced_the_run(tmp_path):
    """Rebuilding an old run's report after an upgrade must not relabel the run."""
    from pfvarx import __version__
    from pfvarx.reporting.build import harvest
    run = _minimal_run(tmp_path, mode='FO001')
    (run / 'snapshot.json').write_text(json.dumps({'files': {}, 'pfvarx_version': '9.8.7',
                                                   'tools_installed': {'samtools': 'samtools 1.24'},
                                                   'packages': {'fgbio': '4.1.1 (hdfd78af_0)'}}))
    software = harvest(tmp_path, run)['software']
    assert software['pfvarx_version'] == '9.8.7' != __version__
    assert software['tools_installed'] == {'samtools': 'samtools 1.24', 'fgbio (package)': '4.1.1 (hdfd78af_0)'}


def test_snapshot_records_scientific_package_versions(tmp_path):
    from pfvarx.provenance.snapshot import package_versions
    meta = tmp_path / 'conda-meta'
    meta.mkdir()
    for name, version, build in (('fgbio', '4.1.1', 'hdfd78af_0'), ('r-moire', '3.7.0', 'r45h_0'), ('zlib', '1.3', 'x')):
        (meta / f'{name}-{version}-{build}.json').write_text(json.dumps({'name': name, 'version': version, 'build': build}))
    (meta / 'broken.json').write_text('{not json')
    assert package_versions(tmp_path) == {'fgbio': '4.1.1 (hdfd78af_0)', 'r-moire': '3.7.0 (r45h_0)'}
    assert package_versions(tmp_path / 'absent') == {}
