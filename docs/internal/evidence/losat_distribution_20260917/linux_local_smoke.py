"""Local candidate acceptance; no public Release identity is claimed."""
import hashlib
import io
import json
from pathlib import Path
import sys
from unittest.mock import patch

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from gbdraw import losat_setup as setup
from gbdraw.analysis import protein_colinearity as protein, collinearity

root = Path('/tmp/losat-release-20260917')
repo = root / 'losat'
sys.path.insert(0, str(repo / 'LOSAT/tests'))
import prepare_release_candidate_v010 as release

output = root / 'evidence/gbdraw-linux-local'
output.mkdir(exist_ok=True)
target = setup.runtime_target()
binary = repo / 'LOSAT/target/release/LOSAT'
candidate = '0dd6c2b83f437a940d5fcbea5ca8a247934b48f3'
archive_root = f'LOSAT-0.1.0-{target}'
filename = archive_root + '.tar.gz'
binary_hash = setup.file_sha256(binary)
metadata = {'candidate_sha': candidate, 'release': 'v0.1.0',
            'artifact': {'target': target}, 'binary': {'sha256': binary_hash}}
payload = release.tar_bytes(release.archive_entries(repo, archive_root, 'LOSAT', binary, metadata))
(output / filename).write_bytes(payload)
entry = {'filename': filename, 'sha256': hashlib.sha256(payload).hexdigest(), 'size': len(payload),
         'binary': 'LOSAT', 'binary_sha256': binary_hash, 'binary_size': binary.stat().st_size}
lock = {'schema_version': 1, 'version': '0.1.0', 'candidate_sha': candidate, 'artifacts': {target: entry}}

def download(url, timeout):
    assert url == f'{setup.RELEASE_URL}/v0.1.0/{filename}'
    response = io.BytesIO(payload)
    response.geturl = lambda: url
    return response

report = {'scope': 'Local Linux candidate archive and real search integration; transport supplied by fixture, not public URL acceptance',
          'candidate_sha': candidate, 'target': target, 'archive': entry,
          'gbdraw_package': str(Path(setup.__file__).resolve()), 'cases': []}
query_path = repo / 'LOSAT/tests/fasta/SicyWSV.faa'
subject_path = repo / 'LOSAT/tests/fasta/PajaWSV.faa'
proteins = list(SeqIO.parse(query_path, 'fasta'))[:6]
records = []
for record_index in range(2):
    record = SeqRecord(Seq('N' * sum(len(p.seq)*3+30 for p in proteins)), id=f'genome{record_index}')
    offset = 0
    for index, p in enumerate(proteins):
        length = len(p.seq) * 3
        record.features.append(SeqFeature(FeatureLocation(offset, offset+length, strand=1), type='CDS',
            qualifiers={'translation': [str(p.seq)], 'protein_id': [f'protein{index}'], 'locus_tag': [f'gene{index}']}))
        offset += length + 30
    records.append(record)

with patch.object(setup, 'read_release_lock', return_value=lock), patch.object(setup, 'cache_root', return_value=output / 'cache'):
    with patch.object(setup, 'urlopen', side_effect=download):
        installed = setup.setup_losat()
    with patch.object(setup, 'urlopen', side_effect=AssertionError('network forbidden after first setup')):
        assert setup.setup_losat() == installed
        assert setup.managed_losat() == installed
        raw = []
        frame = protein.run_losatp_blastp(query_path.read_text(), subject_path.read_text(), threads=1, raw_output_callback=raw.append)
        digest = hashlib.sha256(raw[0].encode()).hexdigest()
        assert digest == 'fd4b010800e32ce6c823cb38b42a10b7845f3342edae892acccc8f554f9edf34', digest
        report['cases'].append({'case': 'release_blastp_smoke_offline', 'output_sha256': digest, 'rows': len(frame), 'status': 'PASS'})
        for name, builder in [('pairwise', protein.build_pairwise_protein_blastp_comparisons),
                              ('similarity_groups', protein.build_rbh_orthogroup_protein_blastp_comparisons),
                              ('collinear', collinearity.build_orthogroup_collinearity_blocks)]:
            managed = builder(records, losatp_threads=1)
            explicit = builder(records, losatp_threads=1, losatp_bin=str(binary))
            actual_frames = collinearity.convert_collinearity_blocks_to_comparisons(managed, records=records) if name == 'collinear' else managed.comparisons
            expected_frames = collinearity.convert_collinearity_blocks_to_comparisons(explicit, records=records) if name == 'collinear' else explicit.comparisons
            assert len(actual_frames) == len(expected_frames)
            for actual, expected in zip(actual_frames, expected_frames):
                assert actual.equals(expected), name
            rows = sum(len(x) for x in actual_frames)
            assert rows > 0, name
            report['cases'].append({'case': name, 'rows': rows, 'managed_equals_explicit': True, 'status': 'PASS'})
report['status'] = 'PASS'
(output / 'result.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
