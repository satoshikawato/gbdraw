from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"


def test_collinearity_popup_uses_display_ids_and_hides_internal_rows(tmp_path: Path) -> None:
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not available")

    feature_utils_path = tmp_path / "feature-utils.mjs"
    feature_utils_path.write_text(
        (WEB_ROOT / "js" / "app" / "feature-utils.js").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    sequence_fasta_path = tmp_path / "feature-sequence-fasta.mjs"
    sequence_fasta_path.write_text(
        (WEB_ROOT / "js" / "app" / "feature-sequence-fasta.js")
        .read_text(encoding="utf-8")
        .replace("./feature-utils.js", "./feature-utils.mjs"),
        encoding="utf-8",
    )
    losat_normalization_path = tmp_path / "losat-normalization.mjs"
    losat_normalization_path.write_text(
        (WEB_ROOT / "js" / "app" / "losat-normalization.js").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    source_path = WEB_ROOT / "js" / "app" / "pairwise-match-popup.js"
    module_path = tmp_path / "pairwise-match-popup.mjs"
    module_path.write_text(
        source_path.read_text(encoding="utf-8")
        .replace("./feature-utils.js", "./feature-utils.mjs")
        .replace("./feature-sequence-fasta.js", "./feature-sequence-fasta.mjs")
        .replace("./losat-normalization.js", "./losat-normalization.mjs"),
        encoding="utf-8",
    )
    check_path = tmp_path / "check-collinearity-popup.mjs"
    check_path.write_text(
        f"""
        import {{ buildPairwiseMatchHoverRows, buildPairwiseMatchPayload }} from {module_path.as_uri()!r};

        const assert = (condition, message) => {{
          if (!condition) throw new Error(message);
        }};
        const rowValue = (section, label) => {{
          const row = section.rows.find((entry) => entry.label === label || entry.label === `${{label}}s`);
          return row ? row.value : '';
        }};

        const attrs = new Map(Object.entries({{
          'data-gbdraw-pairwise-match-id': 'comparison3_match3',
          'data-match-kind': 'collinear',
          'data-collinearity-block-id': 'block_0001',
          'data-collinearity-block-kind': 'cluster',
          'data-collinearity-orientation': 'plus',
          'data-collinearity-color-mode': 'orientation_identity',
          'data-collinearity-block-score': '25796',
          'data-collinearity-anchor-index': '1',
          'data-collinearity-anchor-count': '2',
          'data-orthogroup-id': 'og_1;og_2',
          'data-query-record-id': 'record_a',
          'data-subject-record-id': 'record_b',
          'data-qstart': '10',
          'data-qend': '80',
          'data-sstart': '100',
          'data-send': '180',
          'data-query-feature-svg-id': 'fq1;fq2',
          'data-subject-feature-svg-id': 'fs1;fs2',
          'data-query-protein-id': 'p_record_a_10_40_1_abcdef123456;p_record_a_50_80_1_123456abcdef',
          'data-subject-protein-id': 'gbd_r0002_cds000001;gbd_r0002_cds000002',
          'data-query-unit-id': 'gbd_r0001_unit000001;gbd_r0001_unit000002',
          'data-subject-unit-id': 'gbd_r0002_unit000003;gbd_r0002_unit000004',
          'data-query-locus-id': 'qlocus1;qlocus2',
          'data-subject-locus-id': 'slocus1;slocus2',
          'data-query-display-name': 'qgene1;qgene2',
          'data-subject-display-name': 'subject display 1;subject display 2',
          'data-identity': '76.10',
          'data-alignment-length': '17298'
        }}));
        const element = {{
          style: {{}},
          getAttribute: (name) => attrs.get(name) || ''
        }};
        const featureLookup = new Map([
          ['fq1', {{
            svg_id: 'fq1',
            record_id: 'record_a',
            proteinId: 'p_record_a_10_40_1_abcdef123456',
            sourceProteinId: 'WP_000001.1',
            qualifiers: {{ protein_id: ['WP_000001.1'] }},
            product: 'query product 1'
          }}],
          ['fq2', {{
            svg_id: 'fq2',
            record_id: 'record_a',
            proteinId: 'p_record_a_50_80_1_123456abcdef',
            qualifiers: {{ protein_id: ['WP_000002.1'] }},
            product: 'query product 2'
          }}],
          ['fs1', {{
            svg_id: 'fs1',
            record_id: 'record_b',
            proteinId: 'gbd_r0002_cds000001',
            locus_tag: 'SLOCUS_001',
            product: 'subject product'
          }}],
          ['fs2', {{
            svg_id: 'fs2',
            record_id: 'record_b',
            proteinId: 'gbd_r0002_cds000002',
            locus_tag: 'SLOCUS_002',
            product: 'subject product 2'
          }}]
        ]);

        const payload = buildPairwiseMatchPayload(element, {{
          featureLookup,
          orthogroups: [{{
            id: 'og_1',
            name: 'orthogroup display',
            members: [
              {{ recordId: 'record_a', featureSvgId: 'fq1', sourceProteinId: 'WP_000001.1' }},
              {{ recordId: 'record_b', featureSvgId: 'fs1', proteinId: 'gbd_r0002_cds000001', locusTag: 'SLOCUS_001' }}
            ]
          }}, {{
            id: 'og_2',
            name: 'orthogroup display 2',
            members: [
              {{ recordId: 'record_a', featureSvgId: 'fq2', sourceProteinId: 'WP_000002.1' }},
              {{ recordId: 'record_b', featureSvgId: 'fs2', proteinId: 'gbd_r0002_cds000002', locusTag: 'SLOCUS_002' }}
            ]
          }}]
        }});

        const sectionTitles = payload.sections.map((section) => section.title);
        assert(!sectionTitles.includes('Alignment'), `Alignment section leaked: ${{JSON.stringify(sectionTitles)}}`);
        assert(!sectionTitles.includes('Orthogroup'), `Single-OG section leaked: ${{JSON.stringify(sectionTitles)}}`);
        assert(sectionTitles.includes('Orthogroups covered'), `Block OG section missing: ${{JSON.stringify(sectionTitles)}}`);
        assert(sectionTitles.includes('Query'), `Query section missing: ${{JSON.stringify(sectionTitles)}}`);
        assert(sectionTitles.includes('Subject'), `Subject section missing: ${{JSON.stringify(sectionTitles)}}`);
        assert(!sectionTitles.includes('Query feature'), `Old Query feature title leaked: ${{JSON.stringify(sectionTitles)}}`);
        assert(!sectionTitles.includes('Subject feature'), `Old Subject feature title leaked: ${{JSON.stringify(sectionTitles)}}`);

        const labels = payload.sections.flatMap((section) => section.rows.map((row) => row.label));
        for (const forbidden of ['Match ID', 'Unit ID', 'Query unit', 'Subject unit']) {{
          assert(!labels.includes(forbidden), `${{forbidden}} leaked: ${{JSON.stringify(labels)}}`);
        }}
        assert(!labels.includes('Orthogroup ID'), `Long block Orthogroup ID leaked: ${{JSON.stringify(labels)}}`);
        assert(labels.includes('Number of orthogroups covered'), `Orthogroup count missing: ${{JSON.stringify(labels)}}`);

        const query = payload.sections.find((section) => section.title === 'Query');
        const subject = payload.sections.find((section) => section.title === 'Subject');
        assert(rowValue(query, 'Protein ID') === 'WP_000001.1; WP_000002.1', JSON.stringify(query.rows));
        assert(!rowValue(query, 'Protein ID').includes('p_record_a_'), JSON.stringify(query.rows));
        assert(rowValue(subject, 'Protein ID') === 'SLOCUS_001; SLOCUS_002', JSON.stringify(subject.rows));
        assert(!rowValue(subject, 'Protein ID').includes('gbd_r0002_cds000001'), JSON.stringify(subject.rows));
        assert(query.featureRows.length === 2, JSON.stringify(query.featureRows));
        assert(query.featureRows.every((row) => row.canOpen), JSON.stringify(query.featureRows));
        assert(query.featureRows.map((row) => row.label).join(',') === 'WP_000001.1,WP_000002.1', JSON.stringify(query.featureRows));
        assert(query.featureRows[0].subLabel === 'qlocus1 / qgene1', JSON.stringify(query.featureRows[0]));
        assert(query.featureRows[0].record === 'record_a', JSON.stringify(query.featureRows[0]));
        assert(query.featureRows[0].product === 'query product 1', JSON.stringify(query.featureRows[0]));
        assert(subject.featureRows.length === 2, JSON.stringify(subject.featureRows));
        assert(subject.featureRows.map((row) => row.label).join(',') === 'SLOCUS_001,SLOCUS_002', JSON.stringify(subject.featureRows));
        assert(subject.featureRows[0].subLabel === 'subject display 1', JSON.stringify(subject.featureRows[0]));
        assert(subject.featureRows[0].feature.svg_id === 'fs1', JSON.stringify(subject.featureRows[0]));

        const duplicateAttrs = new Map(Object.entries({{
          'data-gbdraw-pairwise-match-id': 'comparison3_match4',
          'data-match-kind': 'collinear',
          'data-collinearity-block-id': 'block_0002',
          'data-query-record-id': 'AP027131.1',
          'data-qstart': '22946',
          'data-qend': '24703',
          'data-query-feature-svg-id': 'fq_dup',
          'data-query-protein-id': 'BDV02135.1',
          'data-query-locus-id': 'HPAVJP_0240',
          'data-query-display-name': 'HPAVJP_0240'
        }}));
        const duplicatePayload = buildPairwiseMatchPayload({{
          style: {{}},
          getAttribute: (name) => duplicateAttrs.get(name) || ''
        }}, {{
          featureLookup: new Map([['fq_dup', {{
            svg_id: 'fq_dup',
            record_id: 'AP027131.1',
            start: 22945,
            end: 24703,
            strand: '+',
            sourceProteinId: 'BDV02135.1',
            gene: 'HPAVJP_0240',
            locus_tag: 'HPAVJP_0240',
            product: 'DHH family phosphoesterase',
            qualifiers: {{
              protein_id: ['BDV02135.1'],
              gene: ['HPAVJP_0240'],
              locus_tag: ['HPAVJP_0240']
            }}
          }}]])
        }});
        const duplicateQuery = duplicatePayload.sections.find((section) => section.title === 'Query');
        assert(duplicateQuery.featureRows[0].label === 'BDV02135.1', JSON.stringify(duplicateQuery.featureRows[0]));
        assert(duplicateQuery.featureRows[0].subLabel === 'HPAVJP_0240', JSON.stringify(duplicateQuery.featureRows[0]));
        assert(!duplicateQuery.featureRows[0].subLabel.includes(' / HPAVJP_0240'), JSON.stringify(duplicateQuery.featureRows[0]));

        const blockOgSection = payload.sections.find((section) => section.title === 'Orthogroups covered');
        assert(rowValue(blockOgSection, 'Number of orthogroups covered') === '2', JSON.stringify(blockOgSection.rows));
        assert(payload.blockOrthogroupCount === 2, JSON.stringify(payload));
        assert(payload.blockOrthogroups.length === 2, JSON.stringify(payload.blockOrthogroups));
        assert(payload.blockOrthogroups[0].id === 'og_1', JSON.stringify(payload.blockOrthogroups));
        assert(payload.blockOrthogroups[0].displayName === 'orthogroup display', JSON.stringify(payload.blockOrthogroups[0]));
        assert(payload.blockOrthogroups[0].queryMember === 'WP_000001.1', JSON.stringify(payload.blockOrthogroups[0]));
        assert(payload.blockOrthogroups[0].subjectMember === 'SLOCUS_001', JSON.stringify(payload.blockOrthogroups[0]));
        assert(payload.blockOrthogroups[0].detailRows.some((row) => row.label === 'Orthogroup ID' && row.value === 'og_1'), JSON.stringify(payload.blockOrthogroups[0].detailRows));
        assert(payload.blockOrthogroups[1].id === 'og_2', JSON.stringify(payload.blockOrthogroups));
        assert(payload.blockOrthogroups[1].queryMember === 'WP_000002.1', JSON.stringify(payload.blockOrthogroups[1]));
        assert(payload.blockOrthogroups[1].subjectMember === 'SLOCUS_002', JSON.stringify(payload.blockOrthogroups[1]));

        const hoverLabels = buildPairwiseMatchHoverRows(payload).map((row) => row.label);
        assert(hoverLabels.includes('Identity'), JSON.stringify(hoverLabels));
        assert(hoverLabels.includes('Orthogroups'), JSON.stringify(hoverLabels));
        assert(!hoverLabels.includes('Orthogroup'), JSON.stringify(hoverLabels));
        """,
        encoding="utf-8",
    )

    subprocess.run([node, str(check_path)], check=True, cwd=REPO_ROOT)
