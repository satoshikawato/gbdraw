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
    color_utils_path = tmp_path / "color-utils.mjs"
    color_utils_path.write_text(
        (WEB_ROOT / "js" / "app" / "color-utils.js").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    conservation_series_path = tmp_path / "conservation-series.mjs"
    conservation_series_path.write_text(
        (WEB_ROOT / "js" / "app" / "conservation-series.js")
        .read_text(encoding="utf-8")
        .replace("./color-utils.js", "./color-utils.mjs"),
        encoding="utf-8",
    )
    match_sequences_path = tmp_path / "match-sequences.mjs"
    match_sequences_path.write_text(
        (WEB_ROOT / "js" / "app" / "match-sequences.js")
        .read_text(encoding="utf-8")
        .replace("./feature-sequence-fasta.js", "./feature-sequence-fasta.mjs")
        .replace("./conservation-series.js", "./conservation-series.mjs"),
        encoding="utf-8",
    )
    losat_normalization_path = tmp_path / "losat-normalization.mjs"
    losat_normalization_path.write_text(
        (WEB_ROOT / "js" / "app" / "losat-normalization.js").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    feature_identity_path = tmp_path / "feature-identity.mjs"
    feature_identity_path.write_text(
        (WEB_ROOT / "js" / "services" / "feature-identity.js").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    source_path = WEB_ROOT / "js" / "app" / "pairwise-match-popup.js"
    module_path = tmp_path / "pairwise-match-popup.mjs"
    module_path.write_text(
        source_path.read_text(encoding="utf-8")
        .replace("./feature-utils.js", "./feature-utils.mjs")
        .replace("./feature-sequence-fasta.js", "./feature-sequence-fasta.mjs")
        .replace("./match-sequences.js", "./match-sequences.mjs")
        .replace("./losat-normalization.js", "./losat-normalization.mjs")
        .replace("../services/feature-identity.js", "./feature-identity.mjs")
        + "\nexport { buildFallbackOrthogroup, featureOrthogroupId, getFeatureForMember, getRenderedFeatureForMember, getGroupMemberForFeatureSvgId, getOrthogroupById, getOrthogroupForMatch, integerAttr };\n",
        encoding="utf-8",
    )
    check_path = tmp_path / "check-collinearity-popup.mjs"
    check_path.write_text(
        f"""
        import {{
          buildFallbackOrthogroup,
          buildPairwiseMatchHoverRows,
          buildPairwiseMatchPayload,
          featureOrthogroupId,
          getFeatureForMember,
          getRenderedFeatureForMember,
          getGroupMemberForFeatureSvgId,
          getOrthogroupById,
          getOrthogroupForMatch,
          integerAttr
        }} from {module_path.as_uri()!r};

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
          'data-query-record-index': '0',
          'data-subject-record-index': '1',
          'data-qstart': '10',
          'data-qend': '80',
          'data-sstart': '100',
          'data-send': '180',
          'data-query-feature-svg-id': 'fq1;fq2',
          'data-subject-feature-svg-id': 'fs1;fs2',
          'data-query-stable-feature-svg-id': 'fq1;fq2',
          'data-subject-stable-feature-svg-id': 'fs1;fs2',
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
            fileIdx: 0,
            svg_id: 'fq1',
            stable_feature_id: 'fq1',
            record_id: 'record_a',
            proteinId: 'p_record_a_10_40_1_abcdef123456',
            sourceProteinId: 'WP_000001.1',
            qualifiers: {{ protein_id: ['WP_000001.1'] }},
            product: 'query product 1'
          }}],
          ['fq2', {{
            fileIdx: 0,
            svg_id: 'fq2',
            stable_feature_id: 'fq2',
            record_id: 'record_a',
            proteinId: 'p_record_a_50_80_1_123456abcdef',
            qualifiers: {{ protein_id: ['WP_000002.1'] }},
            product: 'query product 2'
          }}],
          ['fs1', {{
            fileIdx: 1,
            svg_id: 'fs1',
            stable_feature_id: 'fs1',
            record_id: 'record_b',
            proteinId: 'gbd_r0002_cds000001',
            locus_tag: 'SLOCUS_001',
            product: 'subject product'
          }}],
          ['fs2', {{
            fileIdx: 1,
            svg_id: 'fs2',
            stable_feature_id: 'fs2',
            record_id: 'record_b',
            proteinId: 'gbd_r0002_cds000002',
            locus_tag: 'SLOCUS_002',
            product: 'subject product 2'
          }}],
          ['display-space-rendered-id', {{
            fileIdx: 4,
            feature_index: 25,
            recordKey: 'record-key-e',
            biologicalFeatureId: 'bio-e',
            svg_id: 'display-space-rendered-id',
            stable_feature_id: 'source-space-stable-id',
            record_id: 'record_e',
            sourceProteinId: 'CAG34720.1'
          }}]
        ]);

        const popupOptions = {{
          featureLookup,
          sourceFeatures: [{{
            fileIdx: 4,
            feature_index: 25,
            recordKey: 'record-key-e',
            biologicalFeatureId: 'bio-e',
            record_id: 'record_e',
            stable_feature_id: 'source-space-stable-id',
            source_protein_id: 'CAG34720.1',
            qualifiers: {{ protein_id: ['CAG34720.1'] }},
            nucleotide_sequence: 'ATGGGGTAA',
            amino_acid_sequence: 'MG'
          }}],
          orthogroups: [{{
            id: 'og_1',
            name: 'orthogroup display',
            members: [
              {{ recordIndex: 0, recordId: 'record_a', featureSvgId: 'fq1', sourceProteinId: 'WP_000001.1' }},
              {{ recordIndex: 1, recordId: 'record_b', featureSvgId: 'fs1', proteinId: 'gbd_r0002_cds000001', locusTag: 'SLOCUS_001' }},
              {{
                recordId: 'record_e',
                recordIndex: 4,
                featureIndex: 25,
                recordKey: 'record-key-e',
                biologicalFeatureId: 'bio-e',
                featureSvgId: 'source-space-stable-id',
                stableFeatureSvgId: 'source-space-stable-id',
                renderedFeatureSvgId: 'display-space-rendered-id',
                sourceProteinId: 'CAG34720.1'
              }}
            ]
          }}, {{
            id: 'og_2',
            name: 'orthogroup display 2',
            members: [
              {{ recordIndex: 0, recordId: 'record_a', featureSvgId: 'fq2', sourceProteinId: 'WP_000002.1' }},
              {{ recordIndex: 1, recordId: 'record_b', featureSvgId: 'fs2', proteinId: 'gbd_r0002_cds000002', locusTag: 'SLOCUS_002' }}
            ]
          }}]
        }};
        const conflictingGroupAliases = {{
          id: 'og_conflict_a',
          orthogroupId: 'og_conflict_b',
          members: [{{ recordIndex: 0, stableFeatureSvgId: 'fq1' }}]
        }};
        assert(
          getOrthogroupById([conflictingGroupAliases], 'og_conflict_a') === null,
          'conflicting group aliases resolved by first value'
        );
        assert(
          getOrthogroupById([
            {{ id: 'og_duplicate' }},
            {{ orthogroupId: 'og_duplicate' }}
          ], 'og_duplicate') === null,
          'duplicate group IDs resolved by first group'
        );
        const scopedGlobalGroup = {{ id: 'og_scoped', scope: 'global' }};
        const scopedLocalGroup = {{
          id: 'og_scoped',
          scope: 'cross_record',
          presentationScope: 'adjacent_local'
        }};
        assert(
          getOrthogroupById(
            [scopedGlobalGroup, scopedLocalGroup],
            'og_scoped',
            'adjacent_local'
          ) === scopedLocalGroup,
          'local group did not resolve inside its metadata namespace'
        );
        assert(
          getOrthogroupById(
            [scopedGlobalGroup, scopedLocalGroup],
            'og_scoped',
            'global'
          ) === scopedGlobalGroup,
          'global group did not resolve inside its metadata namespace'
        );
        assert(
          getOrthogroupById(
            [scopedLocalGroup, {{ ...scopedLocalGroup }}],
            'og_scoped',
            'adjacent_local'
          ) === null,
          'duplicate local group IDs did not fail closed'
        );
        assert(
          getOrthogroupById([{{
            ...scopedLocalGroup,
            collinearGroupScope: 'global_collinear'
          }}], 'og_scoped', 'adjacent_local') === null,
          'conflicting presentation-scope aliases did not fail closed'
        );
        assert(
          getOrthogroupForMatch([conflictingGroupAliases], {{
            orthogroupId: 'og_conflict_a',
            queryFeatureSvgId: 'fq1',
            featureLookup
          }}) === null,
          'explicit conflicting group ID fell back to member lookup'
        );
        const conflictingFeatureGroupAliases = {{
          ...featureLookup.get('fq1'),
          orthogroupId: 'og_conflict_a',
          orthogroup_id: 'og_conflict_b'
        }};
        assert(
          featureOrthogroupId(conflictingFeatureGroupAliases) === '',
          'conflicting feature group aliases resolved by first value'
        );
        assert(
          buildFallbackOrthogroup({{
            orthogroupId: 'og_conflict_a',
            queryFeature: conflictingFeatureGroupAliases,
            subjectFeature: null,
            featureLookup: new Map([['fq1', conflictingFeatureGroupAliases]])
          }}) === null,
          'conflicting feature group aliases built a fallback group'
        );
        const strictIntegerElement = (value) => ({{
          getAttribute: (name) => name === 'data-index' ? value : ''
        }});
        assert(integerAttr(strictIntegerElement('4'), 'data-index') === 4, 'integerAttr rejected 4');
        for (const invalidInteger of [
          '-1',
          '4.5',
          '9007199254740992'
        ]) {{
          assert(
            integerAttr(strictIntegerElement(invalidInteger), 'data-index') === null,
            `integerAttr accepted ${{invalidInteger}}`
          );
        }}

        const invalidSequenceAttrs = new Map(Object.entries({{
          'data-gbdraw-pairwise-match-id': 'invalid-sequence-identity',
          'data-match-kind': 'pairwise',
          'data-query-record-id': 'recA',
          'data-query-record-index': '-1',
          'data-qstart': '1',
          'data-qend': '4'
        }}));
        let invalidQueryResolverCalls = 0;
        const invalidSequencePayload = buildPairwiseMatchPayload({{
          style: {{}},
          getAttribute: (name) => invalidSequenceAttrs.get(name) || ''
        }}, {{
          resolveSequenceSource: (_sourceKey, recordId) => {{
            if (recordId !== 'recA') return null;
            invalidQueryResolverCalls += 1;
            return {{
              key: 'linear:record:0',
              recordId: 'recA',
              sequence: 'AAAA',
              origin: 'linear-record',
              recordIndex: 0
            }};
          }}
        }});
        assert(
          invalidSequencePayload.sequenceBundle.entries[0].available === false,
          JSON.stringify(invalidSequencePayload.sequenceBundle.entries[0])
        );
        assert(
          invalidQueryResolverCalls === 0,
          'invalid record index reached the sequence resolver'
        );

        const endpointPayload = (role, overrides = {{}}) => {{
          const roleTitle = role === 'query' ? 'Query' : 'Subject';
          const endpointAttrs = new Map(Object.entries({{
            'data-gbdraw-pairwise-match-id': `endpoint-${{role}}`,
            'data-match-kind': 'pairwise',
            [`data-${{role}}-record-id`]: 'record_e',
            [`data-${{role}}-record-index`]: '4',
            [`data-${{role}}-feature-svg-id`]: 'display-space-rendered-id',
            [`data-${{role}}-stable-feature-svg-id`]: 'source-space-stable-id',
            [`data-${{role}}-feature-index`]: '25',
            ...(role === 'query'
              ? {{ 'data-qstart': '1', 'data-qend': '9' }}
              : {{ 'data-sstart': '1', 'data-send': '9' }}),
            ...overrides
          }}));
          const endpointResult = buildPairwiseMatchPayload({{
            style: {{}},
            getAttribute: (name) => endpointAttrs.get(name) || ''
          }}, {{ featureLookup, sourceFeatures: popupOptions.sourceFeatures }});
          return endpointResult.sections.find((section) => section.title === roleTitle)?.featureRows?.[0];
        }};
        for (const role of ['query', 'subject']) {{
          const validEndpoint = endpointPayload(role);
          assert(validEndpoint?.canOpen, `valid ${{role}} endpoint did not open`);
          assert(
            validEndpoint.feature?.svg_id === 'display-space-rendered-id',
            `valid ${{role}} endpoint resolved incorrectly`
          );
          const prefix = `data-${{role}}`;
          for (const invalidAttrs of [
            {{
              [`${{prefix}}-record-index`]: '',
              [`${{prefix}}-stable-feature-svg-id`]: '',
              [`${{prefix}}-feature-index`]: ''
            }},
            {{
              [`${{prefix}}-stable-feature-svg-id`]: '',
              [`${{prefix}}-feature-index`]: ''
            }},
            {{ [`${{prefix}}-stable-feature-svg-id`]: 'wrong-stable-id' }},
            {{ [`${{prefix}}-feature-index`]: '24' }},
            {{ [`${{prefix}}-record-index`]: '3' }},
            {{ [`${{prefix}}-record-index`]: '-1' }},
            {{ [`${{prefix}}-record-index`]: '4.5' }},
            {{ [`${{prefix}}-record-index`]: '9007199254740992' }},
            {{ [`${{prefix}}-feature-index`]: '-1' }},
            {{ [`${{prefix}}-feature-index`]: '25.5' }},
            {{ [`${{prefix}}-feature-index`]: '9007199254740992' }}
          ]) {{
            const rejectedEndpoint = endpointPayload(role, invalidAttrs);
            assert(
              rejectedEndpoint && !rejectedEndpoint.canOpen && rejectedEndpoint.feature === null,
              `invalid ${{role}} endpoint opened: ${{JSON.stringify(invalidAttrs)}}`
            );
          }}
        }}
        const canonicalMember = popupOptions.orthogroups[0].members[2];
        const canonicalSource = popupOptions.sourceFeatures[0];
        assert(
          getFeatureForMember(canonicalMember, featureLookup, popupOptions.sourceFeatures) === canonicalSource,
          'complete canonical/source identity did not resolve'
        );
        assert(
          getFeatureForMember({{
            recordKey: 'record-key-e',
            biologicalFeatureId: 'bio-e'
          }}, featureLookup, popupOptions.sourceFeatures) === canonicalSource,
          'schema-3 canonical identity did not resolve'
        );
        const sourceIdentityMember = {{
          recordIndex: 4,
          featureIndex: 25,
          stableFeatureSvgId: 'source-space-stable-id',
          renderedFeatureSvgId: 'display-space-rendered-id',
          sourceProteinId: 'CAG34720.1'
        }};
        assert(
          getFeatureForMember(
            sourceIdentityMember,
            featureLookup,
            popupOptions.sourceFeatures
          ) === canonicalSource,
          'record-scoped LOSATP source identity did not join to the compact catalog feature'
        );
        for (const invalidMember of [
          {{ stableFeatureSvgId: 'source-space-stable-id' }},
          {{ recordIndex: -1, stableFeatureSvgId: 'source-space-stable-id' }},
          {{ recordIndex: '4.5', stableFeatureSvgId: 'source-space-stable-id' }},
          {{ recordIndex: 4, record_index: 3, stableFeatureSvgId: 'source-space-stable-id' }},
          {{
            recordIndex: 4,
            stableFeatureSvgId: 'source-space-stable-id',
            featureSvgId: 'wrong-stable-id'
          }},
          {{
            recordKey: 'record-key-e',
            biologicalFeatureId: 'wrong-bio-id'
          }},
          {{
            recordIndex: 4,
            featureIndex: 24,
            stableFeatureSvgId: 'source-space-stable-id'
          }}
        ]) {{
          assert(
            getFeatureForMember(invalidMember, featureLookup, popupOptions.sourceFeatures) === null,
            `invalid member identity resolved: ${{JSON.stringify(invalidMember)}}`
          );
        }}
        assert(
          getFeatureForMember(
            canonicalMember,
            featureLookup,
            [canonicalSource, {{ ...canonicalSource }}]
          ) === null,
          'duplicate biological candidates did not fail closed'
        );
        assert(
          getFeatureForMember(
            {{ ...canonicalMember, renderedFeatureSvgId: 'missing-rendered-id' }},
            featureLookup,
            popupOptions.sourceFeatures
          ) === null,
          'unknown rendered alias was ignored during source resolution'
        );
        assert(
          getFeatureForMember(
            {{ ...canonicalMember, renderedFeatureSvgId: '' }},
            featureLookup,
            [{{ ...canonicalSource, rendered_svg_id: 'missing-rendered-id' }}]
          ) === null,
          'unknown source rendered alias was ignored'
        );
        assert(
          getRenderedFeatureForMember(canonicalMember, canonicalSource, featureLookup)?.svg_id ===
            'display-space-rendered-id',
          'canonical rendered endpoint did not resolve'
        );
        const duplicateLocationSources = [4, 7].map((sourceFeatureIndex) => ({{
          fileIdx: 0,
          sourceFeatureIndex,
          stable_feature_id: 'same-location-stable-id',
          record_id: 'same-location-record'
        }}));
        const duplicateLocationLookup = new Map([4, 7].map((sourceFeatureIndex) => [
          `same-location-rendered-${{sourceFeatureIndex}}`,
          {{
            fileIdx: 0,
            sourceFeatureIndex,
            stable_feature_id: 'same-location-stable-id',
            svg_id: `same-location-rendered-${{sourceFeatureIndex}}`,
            record_id: 'same-location-record'
          }}
        ]));
        const duplicateLocationMembers = [4, 7].map((featureIndex) => ({{
          recordIndex: 0,
          featureIndex,
          stableFeatureSvgId: 'same-location-stable-id',
          renderedFeatureSvgId: `same-location-rendered-${{featureIndex}}`
        }}));
        duplicateLocationMembers.forEach((member, index) => {{
          assert(
            getFeatureForMember(
              member,
              duplicateLocationLookup,
              duplicateLocationSources
            ) === duplicateLocationSources[index],
            'sourceFeatureIndex did not survive duplicate-location source resolution'
          );
          assert(
            getRenderedFeatureForMember(
              member,
              duplicateLocationSources[index],
              duplicateLocationLookup
            )?.svg_id === `same-location-rendered-${{member.featureIndex}}`,
            'sourceFeatureIndex did not resolve the exact duplicate-location DOM feature'
          );
        }});
        assert(
          getFeatureForMember(
            {{ recordIndex: 0, stableFeatureSvgId: 'same-location-stable-id' }},
            duplicateLocationLookup,
            duplicateLocationSources
          ) === null,
          'duplicate-location member without a source index did not fail closed'
        );
        assert(
          getRenderedFeatureForMember(
            {{ ...canonicalMember, renderedFeatureSvgId: 'missing-rendered-id' }},
            canonicalSource,
            featureLookup
          ) === null,
          'unknown rendered alias was ignored'
        );
        const duplicateRenderedLookup = new Map(featureLookup);
        duplicateRenderedLookup.set('duplicate-rendered-object', {{
          ...featureLookup.get('display-space-rendered-id')
        }});
        assert(
          getRenderedFeatureForMember(canonicalMember, canonicalSource, duplicateRenderedLookup) === null,
          'duplicate rendered candidates did not fail closed'
        );
        const incompleteRenderedLookup = new Map(featureLookup);
        incompleteRenderedLookup.set('display-space-rendered-id', {{
          svg_id: 'display-space-rendered-id',
          recordKey: 'record-key-e',
          biologicalFeatureId: 'bio-e'
        }});
        assert(
          getRenderedFeatureForMember(canonicalMember, canonicalSource, incompleteRenderedLookup) === null,
          'rendered endpoint without record/stable/source agreement resolved'
        );
        assert(
          getGroupMemberForFeatureSvgId(
            popupOptions.orthogroups[0],
            'display-space-rendered-id',
            featureLookup
          ) === canonicalMember,
          'rendered endpoint did not map to its unique member'
        );
        assert(
          getGroupMemberForFeatureSvgId(
            {{
              members: [canonicalMember, {{ ...canonicalMember }}]
            }},
            'display-space-rendered-id',
            featureLookup
          ) === null,
          'duplicate group members did not fail closed'
        );
        const payload = buildPairwiseMatchPayload(element, popupOptions);

        const sectionTitles = payload.sections.map((section) => section.title);
        assert(!sectionTitles.includes('Alignment'), `Alignment section leaked: ${{JSON.stringify(sectionTitles)}}`);
        assert(!sectionTitles.includes('Similarity group'), `Single-group section leaked: ${{JSON.stringify(sectionTitles)}}`);
        assert(sectionTitles.includes('Similarity groups covered'), `Block group section missing: ${{JSON.stringify(sectionTitles)}}`);
        assert(sectionTitles.includes('Query'), `Query section missing: ${{JSON.stringify(sectionTitles)}}`);
        assert(sectionTitles.includes('Subject'), `Subject section missing: ${{JSON.stringify(sectionTitles)}}`);
        assert(!sectionTitles.includes('Query feature'), `Old Query feature title leaked: ${{JSON.stringify(sectionTitles)}}`);
        assert(!sectionTitles.includes('Subject feature'), `Old Subject feature title leaked: ${{JSON.stringify(sectionTitles)}}`);

        const labels = payload.sections.flatMap((section) => section.rows.map((row) => row.label));
        for (const forbidden of ['Match ID', 'Unit ID', 'Query unit', 'Subject unit']) {{
          assert(!labels.includes(forbidden), `${{forbidden}} leaked: ${{JSON.stringify(labels)}}`);
        }}
        assert(!labels.includes('Similarity group ID'), `Long block similarity group ID leaked: ${{JSON.stringify(labels)}}`);
        assert(labels.includes('Number of similarity groups covered'), `Similarity-group count missing: ${{JSON.stringify(labels)}}`);

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
          'data-query-record-index': '0',
          'data-qstart': '22946',
          'data-qend': '24703',
          'data-query-feature-svg-id': 'fq_dup',
          'data-query-stable-feature-svg-id': 'fq_dup',
          'data-query-protein-id': 'BDV02135.1',
          'data-query-locus-id': 'HPAVJP_0240',
          'data-query-display-name': 'HPAVJP_0240'
        }}));
        const duplicatePayload = buildPairwiseMatchPayload({{
          style: {{}},
          getAttribute: (name) => duplicateAttrs.get(name) || ''
        }}, {{
          featureLookup: new Map([['fq_dup', {{
            fileIdx: 0,
            svg_id: 'fq_dup',
            stable_feature_id: 'fq_dup',
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

        const runtimeOnlyHandle = 'h_aaaaaaaaaaaaaaaaaaaaaaaaaa';
        const runtimeOnlyAttrs = new Map(Object.entries({{
          'data-gbdraw-pairwise-match-id': 'comparison3_match5',
          'data-match-kind': 'pairwise',
          'data-query-record-id': 'record_runtime',
          'data-qstart': '1',
          'data-qend': '30',
          'data-query-protein-id': runtimeOnlyHandle
        }}));
        const runtimeOnlyPayload = buildPairwiseMatchPayload({{
          style: {{}},
          getAttribute: (name) => runtimeOnlyAttrs.get(name) || ''
        }}, {{ featureLookup: new Map() }});
        const runtimeOnlyQuery = runtimeOnlyPayload.sections.find((section) => section.title === 'Query');
        assert(runtimeOnlyQuery, JSON.stringify(runtimeOnlyPayload.sections));
        assert(!JSON.stringify(runtimeOnlyQuery).includes(runtimeOnlyHandle), JSON.stringify(runtimeOnlyQuery));

        const blockOgSection = payload.sections.find((section) => section.title === 'Similarity groups covered');
        assert(rowValue(blockOgSection, 'Number of similarity groups covered') === '2', JSON.stringify(blockOgSection.rows));
        assert(payload.blockOrthogroupCount === 2, JSON.stringify(payload));
        assert(payload.blockOrthogroups.length === 2, JSON.stringify(payload.blockOrthogroups));
        assert(payload.blockOrthogroups[0].id === 'og_1', JSON.stringify(payload.blockOrthogroups));
        assert(payload.blockOrthogroups[0].displayName === 'orthogroup display', JSON.stringify(payload.blockOrthogroups[0]));
        assert(payload.blockOrthogroups[0].queryMember === 'WP_000001.1', JSON.stringify(payload.blockOrthogroups[0]));
        assert(payload.blockOrthogroups[0].subjectMember === 'SLOCUS_001', JSON.stringify(payload.blockOrthogroups[0]));
        assert(payload.blockOrthogroups[0].detailRows.some((row) => row.label === 'Similarity group ID' && row.value === 'og_1'), JSON.stringify(payload.blockOrthogroups[0].detailRows));
        const reverseMember = payload.blockOrthogroups[0].memberRows.find((row) => row.proteinId === 'CAG34720.1');
        assert(reverseMember && reverseMember.aaFasta.includes('MG'), JSON.stringify(payload.blockOrthogroups[0].memberRows));
        assert(reverseMember.canOpen, JSON.stringify(reverseMember));
        assert(reverseMember.feature.svg_id === 'display-space-rendered-id', JSON.stringify(reverseMember));
        assert(payload.blockOrthogroups[0].memberAaFasta.includes('>CAG34720.1'), payload.blockOrthogroups[0].memberAaFasta);
        featureLookup.get('display-space-rendered-id').stable_feature_id = 'wrong-source-id';
        const mismatchedPayload = buildPairwiseMatchPayload(element, popupOptions);
        const mismatchedReverseMember = mismatchedPayload.blockOrthogroups[0].memberRows
          .find((row) => row.proteinId === 'CAG34720.1');
        assert(!mismatchedReverseMember.canOpen, JSON.stringify(mismatchedReverseMember));
        featureLookup.get('display-space-rendered-id').stable_feature_id = 'source-space-stable-id';
        featureLookup.set('source-space-stable-id', popupOptions.sourceFeatures[0]);
        popupOptions.orthogroups[0].members[2].recordIndex = 3;
        const wrongRecordPayload = buildPairwiseMatchPayload(element, popupOptions);
        const wrongRecordReverseMember = wrongRecordPayload.blockOrthogroups[0].memberRows
          .find((row) => row.proteinId === 'CAG34720.1');
        assert(!wrongRecordReverseMember.aaFasta, JSON.stringify(wrongRecordReverseMember));
        assert(!wrongRecordReverseMember.canOpen, JSON.stringify(wrongRecordReverseMember));
        popupOptions.orthogroups[0].members[2].recordIndex = 4;
        featureLookup.delete('source-space-stable-id');
        assert(payload.blockOrthogroups[1].id === 'og_2', JSON.stringify(payload.blockOrthogroups));
        assert(payload.blockOrthogroups[1].queryMember === 'WP_000002.1', JSON.stringify(payload.blockOrthogroups[1]));
        assert(payload.blockOrthogroups[1].subjectMember === 'SLOCUS_002', JSON.stringify(payload.blockOrthogroups[1]));

        const hoverLabels = buildPairwiseMatchHoverRows(payload).map((row) => row.label);
        assert(hoverLabels.includes('Identity'), JSON.stringify(hoverLabels));
        assert(hoverLabels.includes('Similarity groups'), JSON.stringify(hoverLabels));
        assert(!hoverLabels.includes('Similarity group'), JSON.stringify(hoverLabels));
        """,
        encoding="utf-8",
    )

    subprocess.run([node, str(check_path)], check=True, cwd=REPO_ROOT)
