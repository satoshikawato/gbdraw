import {
  DIAGRAM_HELPER_OPERATIONS,
  runDiagramHelperOperation
} from '../services/diagram-generation.js';
import {
  cloneFileBytesForTransfer,
  getSessionResourceSource,
  readFileText
} from '../services/file-content-cache.js';

const normalizeRecordLength = (value) => {
  const numeric = Number(value);
  return Number.isInteger(numeric) && numeric > 0 ? numeric : null;
};

const escapeRegex = (string) => String(string ?? '').replace(/[.*+?^${}()|[\]\\]/g, '\\$&');

const NON_ORGANISM_PATTERN = /^(?:synthetic construct|artificial sequence|unidentified(?: organism)?|unknown(?: organism)?|vector|(?:unidentified )?cloning vector|expression vector)$/i;

export const formatInferredOrganismStrain = ({ organism = '', strain = '' } = {}) => {
  const rawOrganism = String(organism || '').trim();
  const rawStrain = String(strain || '').trim();
  if (!rawOrganism || NON_ORGANISM_PATTERN.test(rawOrganism)) {
    return rawStrain;
  }

  const isCandidatus = /^Candidatus\s+/i.test(rawOrganism);
  const nameWithoutCand = isCandidatus ? rawOrganism.replace(/^Candidatus\s+/i, '').trim() : rawOrganism;
  const words = nameWithoutCand.split(/\s+/).filter(Boolean);

  let speciesPart = '';
  let rest = '';

  if (words.length >= 2) {
    const binomial = `${words[0]} ${words[1]}`;
    speciesPart = `${isCandidatus ? 'Candidatus ' : ''}<i>${binomial}</i>`;
    rest = words.slice(2).join(' ');
  } else if (words.length === 1) {
    speciesPart = `${isCandidatus ? 'Candidatus ' : ''}<i>${words[0]}</i>`;
  } else {
    speciesPart = rawOrganism;
  }

  if (rawStrain && !rest.toLowerCase().includes(rawStrain.toLowerCase())) {
    rest = rest ? `${rest} ${rawStrain}` : rawStrain;
  }

  return rest ? `${speciesPart} ${rest}`.trim() : speciesPart.trim();
};

export const formatInferredSubtitle = ({ definition = '', plasmid = '', chromosome = '', organism = '' } = {}) => {
  const rawPlasmid = String(plasmid || '').trim();
  if (rawPlasmid) {
    const clean = rawPlasmid.replace(/^plasmid\s+/i, '').trim();
    return `Plasmid ${clean}`;
  }
  const rawChromosome = String(chromosome || '').trim();
  if (rawChromosome) {
    const clean = rawChromosome.replace(/^chromosome\s+/i, '').trim();
    return `Chromosome ${clean}`;
  }
  const rawDef = String(definition || '').replace(/\s+/g, ' ').trim().replace(/\.$/, '');
  if (!rawDef) return '';

  if (/plasmid\b/i.test(rawDef)) {
    const plasmidMatch = rawDef.match(/(?:plasmid\s+([A-Za-z0-9_-]+)|(p[A-Za-z0-9_-]+)\b)/i);
    if (plasmidMatch) {
      const pName = (plasmidMatch[1] || plasmidMatch[2] || '').trim();
      return pName.toLowerCase().startsWith('plasmid') ? pName : `Plasmid ${pName}`;
    }
  }

  if (/complete\s+genome/i.test(rawDef)) {
    if (/mitochondri/i.test(rawDef)) return 'Mitochondrion, complete genome';
    if (/chloroplast/i.test(rawDef)) return 'Chloroplast, complete genome';
    return 'Complete genome';
  }
  if (/complete\s+sequence/i.test(rawDef)) {
    return 'Complete sequence';
  }

  if (organism && !NON_ORGANISM_PATTERN.test(organism)) {
    const orgPrefix = new RegExp(`^${escapeRegex(organism)}[,\\s]*`, 'i');
    const stripped = rawDef.replace(orgPrefix, '').replace(/^(?:DNA|genomic DNA|cDNA)[,\\s]*/i, '').trim();
    if (/(?:gene cluster|biosynthetic gene cluster|cluster|operon)/i.test(stripped)) {
      return stripped.charAt(0).toUpperCase() + stripped.slice(1);
    }
  }

  const clusterMatch = rawDef.match(/([A-Za-z0-9_-]+(?:\s+[A-Za-z0-9_-]+)*\s+(?:gene cluster|biosynthetic gene cluster|cluster|operon))/i);
  if (clusterMatch) {
    const res = clusterMatch[1].trim();
    return res.charAt(0).toUpperCase() + res.slice(1);
  }

  return '';
};

export const extractGenBankMetadata = (chunk) => {
  const text = String(chunk || '');
  const defMatch = text.match(/^DEFINITION\s+([\s\S]*?)(?=^[A-Z]|\/\/)/m);
  const definition = defMatch ? defMatch[1].replace(/\r?\n\s+/g, ' ').trim() : '';

  const sourceMatch = text.match(/^ {5}source\s+[\s\S]*?(?=^ {5}[a-z]|\/\/)/mi);
  const sourceBlock = sourceMatch ? sourceMatch[0] : text;

  const extractQualifier = (key, block) => {
    const regex = new RegExp(`/${key}="([^"]*(?:\\r?\\n {21}[^"]*)*)"`, 'i');
    const match = block.match(regex);
    if (!match) return '';
    return match[1].replace(/\r?\n\s+/g, ' ').trim();
  };

  let organism = extractQualifier('organism', sourceBlock);
  if (!organism) {
    const orgLine = text.match(/^ {2}ORGANISM\s+([^\r\n]+)/m);
    if (orgLine) organism = orgLine[1].trim();
  }

  const strain = extractQualifier('strain', sourceBlock);
  const isolate = extractQualifier('isolate', sourceBlock);
  const plasmid = extractQualifier('plasmid', sourceBlock);
  const chromosome = extractQualifier('chromosome', sourceBlock);

  const effectiveStrain = strain || isolate || '';
  const inferredDefinition = formatInferredOrganismStrain({
    organism,
    strain: effectiveStrain
  });
  const inferredSubtitle = formatInferredSubtitle({
    definition,
    plasmid,
    chromosome,
    organism
  });

  return {
    organism,
    strain: effectiveStrain,
    plasmid,
    chromosome,
    definition,
    inferredDefinition,
    inferredSubtitle
  };
};

export const normalizeSequenceRecords = (payload) => {
  if (payload?.error) throw new Error(String(payload.error));
  if (!Array.isArray(payload?.records)) throw new Error('Record list response is invalid.');

  const records = [];
  const seenSelectors = new Set();
  payload.records.forEach((entry, index) => {
    const selector = String(entry?.selector ?? `#${index + 1}`).trim();
    if (!selector || seenSelectors.has(selector)) return;
    seenSelectors.add(selector);
    records.push({
      selector,
      recordId: String(entry?.record_id ?? entry?.recordId ?? '').trim() || `Record_${index + 1}`,
      recordLength: normalizeRecordLength(entry?.record_length ?? entry?.recordLength),
      detectedTopology: ['circular', 'linear'].includes(entry?.topology ?? entry?.detectedTopology)
        ? (entry.topology ?? entry.detectedTopology)
        : 'unknown',
      organism: String(entry?.organism ?? '').trim(),
      strain: String(entry?.strain ?? '').trim(),
      inferredDefinition: String(entry?.inferredDefinition ?? entry?.inferred_definition ?? '').trim(),
      inferredSubtitle: String(entry?.inferredSubtitle ?? entry?.inferred_subtitle ?? '').trim()
    });
  });

  if (records.length === 0) throw new Error('No records found.');
  return records;
};

const parseGenBankRecordText = (text) => {
  const records = String(text || '')
    .split(/^\/\/\s*$/m)
    .map((chunk, index) => {
      const locus = chunk.match(/^LOCUS\s+(\S+)(?:\s+(\d+)\s+(?:bp|aa)\b)?/m);
      if (!locus) return null;
      const accession = chunk.match(/^ACCESSION\s+(\S+)/m)?.[1];
      const version = chunk.match(/^VERSION\s+(\S+)/m)?.[1];
      const metadata = extractGenBankMetadata(chunk);
      return {
        selector: `#${index + 1}`,
        record_id: version || accession || locus[1],
        record_length: locus[2] ? Number(locus[2]) : null,
        topology: chunk.match(/^LOCUS\s+\S+\s+\d+\s+(?:bp|aa)\b[^\r\n]*\s(circular|linear)(?:\s|$)/m)?.[1] || 'unknown',
        organism: metadata.organism,
        strain: metadata.strain,
        inferredDefinition: metadata.inferredDefinition,
        inferredSubtitle: metadata.inferredSubtitle
      };
    })
    .filter(Boolean)
    .map((record, index) => ({ ...record, selector: `#${index + 1}` }));
  return normalizeSequenceRecords({ records });
};

const parseFastaRecordText = (text) => {
  const records = [];
  let current = null;
  String(text || '').split(/\r?\n/).forEach((line) => {
    if (line.startsWith('>')) {
      if (current) records.push(current);
      current = {
        selector: `#${records.length + 1}`,
        record_id: line.slice(1).trim().split(/\s+/)[0],
        record_length: 0
      };
    } else if (current) {
      current.record_length += line.replace(/\s+/g, '').length;
    }
  });
  if (current) records.push(current);
  return normalizeSequenceRecords({ records });
};

export const parseSequenceRecordText = (text, format) => {
  if (format === 'genbank') return parseGenBankRecordText(text);
  if (format === 'fasta') return parseFastaRecordText(text);
  throw new Error(`Unsupported format: ${String(format)}.`);
};

const defaultTextReader = (file) => (
  typeof file?.text === 'function' || getSessionResourceSource(file)
    ? () => readFileText(file)
    : null
);

export const discoverSequenceRecords = async ({
  file,
  format,
  readText = null,
  runHelperOperation = runDiagramHelperOperation
}) => {
  if (!file) throw new Error('A sequence file is required.');
  const readSourceText = typeof readText === 'function'
    ? () => readText(file)
    : defaultTextReader(file);
  if (readSourceText) {
    try {
      return parseSequenceRecordText(await readSourceText(), format);
    } catch {
      // The packaged Worker parser handles variants beyond the lightweight text fast path.
    }
  }
  const response = await runHelperOperation(
    DIAGRAM_HELPER_OPERATIONS.LIST_SEQUENCE_RECORDS,
    {
      format,
      files: [{ role: 'source', bytes: await cloneFileBytesForTransfer(file) }]
    }
  );
  return normalizeSequenceRecords(response.result);
};

export const discoverGffFastaRecords = async ({
  gffFile,
  fastaFile,
  readText = null,
  runHelperOperation = runDiagramHelperOperation
}) => {
  if (!gffFile || !fastaFile) throw new Error('GFF3 and FASTA files are required.');
  const readSourceText = typeof readText === 'function'
    ? () => readText(fastaFile)
    : defaultTextReader(fastaFile);
  if (readSourceText) {
    try {
      return parseSequenceRecordText(await readSourceText(), 'fasta');
    } catch {
      // Let the Worker validate the paired GFF3/FASTA record set together.
    }
  }
  const response = await runHelperOperation(
    DIAGRAM_HELPER_OPERATIONS.LIST_GFF_FASTA_RECORDS,
    {
      files: [
        { role: 'gff', bytes: await cloneFileBytesForTransfer(gffFile) },
        { role: 'fasta', bytes: await cloneFileBytesForTransfer(fastaFile) }
      ]
    }
  );
  return normalizeSequenceRecords(response.result);
};
