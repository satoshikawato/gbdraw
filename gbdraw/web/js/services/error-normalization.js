const DEFAULT_SUMMARY_LIMIT = 1000;
const DEFAULT_DETAIL_LIMIT = 4000;
const OMITTED_DETAIL_KEYS = /^(?:base64|bytes?|content|data|file|files|sequence|sequences|stack|traceback)$/i;

const truncateText = (value, limit) => {
  const text = String(value ?? '').trim();
  if (text.length <= limit) return text;
  const suffix = '\n[details truncated]';
  return `${text.slice(0, Math.max(0, limit - suffix.length)).trimEnd()}${suffix}`;
};

const stripTraceback = (value) => {
  const lines = String(value ?? '').split(/\r?\n/);
  const kept = [];
  let inTraceback = false;
  for (const line of lines) {
    const trimmed = line.trim();
    if (/^Traceback(?:\s|\(|:)/i.test(trimmed)) {
      inTraceback = true;
      continue;
    }
    if (
      inTraceback &&
      (
        /^\s*File\s+"/.test(line) ||
        /^\s*(?:at\s+|[\^~]+\s*$)/.test(line) ||
        /^\s+/.test(line)
      )
    ) {
      continue;
    }
    if (inTraceback && trimmed) {
      inTraceback = false;
      continue;
    }
    if (/^\s*at\s+\S+/.test(line)) continue;
    kept.push(line);
  }
  return kept.join('\n').trim();
};

const primitiveText = (value) => {
  if (value === null || value === undefined) return '';
  if (typeof value === 'string') return value;
  if (typeof value === 'number' || typeof value === 'bigint' || typeof value === 'boolean') {
    return String(value);
  }
  return '';
};

export const safeErrorText = (value, { limit = DEFAULT_DETAIL_LIMIT } = {}) => {
  const primitive = primitiveText(value);
  if (primitive) return truncateText(stripTraceback(primitive), limit);
  if (!value || typeof value !== 'object') {
    return truncateText(stripTraceback(String(value ?? '')), limit);
  }

  if (value instanceof Error) {
    return truncateText(stripTraceback(value.message || value.name || 'Error'), limit);
  }

  const preferred = [
    value.summary,
    value.message,
    value.reason,
    value.error_description,
    value.error
  ];
  for (const candidate of preferred) {
    const text = primitiveText(candidate);
    if (text) return truncateText(stripTraceback(text), limit);
  }

  const fields = [];
  for (const [key, candidate] of Object.entries(value)) {
    if (OMITTED_DETAIL_KEYS.test(key)) continue;
    const text = primitiveText(candidate);
    if (!text) continue;
    fields.push(`${key}: ${stripTraceback(text)}`);
    if (fields.length >= 8) break;
  }
  return truncateText(
    fields.filter(Boolean).join('\n') || 'Unexpected structured error.',
    limit
  );
};

const lastNonemptyLine = (value) => {
  const lines = stripTraceback(value).split(/\r?\n/).map((line) => line.trim()).filter(Boolean);
  return lines[lines.length - 1] || '';
};

const normalizeDetailSections = (sections, limit) => {
  if (!Array.isArray(sections)) return [];
  const normalized = [];
  for (const section of sections) {
    const label = safeErrorText(section?.label || 'Details', { limit: 80 }) || 'Details';
    const text = safeErrorText(section?.text ?? section, { limit });
    if (!text || /^(?:traceback|stack)$/i.test(label)) continue;
    normalized.push({ label, text });
    if (normalized.length >= 8) break;
  }
  return normalized;
};

export const normalizeUserFacingError = (
  value,
  {
    summaryLimit = DEFAULT_SUMMARY_LIMIT,
    detailLimit = DEFAULT_DETAIL_LIMIT
  } = {}
) => {
  if (!value) return null;

  let error = value;
  if (typeof error === 'string') {
    const trimmed = error.trim();
    if (trimmed.startsWith('{') && trimmed.endsWith('}')) {
      try {
        error = JSON.parse(trimmed);
      } catch {
        error = trimmed;
      }
    }
  }

  if (!error || typeof error !== 'object') {
    const text = safeErrorText(error, { limit: detailLimit });
    return {
      summary: truncateText(lastNonemptyLine(text) || 'Unknown error', summaryLimit),
      details: text ? [{ label: 'Details', text }] : []
    };
  }

  const type = safeErrorText(error.type || error.name || error.stage, { limit: 120 });
  const message = safeErrorText(
    error.summary ||
      error.message ||
      lastNonemptyLine(error.stderr) ||
      lastNonemptyLine(error.stdout) ||
      error.reason ||
      safeErrorText(error, { limit: summaryLimit }) ||
      'Unknown error',
    { limit: summaryLimit }
  );
  const summary = (
    type &&
    !/^(?:error|unknown)$/i.test(type) &&
    message &&
    !message.toLowerCase().startsWith(type.toLowerCase())
  )
    ? `${type}: ${message}`
    : message;

  const details = normalizeDetailSections(error.details, detailLimit);
  for (const [label, raw] of [['STDERR', error.stderr], ['STDOUT', error.stdout]]) {
    const text = safeErrorText(raw, { limit: detailLimit });
    if (text) details.push({ label, text });
  }
  const notes = Array.isArray(error.notes) ? error.notes : [];
  for (const note of notes) {
    const text = safeErrorText(note, { limit: detailLimit });
    if (text) details.push({ label: 'Cleanup note', text });
  }

  return {
    summary: truncateText(summary || 'Unknown error', summaryLimit),
    details: details.slice(0, 8)
  };
};
