const SAFE_SHELL_TOKEN_RE = /^[A-Za-z0-9_@%+=:,./-]+$/;
const FORMAT_FLAGS = new Set(['-f', '--format']);

const normalizePath = (value) => String(value ?? '').trim();

const fallbackNameFromPath = (path) => {
  const name = normalizePath(path).split('/').filter(Boolean).pop();
  return name || 'file';
};

const normalizeMetadataEntries = (fileMetadata) => {
  if (!fileMetadata) return [];
  if (fileMetadata instanceof Map) return Array.from(fileMetadata.entries());
  if (Array.isArray(fileMetadata)) {
    return fileMetadata
      .map((entry) => {
        if (Array.isArray(entry) && entry.length >= 2) return [entry[0], entry[1]];
        if (entry && typeof entry === 'object') return [entry.path, entry];
        return null;
      })
      .filter(Boolean);
  }
  if (typeof fileMetadata === 'object') return Object.entries(fileMetadata);
  return [];
};

const normalizeFileMetadata = (fileMetadata) => {
  const map = new Map();
  normalizeMetadataEntries(fileMetadata).forEach(([pathRaw, entryRaw]) => {
    const path = normalizePath(pathRaw || entryRaw?.path);
    if (!path) return;
    const entry = entryRaw && typeof entryRaw === 'object' ? entryRaw : {};
    const name = String(entry.name || fallbackNameFromPath(path)).trim() || fallbackNameFromPath(path);
    const kind = String(entry.kind || '').trim() === 'generated' ? 'generated' : 'uploaded';
    const slot = String(entry.slot || '').trim();
    map.set(path, {
      path,
      name,
      kind,
      slot
    });
  });
  return map;
};

const hasRenderFormatArg = (args) => {
  for (let idx = 0; idx < args.length; idx += 1) {
    const token = String(args[idx] ?? '');
    if (FORMAT_FLAGS.has(token)) return true;
    if (token.startsWith('--format=')) return true;
  }
  return false;
};

const withSvgFormat = (args) => (hasRenderFormatArg(args) ? [...args] : [...args, '-f', 'svg']);

export const quoteShellArg = (value) => {
  const token = String(value ?? '');
  if (token.length > 0 && SAFE_SHELL_TOKEN_RE.test(token)) return token;
  return `'${token.replace(/'/g, `'\\''`)}'`;
};

export const buildShellCommand = (tokens) => (Array.isArray(tokens) ? tokens : [])
  .map((token) => quoteShellArg(token))
  .join(' ');

export const formatElapsedMs = (elapsedMs) => {
  const ms = Number(elapsedMs);
  if (!Number.isFinite(ms) || ms < 0) return 'n/a';
  if (ms < 1000) return `${Math.round(ms)} ms`;
  if (ms < 60_000) return `${(ms / 1000).toFixed(ms < 10_000 ? 2 : 1)} s`;
  const minutes = Math.floor(ms / 60_000);
  const seconds = (ms - minutes * 60_000) / 1000;
  return `${minutes} min ${seconds.toFixed(1)} s`;
};

export const reproducibilityLabel = (level) => {
  const normalized = String(level || '').trim();
  if (normalized === 'exact-uploaded-files') return 'Uploaded files';
  if (normalized === 'requires-helper-files') return 'Helper files needed';
  if (normalized === 'session-recommended') return 'Session recommended';
  return 'Pseudo command';
};

export const isCliInvocationSessionExportable = (invocation) => {
  if (!invocation || typeof invocation !== 'object') return false;
  if (invocation.sessionExportable === false) return false;
  const bindings = Array.isArray(invocation.fileBindings) ? invocation.fileBindings : [];
  return bindings.every((binding) => String(binding?.slot || '').startsWith('files.'));
};

export const buildRunInfo = ({
  mode,
  args,
  fileMetadata,
  elapsedMs,
  resultCount,
  startedAtIso,
  generatedBy = 'gbdraw-web'
} = {}) => {
  const normalizedMode = String(mode || '').trim() === 'linear' ? 'linear' : 'circular';
  const metadata = normalizeFileMetadata(fileMetadata);
  const helperFiles = [];
  const unresolvedFileArgs = [];
  const displayArgs = (Array.isArray(args) ? args : []).map((arg, argIndex) => {
    const token = String(arg ?? '');
    const meta = metadata.get(token);
    if (!meta) {
      if (token.startsWith('/')) {
        unresolvedFileArgs.push({ argIndex, path: token });
      }
      return token;
    }
    const displayName = meta.name || fallbackNameFromPath(token);
    if (meta.kind === 'generated') {
      helperFiles.push({
        path: meta.path,
        name: displayName,
        slot: meta.slot || ''
      });
    }
    return displayName;
  });

  const argsWithFormat = withSvgFormat(displayArgs);
  const bindingArgIndexesByName = new Map();
  argsWithFormat.forEach((token, index) => {
    if (!bindingArgIndexesByName.has(token)) bindingArgIndexesByName.set(token, []);
    bindingArgIndexesByName.get(token).push(index);
  });

  const fixedFileBindings = [];
  const seenBindingKeys = new Set();
  (Array.isArray(args) ? args : []).forEach((arg, originalIndex) => {
    const meta = metadata.get(String(arg ?? ''));
    if (!meta?.slot) return;
    const displayName = meta.name || fallbackNameFromPath(meta.path);
    const indexes = bindingArgIndexesByName.get(displayName) || [];
    const argIndex = indexes.shift();
    if (!Number.isInteger(argIndex)) return;
    const key = `${argIndex}:${meta.slot}`;
    if (seenBindingKeys.has(key)) return;
    seenBindingKeys.add(key);
    fixedFileBindings.push({
      argIndex,
      slot: meta.slot,
      name: displayName,
      originalArgIndex: originalIndex
    });
  });

  const invocation = {
    schema: 1,
    mode: normalizedMode,
    args: argsWithFormat,
    renderFormats: ['svg'],
    fileBindings: fixedFileBindings.map(({ originalArgIndex: _originalArgIndex, ...binding }) => binding),
    generatedBy,
    sessionExportable: fixedFileBindings.every((binding) => String(binding.slot).startsWith('files.')) &&
      unresolvedFileArgs.length === 0
  };

  const hasGeneratedBindings = fixedFileBindings.some((binding) => !String(binding.slot).startsWith('files.'));
  const hasLosatHelpers = helperFiles.some((helper) => /losat|blast/i.test(`${helper.slot} ${helper.name}`));
  const notes = [];
  let level = 'exact-uploaded-files';
  if (helperFiles.length > 0) {
    level = hasLosatHelpers ? 'session-recommended' : 'requires-helper-files';
    notes.push('The displayed command assumes browser-generated helper files have been saved next to the input files.');
  }
  if (hasGeneratedBindings) {
    notes.push('Generated helper files are not embedded in exported sessions in this first pass.');
  }
  if (unresolvedFileArgs.length > 0) {
    level = 'pseudo';
    notes.push('Some browser virtual paths could not be mapped to uploaded or generated file names.');
  } else if (fixedFileBindings.length === 0 && metadata.size > 0) {
    level = 'pseudo';
  }
  if (notes.length === 0) {
    notes.push('The command uses the uploaded file names shown here.');
  }

  const commandArgs = ['gbdraw', normalizedMode, ...argsWithFormat];
  const sessionCommand = isCliInvocationSessionExportable(invocation)
    ? buildShellCommand(['gbdraw', normalizedMode, '--session', 'session.gbdraw-session.json', '-f', 'svg'])
    : '';

  return {
    schema: 1,
    mode: normalizedMode,
    startedAtIso: startedAtIso || new Date().toISOString(),
    elapsedMs: Number.isFinite(Number(elapsedMs)) ? Number(elapsedMs) : 0,
    resultCount: Number.isFinite(Number(resultCount)) ? Number(resultCount) : 0,
    command: buildShellCommand(commandArgs),
    commandArgs,
    sessionCommand,
    invocation,
    helperFiles,
    reproducibility: {
      level,
      notes
    }
  };
};
