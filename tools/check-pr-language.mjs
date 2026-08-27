#!/usr/bin/env node

import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { pathToFileURL } from 'node:url';

const OPAQUE_LANGUAGE = Object.freeze([
  { label: 'steady-state', pattern: /\bsteady[ -]state\b/i },
  { label: 'CI topology', pattern: /\bCI[ -]+topolog(?:y|ies)\b/i },
  { label: 'main admission', pattern: /\bmain[ -]+admission\b/i },
  { label: 'release admission', pattern: /\brelease[ -]+admission\b/i },
  { label: 'admission-policy', pattern: /\badmission[ -]+policy\b/i },
  { label: 'transition-only', pattern: /\btransition[ -]+only\b/i },
  { label: 'cutover', pattern: /\bcut[ -]?over\b/i }
]);

const SUMMARY_HEADING = /^[ \t]{0,3}##[ \t]+Plain-language summary[ \t]*$/i;
const MARKDOWN_HEADING = /^[ \t]{0,3}#{1,6}[ \t]+\S/;

const fenceAt = (line) => {
  const match = line.match(/^[ \t]{0,3}(`{3,}|~{3,})/);
  return match ? { character: match[1][0], length: match[1].length } : null;
};

const closesFence = (line, fence) => {
  const trimmed = line.trim();
  return trimmed.length >= fence.length
    && [...trimmed].every((character) => character === fence.character);
};

export const extractPlainLanguageSummary = (body) => {
  if (typeof body !== 'string') return null;

  const lines = body.split(/\r?\n/);
  const summary = [];
  let fence = null;
  let found = false;

  for (const line of lines) {
    if (fence) {
      if (found) summary.push(line);
      if (closesFence(line, fence)) fence = null;
      continue;
    }

    const openingFence = fenceAt(line);
    if (openingFence) {
      fence = openingFence;
      if (found) summary.push(line);
      continue;
    }

    if (!found) {
      if (SUMMARY_HEADING.test(line)) found = true;
      continue;
    }

    if (MARKDOWN_HEADING.test(line)) break;
    summary.push(line);
  }

  return found ? summary.join('\n').trim() : null;
};

const stripFencedCode = (text) => {
  const kept = [];
  let fence = null;

  for (const line of text.split(/\r?\n/)) {
    if (fence) {
      if (closesFence(line, fence)) fence = null;
      kept.push('');
      continue;
    }
    const openingFence = fenceAt(line);
    if (openingFence) {
      fence = openingFence;
      kept.push('');
      continue;
    }
    kept.push(line);
  }

  return kept.join('\n');
};

const stripInlineCode = (text) => {
  let result = '';
  let index = 0;

  while (index < text.length) {
    if (text[index] !== '`') {
      result += text[index];
      index += 1;
      continue;
    }

    let length = 1;
    while (text[index + length] === '`') length += 1;
    const delimiter = '`'.repeat(length);
    const end = text.indexOf(delimiter, index + length);
    if (end === -1) {
      result += delimiter;
      index += length;
      continue;
    }
    result += ' ';
    index = end + length;
  }

  return result;
};

const explanatoryText = (text) => stripInlineCode(stripFencedCode(text))
  .replace(/<!--[\s\S]*?-->/g, ' ');

export const findOpaqueLanguage = (text, { ignoreCode = true } = {}) => {
  const source = ignoreCode ? explanatoryText(text) : text;
  return OPAQUE_LANGUAGE.flatMap(({ label, pattern }) => {
    const match = source.match(pattern);
    return match ? [{ label, phrase: match[0] }] : [];
  });
};

const isPlaceholder = (text) => {
  const visible = explanatoryText(text)
    .replace(/\[([^\]]+)\]\([^)]*\)/g, '$1')
    .replace(/[\s*_>#\[\]{}()!.-]+/g, ' ')
    .trim();
  if (!visible) return true;
  return /^(?:todo|tbd|n\/?a|none|placeholder|coming soon)$/i.test(visible)
    || /^(?:add|describe|write|fill(?: in)?|insert|replace)\b.*\b(?:summary|change|details|text)\b/i.test(visible);
};

const opaqueErrors = (scope, text, options) => findOpaqueLanguage(text, options).map(({ phrase }) => (
  `${scope} uses the opaque phrase "${phrase}". Replace it with the concrete action, object, or post-merge effect.`
));

export const validateTitle = (title) => {
  if (typeof title !== 'string' || !title.trim()) {
    return ['PR title is missing. Provide a concrete action and object with --title.'];
  }
  if (isPlaceholder(title)) {
    return ['PR title is placeholder text. Replace it with a concrete action and object.'];
  }
  return opaqueErrors('PR title', title, { ignoreCode: false });
};

export const validateBody = (body) => {
  const summary = extractPlainLanguageSummary(body);
  if (summary === null) {
    return ['PR body is missing "## Plain-language summary". Add it first and explain what changes, why, and what differs after merge.'];
  }
  if (isPlaceholder(summary)) {
    return ['Plain-language summary is empty or placeholder text. Write 2-4 concrete sentences covering what changes, why, and the post-merge result.'];
  }
  return opaqueErrors('Plain-language summary', summary);
};

export const validateProposal = ({ title, body }) => [
  ...validateTitle(title),
  ...validateBody(body)
];

export const splitShellCommand = (command) => {
  const segments = [];
  let current = [];
  let value = '';
  let started = false;
  let dynamic = false;
  let quote = null;
  let parseError = '';

  const finishToken = () => {
    if (!started) return;
    current.push({ value, dynamic });
    value = '';
    started = false;
    dynamic = false;
  };
  const finishSegment = () => {
    finishToken();
    if (current.length) segments.push(current);
    current = [];
  };

  for (let index = 0; index < command.length; index += 1) {
    const character = command[index];
    const next = command[index + 1];

    if (quote === "'") {
      if (character === "'") quote = null;
      else value += character;
      started = true;
      continue;
    }

    if (quote === '"') {
      if (character === '"') {
        quote = null;
      } else if (character === '\\' && next !== undefined) {
        if (next !== '\n') value += next;
        index += 1;
      } else {
        if (character === '$' || character === '`') dynamic = true;
        value += character;
      }
      started = true;
      continue;
    }

    if (character === '\\' && next !== undefined) {
      if (next !== '\n') {
        value += next;
        started = true;
      }
      index += 1;
      continue;
    }
    if (character === "'" || character === '"') {
      quote = character;
      started = true;
      continue;
    }
    if (character === '#' && !started) {
      while (index + 1 < command.length && command[index + 1] !== '\n') index += 1;
      continue;
    }
    if (character === '\n') {
      finishSegment();
      continue;
    }
    if (/\s/.test(character)) {
      finishToken();
      continue;
    }
    if (';|&()'.includes(character)) {
      finishSegment();
      if ((character === '&' || character === '|') && next === character) index += 1;
      continue;
    }
    if (character === '$' || character === '`'
        || ((character === '<' || character === '>') && next === '(')) {
      dynamic = true;
    }
    value += character;
    started = true;
  }

  finishSegment();
  if (quote) parseError = `unclosed ${quote} quote`;
  return { segments, parseError };
};

const isAssignment = ({ value }) => /^[A-Za-z_][A-Za-z0-9_]*=/.test(value);

const findPrCommands = (segments) => segments.flatMap((tokens, segmentIndex) => {
  let commandIndex = 0;
  while (tokens[commandIndex] && isAssignment(tokens[commandIndex])) commandIndex += 1;
  if (tokens[commandIndex]?.value !== 'gh'
      || tokens[commandIndex + 1]?.value !== 'pr'
      || !['create', 'edit'].includes(tokens[commandIndex + 2]?.value)) {
    return [];
  }
  return [{
    action: tokens[commandIndex + 2].value,
    args: tokens.slice(commandIndex + 3),
    segmentIndex
  }];
});

const optionOccurrences = (args, name) => {
  const occurrences = [];
  for (let index = 0; index < args.length; index += 1) {
    const token = args[index];
    if (token.value === name) {
      const supplied = args[index + 1];
      occurrences.push(supplied && !supplied.value.startsWith('--')
        ? supplied
        : { value: '', dynamic: false, missing: true });
      index += supplied && !supplied.value.startsWith('--') ? 1 : 0;
    } else if (token.value.startsWith(`${name}=`)) {
      occurrences.push({
        value: token.value.slice(name.length + 1),
        dynamic: token.dynamic,
        missing: token.value.length === name.length + 1
      });
    }
  }
  return occurrences;
};

const oneInspectableOption = (args, name, label, { path = false } = {}) => {
  const occurrences = optionOccurrences(args, name);
  if (occurrences.length !== 1 || occurrences[0].missing || !occurrences[0].value) {
    return { error: `Provide exactly one explicit ${name} value so the hook can inspect the ${label}.` };
  }
  const token = occurrences[0];
  if (token.dynamic || (path && /[*?{}\[]/.test(token.value))) {
    return { error: `${name} uses an unresolved shell expansion. Provide a literal, inspectable ${label}.` };
  }
  return { value: token.value };
};

const bodyFromFile = (args, readFile) => {
  const bodyFile = oneInspectableOption(args, '--body-file', 'body file path', { path: true });
  if (bodyFile.error) return bodyFile;
  if (bodyFile.value === '-') {
    return { error: '--body-file - cannot be inspected. Save the complete PR body to a readable file first.' };
  }
  try {
    return { value: readFile(bodyFile.value, 'utf8') };
  } catch (_error) {
    return { error: `Cannot read PR body file "${bodyFile.value}". Create a readable file first, then run the PR command separately.` };
  }
};

const denied = (reason) => ({ allowed: false, reason });
const allowed = () => ({ allowed: true });

const validationDecision = (errors) => errors.length ? denied(errors.join(' ')) : allowed();

const evaluateCreate = (args, readFile) => {
  const fillFlag = args.find(({ value }) => (
    ['--fill', '--fill-first', '--fill-verbose'].includes(value.split('=')[0])
  ));
  if (fillFlag) {
    return denied(`${fillFlag.value.split('=')[0]} hides the proposed wording. Supply an explicit --title and a readable --body-file instead.`);
  }

  const title = oneInspectableOption(args, '--title', 'PR title');
  if (title.error) return denied(title.error);
  const body = bodyFromFile(args, readFile);
  if (body.error) return denied(body.error);
  return validationDecision(validateProposal({ title: title.value, body: body.value }));
};

const evaluateEdit = (args, readFile) => {
  const titleOptions = optionOccurrences(args, '--title');
  const bodyOptions = optionOccurrences(args, '--body');
  const bodyFileOptions = optionOccurrences(args, '--body-file');
  if (!titleOptions.length && !bodyOptions.length && !bodyFileOptions.length) return allowed();

  const errors = [];
  if (titleOptions.length) {
    const title = oneInspectableOption(args, '--title', 'PR title');
    if (title.error) errors.push(title.error);
    else errors.push(...validateTitle(title.value));
  }

  if (bodyOptions.length && bodyFileOptions.length) {
    errors.push('Use either --body or --body-file for a wording change, not both.');
  } else if (bodyOptions.length) {
    const body = oneInspectableOption(args, '--body', 'PR body');
    if (body.error) errors.push(body.error);
    else errors.push(...validateBody(body.value));
  } else if (bodyFileOptions.length) {
    const body = bodyFromFile(args, readFile);
    if (body.error) errors.push(body.error);
    else errors.push(...validateBody(body.value));
  }

  return validationDecision(errors);
};

export const evaluateHookCommand = (command, { readFile = readFileSync } = {}) => {
  if (typeof command !== 'string') return allowed();
  const parsed = splitShellCommand(command);
  const prCommands = findPrCommands(parsed.segments);
  if (!prCommands.length) return allowed();
  if (parsed.parseError) {
    return denied(`Cannot inspect this PR command because it has an ${parsed.parseError}. Use ordinary explicit gh pr arguments.`);
  }

  for (const prCommand of prCommands) {
    const changesWording = prCommand.action === 'create'
      || ['--title', '--body', '--body-file'].some((name) => (
        optionOccurrences(prCommand.args, name).length > 0
      ));
    if (changesWording && parsed.segments.length !== 1) {
      return denied(`Run gh pr ${prCommand.action} as a separate command after the complete body file exists so the hook can inspect the final wording.`);
    }
    const decision = prCommand.action === 'create'
      ? evaluateCreate(prCommand.args, readFile)
      : evaluateEdit(prCommand.args, readFile);
    if (!decision.allowed) return decision;
  }
  return allowed();
};

export const hookDenialOutput = (reason) => ({
  hookSpecificOutput: {
    hookEventName: 'PreToolUse',
    permissionDecision: 'deny',
    permissionDecisionReason: reason
  }
});

export const evaluateHookEvent = (event, options) => evaluateHookCommand(
  event?.tool_input?.command,
  options
);

const readManualArguments = (args) => {
  const values = new Map();
  for (let index = 0; index < args.length; index += 1) {
    const name = args[index];
    if (!['--title', '--body-file'].includes(name) || index + 1 >= args.length) {
      return { error: 'Usage: node tools/check-pr-language.mjs --title "..." --body-file path/to/body.md' };
    }
    if (values.has(name)) return { error: `${name} may be supplied only once.` };
    values.set(name, args[index + 1]);
    index += 1;
  }
  if (!values.has('--title') || !values.has('--body-file')) {
    return { error: 'Manual mode requires both --title and --body-file.' };
  }
  return { title: values.get('--title'), bodyFile: values.get('--body-file') };
};

const runManual = (args) => {
  const parsed = readManualArguments(args);
  if (parsed.error) {
    process.stderr.write(`${parsed.error}\n`);
    return 1;
  }

  let body;
  try {
    body = readFileSync(parsed.bodyFile, 'utf8');
  } catch (_error) {
    process.stderr.write(`Cannot read PR body file "${parsed.bodyFile}".\n`);
    return 1;
  }

  const errors = validateProposal({ title: parsed.title, body });
  if (errors.length) {
    process.stderr.write(`PR language check failed:\n- ${errors.join('\n- ')}\n`);
    return 1;
  }
  process.stdout.write('PR language check passed.\n');
  return 0;
};

const runHook = () => {
  let event;
  try {
    event = JSON.parse(readFileSync(0, 'utf8'));
  } catch (_error) {
    process.stdout.write(`${JSON.stringify(hookDenialOutput(
      'Cannot inspect the Codex hook event. Retry the PR command with ordinary explicit arguments.'
    ))}\n`);
    return 0;
  }
  const decision = evaluateHookEvent(event);
  if (!decision.allowed) process.stdout.write(`${JSON.stringify(hookDenialOutput(decision.reason))}\n`);
  return 0;
};

const run = (args) => args.length === 1 && args[0] === '--hook'
  ? runHook()
  : runManual(args);

if (process.argv[1]
    && import.meta.url === pathToFileURL(resolve(process.argv[1])).href) {
  process.exitCode = run(process.argv.slice(2));
}
