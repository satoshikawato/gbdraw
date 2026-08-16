export const maskJavaScript = (source = '', { strings = true } = {}) => {
  const masked = source.split('');
  let state = 'code';
  let escaped = false;
  let regexCharacterClass = false;
  const hide = (index) => {
    if (source[index] !== '\n' && source[index] !== '\r') masked[index] = ' ';
  };

  for (let index = 0; index < source.length; index += 1) {
    const character = source[index];
    const next = source[index + 1];

    if (state === 'line-comment') {
      if (character === '\n' || character === '\r') state = 'code';
      else hide(index);
      continue;
    }
    if (state === 'block-comment') {
      hide(index);
      if (character === '*' && next === '/') {
        hide(index + 1);
        index += 1;
        state = 'code';
      }
      continue;
    }
    if (state === 'regex') {
      if (strings) hide(index);
      if (escaped) {
        escaped = false;
        continue;
      }
      if (character === '\\') {
        escaped = true;
        continue;
      }
      if (character === '[') regexCharacterClass = true;
      else if (character === ']') regexCharacterClass = false;
      else if (character === '/' && !regexCharacterClass) state = 'code';
      continue;
    }
    if (state !== 'code') {
      if (strings) hide(index);
      if (escaped) {
        escaped = false;
        continue;
      }
      if (character === '\\') {
        escaped = true;
        continue;
      }
      const delimiter = state === 'single-quote' ? "'" : state === 'double-quote' ? '"' : '`';
      if (character === delimiter) state = 'code';
      continue;
    }

    if (character === '/' && next === '/') {
      hide(index);
      hide(index + 1);
      index += 1;
      state = 'line-comment';
    } else if (character === '/' && next === '*') {
      hide(index);
      hide(index + 1);
      index += 1;
      state = 'block-comment';
    } else if (character === '/') {
      let previousIndex = index - 1;
      while (previousIndex >= 0 && /\s/.test(source[previousIndex])) previousIndex -= 1;
      const previous = source[previousIndex] || '';
      const followsOperator = !previous || /[({\[=,:;!?&|+\-*%^~<>]/.test(previous);
      const prefix = source.slice(Math.max(0, index - 24), index);
      const followsKeyword = /\b(?:return|case|throw|yield|await|else|do|typeof|instanceof|in|of)\s*$/.test(prefix);
      if (followsOperator || followsKeyword) {
        if (strings) hide(index);
        regexCharacterClass = false;
        state = 'regex';
      }
    } else if (character === "'" || character === '"' || character === '`') {
      if (strings) hide(index);
      state = character === "'" ? 'single-quote' : character === '"' ? 'double-quote' : 'template';
    }
  }
  return masked.join('');
};

export const literalImportSpecifiers = (source = '', { dynamic = true } = {}) => {
  const specifiers = [];
  const commentsMasked = maskJavaScript(source, { strings: false });
  const code = maskJavaScript(source);
  const patterns = [
    /(?:^|\n)\s*(?:import|export)\s+[^;]*?\s+from\s+(['"])([^'"]+)\1\s*;?/g,
    /(?:^|\n)\s*import\s*(['"])([^'"]+)\1\s*;?/g
  ];
  if (dynamic) patterns.push(/\bimport\s*\(\s*(['"])([^'"]+)\1\s*\)/g);
  patterns.forEach((pattern) => {
    for (const match of commentsMasked.matchAll(pattern)) {
      const keywordOffset = match[0].search(/\b(?:import|export)\b/);
      const keywordIndex = match.index + keywordOffset;
      if (!/^(?:import|export)\b/.test(code.slice(keywordIndex))) continue;
      specifiers.push(match[2]);
    }
  });
  return specifiers;
};
