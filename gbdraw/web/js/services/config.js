import { state } from '../state.js';
import { resolveColorToHex } from '../app/color-utils.js';

const safeDeepMerge = (target, source) => {
  if (!source || typeof source !== 'object') return;

  Object.keys(source).forEach((key) => {
    // 1. Prevent prototype pollution
    if (['__proto__', 'constructor', 'prototype'].includes(key)) return;

    // 2. Ignore keys not present in target (whitelisting effect)
    if (!Object.prototype.hasOwnProperty.call(target, key)) return;

    const targetValue = target[key];
    const sourceValue = source[key];

    // 3. Recursive merge for objects
    if (
      targetValue &&
      typeof targetValue === 'object' &&
      !Array.isArray(targetValue) &&
      sourceValue &&
      typeof sourceValue === 'object' &&
      !Array.isArray(sourceValue)
    ) {
      safeDeepMerge(targetValue, sourceValue);
      return;
    }

    // 4. For arrays, intentionally overwrite (replacing lists of settings is natural)
    if (Array.isArray(targetValue) && Array.isArray(sourceValue)) {
      target[key].splice(0, target[key].length, ...sourceValue);
      return;
    }

    // 5. Update value only if types match or initial value is null
    if (typeof targetValue === typeof sourceValue || targetValue === null) {
      target[key] = sourceValue;
    }
  });
};

export const exportConfig = () => {
  const configData = {
    form: state.form,
    adv: state.adv,
    colors: state.currentColors.value,
    palette: state.selectedPalette.value,
    rules: state.manualSpecificRules,
    filterMode: state.filterMode.value,
    whitelist: state.manualWhitelist,
    blacklistText: state.manualBlacklist.value
  };
  const blob = new Blob([JSON.stringify(configData, null, 2)], {
    type: 'application/json'
  });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = 'gbdraw_config.json';
  a.click();
  URL.revokeObjectURL(url);
};

export const importConfig = (e) => {
  const file = e.target.files[0];
  if (!file) return;

  if (file.size > 10 * 1024 * 1024) {
    alert('Config file is too large.');
    return;
  }

  const reader = new FileReader();
  reader.onload = (event) => {
    try {
      const data = JSON.parse(event.target.result, (key, value) => {
        if (key === '__proto__' || key === 'constructor' || key === 'prototype') {
          return undefined;
        }
        return value;
      });

      if (data.form) safeDeepMerge(state.form, data.form);
      if (data.adv) safeDeepMerge(state.adv, data.adv);
      if (data.colors) {
        const normalized = {};
        Object.entries(data.colors).forEach(([key, value]) => {
          normalized[key] = resolveColorToHex(String(value || '').trim());
        });
        state.currentColors.value = normalized;
      }
      if (data.palette) state.selectedPalette.value = data.palette;

      if (data.rules && Array.isArray(data.rules)) {
        state.manualSpecificRules.length = 0;
        data.rules.forEach((r) => {
          state.manualSpecificRules.push({
            feat: String(r.feat || ''),
            qual: String(r.qual || ''),
            val: String(r.val || ''),
            color: resolveColorToHex(String(r.color || '#000000')),
            cap: String(r.cap || '')
          });
        });
      }
      if (data.filterMode) state.filterMode.value = data.filterMode;
      if (data.whitelist && Array.isArray(data.whitelist)) {
        state.manualWhitelist.length = 0;
        data.whitelist.forEach((w) => {
          state.manualWhitelist.push({
            feat: String(w.feat || ''),
            qual: String(w.qual || ''),
            key: String(w.key || '')
          });
        });
      }
      if (data.blacklistText) state.manualBlacklist.value = String(data.blacklistText);

      alert('Configuration loaded successfully!');
    } catch (err) {
      console.error(err);
      alert('Failed to load config: Invalid JSON structure.');
    }
  };
  reader.readAsText(file);
  e.target.value = '';
};
