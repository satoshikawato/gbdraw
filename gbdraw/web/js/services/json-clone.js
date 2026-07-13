export const cloneJsonData = (value) => {
  if (value === null || value === undefined) return value;
  return JSON.parse(JSON.stringify(value));
};

export const cloneJsonValue = (value, fallback) => {
  if (value === undefined) return fallback;
  try {
    return cloneJsonData(value);
  } catch (_err) {
    return fallback;
  }
};
