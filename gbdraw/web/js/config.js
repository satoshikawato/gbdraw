export const GBDRAW_WHEEL_NAME = "gbdraw-0.9.0b0-py3-none-any.whl";
export const GBDRAW_WHEEL_CACHE_BUST = "20260222-1600";

export const DEBUG = false;

export const debugLog = (...args) => {
  if (DEBUG) {
    console.log('[DEBUG]', ...args);
  }
};
