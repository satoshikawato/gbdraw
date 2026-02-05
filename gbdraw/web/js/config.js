export const GBDRAW_WHEEL_NAME = "gbdraw-0.9.0b0-py3-none-any.whl";

export const DEBUG = false;

export const debugLog = (...args) => {
  if (DEBUG) {
    console.log('[DEBUG]', ...args);
  }
};
