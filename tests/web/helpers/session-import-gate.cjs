// Delay the actual File handoff to the import Worker; retain the client lifecycle.
const installImportReadGate = async (page, filename) => {
  await page.evaluate(gatedFilename => {
    const nativePost = Worker.prototype.postMessage;
    let release;
    const gate = new Promise(resolve => { release = resolve; });
    window.__GBDRAW_SESSION_IMPORT_GATE__ = { filename: gatedFilename, streamInvocations: 0, release };
    Worker.prototype.postMessage = function(message, ...args) {
      if (message.file?.name !== gatedFilename) return nativePost.call(this, message, ...args);
      window.__GBDRAW_SESSION_IMPORT_GATE__.streamInvocations += 1;
      gate.then(() => nativePost.call(this, message, ...args));
    };
  }, filename);
};
module.exports = { installImportReadGate };
