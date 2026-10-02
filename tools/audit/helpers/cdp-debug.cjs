// Chrome DevTools Protocol helpers for investigating a finding in a Playwright page.
// Use them from a throwaway spec; they are not used by the sweeps.

// Record every exception thrown in the main thread while armed, caught or not.
// Caught exceptions are invisible to page.on('pageerror'); this finds swallowed failures.
const captureExceptions = async (page) => {
  const client = await page.context().newCDPSession(page);
  const captured = [];
  let armed = false;
  await client.send('Debugger.enable');
  client.on('Debugger.paused', async (event) => {
    try {
      if (armed && event.reason === 'exception') {
        const frame = event.callFrames?.[0];
        captured.push({
          description: String(event.data?.description || event.data?.value || '').split('\n').slice(0, 3).join(' | '),
          at: frame ? `${frame.url.split('/').slice(-2).join('/')}:${frame.location.lineNumber + 1}` : ''
        });
      }
    } finally {
      await client.send('Debugger.resume').catch(() => {});
    }
  });
  return {
    arm: async () => { armed = true; await client.send('Debugger.setPauseOnExceptions', { state: 'all' }); },
    disarm: async () => { armed = false; await client.send('Debugger.setPauseOnExceptions', { state: 'none' }); },
    captured
  };
};

// Evaluate an expression in the paused frame each time a line is hit, then resume.
// lineNumber is 1-based, as shown in an editor.
const captureAtBreakpoint = async (page, urlRegex, lineNumber, expression) => {
  const client = await page.context().newCDPSession(page);
  const values = [];
  await client.send('Debugger.enable');
  await client.send('Debugger.setBreakpointByUrl', { urlRegex, lineNumber: lineNumber - 1 });
  client.on('Debugger.paused', async (event) => {
    try {
      if (event.reason !== 'exception' && event.reason !== 'promiseRejection') {
        const frame = event.callFrames?.[0];
        const result = await client.send('Debugger.evaluateOnCallFrame', {
          callFrameId: frame.callFrameId, expression, returnByValue: true
        });
        values.push(result.result?.value ?? result.exceptionDetails?.text);
      }
    } catch (error) {
      values.push(String(error));
    } finally {
      await client.send('Debugger.resume').catch(() => {});
    }
  });
  return values;
};

module.exports = { captureAtBreakpoint, captureExceptions };
