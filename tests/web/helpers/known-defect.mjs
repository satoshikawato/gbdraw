import assert from 'node:assert/strict';

// Node counterpart of Playwright test.fail(true, '<ID>') and pytest
// xfail(strict=True) for the Web GUI audit 2026-09-30. `check` asserts the
// correct behavior. While the defect exists the assertion fails and the test
// passes; once the behavior is fixed this helper fails so the fixing PR removes
// the mark. Errors other than assertion failures still fail the test.
export const assertKnownDefect = async (id, check) => {
  assert.match(id, /^(?:[A-Z]{1,3}-\d{2}|N-\d{2})$/, 'a known defect names its audit ID');
  try {
    await check();
  } catch (error) {
    if (error instanceof assert.AssertionError) return error;
    throw error;
  }
  assert.fail(`${id} no longer reproduces; remove its known-defect mark.`);
};
