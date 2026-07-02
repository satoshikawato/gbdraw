const TEXT_INPUT_TYPES = new Set([
  '',
  'date',
  'datetime-local',
  'email',
  'month',
  'number',
  'password',
  'search',
  'tel',
  'text',
  'time',
  'url',
  'week'
]);

export const isIgnoredTarget = (target) =>
  Boolean(
    target?.closest?.('[data-history-ignore], [data-history-managed], [data-history-scope="transient"]')
  );

const isEditableControl = (element) => {
  if (!element) return false;
  const tag = String(element.tagName || '').toLowerCase();
  if (tag === 'textarea') return true;
  if (tag !== 'input') return false;
  return TEXT_INPUT_TYPES.has(String(element.type || '').toLowerCase());
};

const controlLabel = (element) => {
  if (!element) return 'Edit';
  const tag = String(element.tagName || '').toLowerCase();
  if (tag === 'button') return 'Change setting';
  if (tag === 'select') return 'Change setting';
  const type = String(element.type || '').toLowerCase();
  if (type === 'file') return 'Change uploaded file';
  if (type === 'checkbox' || type === 'radio') return 'Change setting';
  if (type === 'color') return 'Change color';
  return 'Edit setting';
};

export const setupHistoryInputs = ({ root, history, nextTick }) => {
  const appRoot = root || document.getElementById('app');
  if (!appRoot || !history) return () => {};

  const txByElement = new WeakMap();
  let activeButtonTx = null;

  const beginForElement = async (element, source = 'input-adapter') => {
    if (!element || element.disabled || isIgnoredTarget(element)) return null;
    const existing = txByElement.get(element);
    if (existing && !existing.closed) return existing;
    const tx = await history.begin(controlLabel(element), { source });
    if (tx) txByElement.set(element, tx);
    return tx;
  };

  const commitElement = async (element) => {
    const tx = txByElement.get(element);
    if (!tx) return;
    if (tx.deferAdapterCommit) return;
    txByElement.delete(element);
    if (tx.closed) return;
    await nextTick();
    await history.commit(tx);
  };

  const findControl = (eventTarget) =>
    eventTarget?.closest?.('input, textarea, select, button, [contenteditable="true"], .upload-zone') || null;

  const onPointerDown = (event) => {
    const target = findControl(event.target);
    if (!target || isIgnoredTarget(target)) return;
    const button = target.closest?.('button');
    if (button) {
      void beginForElement(button).then((tx) => {
        activeButtonTx = tx;
      });
      return;
    }
    if (target.classList?.contains('upload-zone')) {
      const input = target.querySelector?.('input[type="file"]');
      if (input) void beginForElement(input);
      return;
    }
    const tag = String(target.tagName || '').toLowerCase();
    const type = String(target.type || '').toLowerCase();
    if (tag === 'select' || type === 'checkbox' || type === 'radio' || type === 'color' || type === 'file') {
      void beginForElement(target);
    }
  };

  const onFocusIn = (event) => {
    const target = findControl(event.target);
    if (!target || isIgnoredTarget(target)) return;
    if (isEditableControl(target) || target?.isContentEditable) {
      void beginForElement(target);
    }
  };

  const onKeyDown = (event) => {
    const target = findControl(event.target);
    if (!target || isIgnoredTarget(target)) return;
    if (isEditableControl(target) || target?.isContentEditable) {
      void beginForElement(target);
    }
  };

  const onChange = (event) => {
    const target = findControl(event.target);
    if (!target || isIgnoredTarget(target)) return;
    if (!txByElement.has(target)) {
      void beginForElement(target).then(() => commitElement(target));
      return;
    }
    const tag = String(target.tagName || '').toLowerCase();
    if (tag === 'select' || tag === 'input' || tag === 'textarea') {
      void commitElement(target);
    }
  };

  const onFocusOut = (event) => {
    const target = findControl(event.target);
    if (!target || isIgnoredTarget(target)) return;
    if (isEditableControl(target) || target?.isContentEditable) {
      void commitElement(target);
    }
  };

  const onClick = (event) => {
    const button = event.target?.closest?.('button');
    if (!button || isIgnoredTarget(button)) return;
    setTimeout(() => {
      const tx = activeButtonTx || txByElement.get(button);
      activeButtonTx = null;
      if (!tx) return;
      if (tx.deferAdapterCommit) return;
      txByElement.delete(button);
      if (tx.closed) return;
      void nextTick().then(() => history.commit(tx));
    }, 0);
  };

  const onClickCapture = (event) => {
    const button = event.target?.closest?.('button');
    if (!button || isIgnoredTarget(button)) return;
    if (txByElement.has(button) || activeButtonTx) return;
    void beginForElement(button).then((tx) => {
      activeButtonTx = tx;
    });
  };

  appRoot.addEventListener('pointerdown', onPointerDown, true);
  appRoot.addEventListener('focusin', onFocusIn, true);
  appRoot.addEventListener('keydown', onKeyDown, true);
  appRoot.addEventListener('change', onChange, false);
  appRoot.addEventListener('focusout', onFocusOut, false);
  appRoot.addEventListener('click', onClickCapture, true);
  appRoot.addEventListener('click', onClick, false);

  return () => {
    appRoot.removeEventListener('pointerdown', onPointerDown, true);
    appRoot.removeEventListener('focusin', onFocusIn, true);
    appRoot.removeEventListener('keydown', onKeyDown, true);
    appRoot.removeEventListener('change', onChange, false);
    appRoot.removeEventListener('focusout', onFocusOut, false);
    appRoot.removeEventListener('click', onClickCapture, true);
    appRoot.removeEventListener('click', onClick, false);
  };
};
