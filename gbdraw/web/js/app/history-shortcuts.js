const isTextEditingTarget = (target) => {
  if (!target) return false;
  if (target.isContentEditable) return true;
  const editable = target.closest?.('input, textarea, select, [contenteditable="true"]');
  if (!editable) return false;
  const tag = String(editable.tagName || '').toLowerCase();
  if (tag === 'textarea' || tag === 'select') return true;
  if (tag !== 'input') return Boolean(editable.isContentEditable);
  const type = String(editable.type || 'text').toLowerCase();
  return !['checkbox', 'radio', 'button', 'submit', 'reset', 'file', 'color'].includes(type);
};

export const setupHistoryShortcuts = ({ history, onMounted, onUnmounted }) => {
  const handleKeyDown = (event) => {
    if (!history) return;
    if (!(event.ctrlKey || event.metaKey)) return;
    if (event.altKey) return;
    if (isTextEditingTarget(event.target)) return;

    const key = String(event.key || '').toLowerCase();
    const wantsUndo = key === 'z' && !event.shiftKey;
    const wantsRedo = (key === 'z' && event.shiftKey) || key === 'y';
    if (!wantsUndo && !wantsRedo) return;

    event.preventDefault();
    if (wantsUndo) {
      void history.undo();
    } else {
      void history.redo();
    }
  };

  onMounted(() => {
    document.addEventListener('keydown', handleKeyDown);
  });
  onUnmounted(() => {
    document.removeEventListener('keydown', handleKeyDown);
  });
};
