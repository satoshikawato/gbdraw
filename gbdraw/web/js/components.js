const { ref, reactive, computed, nextTick } = window.Vue;

export const HelpTip = {
  template: '#help-tip-template',
  props: ['text'],
  setup() {
    const visible = ref(false);
    const style = reactive({ top: '0px', left: '0px' });
    const trigger = ref(null);
    const show = () => {
      if (!trigger.value) return;
      const rect = trigger.value.getBoundingClientRect();

      // Get viewport width and tooltip max width
      const viewportWidth = window.innerWidth;
      const tooltipMaxWidth = 260; // max-w-250px + extra margin
      const halfWidth = tooltipMaxWidth / 2;
      const gap = 12; // gap between icon and tooltip

      // Centered on the icon
      let left = rect.left + rect.width / 2;
      let top = rect.top - gap;
      let transform = 'translate(-50%, -100%)'; // above the icon

      // Prevent overflow on left/right (keep at least 10px from viewport edges)
      if (left < halfWidth + 10) {
        left = halfWidth + 10;
      } else if (left > viewportWidth - halfWidth - 10) {
        left = viewportWidth - halfWidth - 10;
      }

      // Prevent overflow at the top (if too close to top edge, show below)
      // Considering header and browser frame, if y < 60px, show below
      if (rect.top < 60) {
        top = rect.bottom + gap;
        transform = 'translate(-50%, 0)'; // transform for below display
      }

      style.top = `${top}px`;
      style.left = `${left}px`;
      style.transform = transform;

      visible.value = true;
    };
    const hide = () => {
      visible.value = false;
    };
    return { visible, style, trigger, show, hide };
  }
};

export const AutoValueField = {
  props: ['visible', 'text'],
  template: `
    <div class="auto-value-field">
      <span v-if="visible" class="auto-value-placeholder">{{ text }}</span>
      <slot></slot>
    </div>
  `
};

export const FileUploader = {
  template: '#file-uploader-template',
  props: ['label', 'accept', 'modelValue', 'small', 'multiple'],
  emits: ['update:modelValue'],
  setup(props, { emit }) {
    const input = ref(null);
    const selectedFiles = computed(() => {
      if (Array.isArray(props.modelValue)) return props.modelValue.filter(Boolean);
      return props.modelValue ? [props.modelValue] : [];
    });
    const hasSelection = computed(() => selectedFiles.value.length > 0);
    const selectedLabel = computed(() => {
      const items = selectedFiles.value;
      if (items.length === 0) return '';
      if (items.length === 1) return items[0]?.name || 'Selected file';
      const firstNames = items.slice(0, 2).map((file) => file?.name || 'file').join(', ');
      const suffix = items.length > 2 ? ` +${items.length - 2}` : '';
      return `${items.length} files: ${firstNames}${suffix}`;
    });
    const handleFile = (e) => {
      const nextFiles = Array.from(e.target.files || []);
      const update = () => {
        if (props.multiple) {
          emit('update:modelValue', nextFiles);
        } else if (nextFiles[0]) {
          emit('update:modelValue', nextFiles[0]);
        }
      };
      const history = window.__GBDRAW_HISTORY__;
      if (history?.runUndoable) {
        void history.runUndoable('Change uploaded file', async () => {
          update();
          await nextTick();
        });
      } else {
        update();
      }
      e.target.value = '';
    };
    const clearFile = () => {
      const history = window.__GBDRAW_HISTORY__;
      if (history?.runUndoable) {
        void history.runUndoable('Change uploaded file', async () => {
          emit('update:modelValue', props.multiple ? [] : null);
          await nextTick();
        });
      } else {
        emit('update:modelValue', props.multiple ? [] : null);
      }
    };
    return { input, handleFile, clearFile, hasSelection, selectedLabel };
  }
};
