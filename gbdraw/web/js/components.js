const { ref, reactive } = window.Vue;

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

export const FileUploader = {
  template: '#file-uploader-template',
  props: ['label', 'accept', 'modelValue', 'small'],
  emits: ['update:modelValue'],
  setup(props, { emit }) {
    const input = ref(null);
    const handleFile = (e) => {
      if (e.target.files[0]) emit('update:modelValue', e.target.files[0]);
      e.target.value = '';
    };
    const clearFile = () => {
      emit('update:modelValue', null);
    };
    return { input, handleFile, clearFile };
  }
};
