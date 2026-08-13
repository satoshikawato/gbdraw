import {
  colorValueForMode,
  colorValueMode,
  toNativeColorInputValue
} from './app/color-utils.js';

const { ref, reactive, computed, nextTick, watch } = window.Vue;

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

export const ColorValueControl = {
  props: {
    modelValue: { default: null },
    fallback: { type: String, default: '#000000' },
    allowNone: { type: Boolean, default: true },
    ariaLabel: { type: String, default: 'Color value' }
  },
  emits: ['update:modelValue'],
  setup(props, { emit }) {
    const mode = computed(() => colorValueMode(props.modelValue));
    const pickerValue = computed(() => (
      toNativeColorInputValue(props.modelValue, props.fallback)
    ));
    const updateMode = (event) => {
      const nextMode = String(event?.target?.value || 'auto');
      emit(
        'update:modelValue',
        colorValueForMode(nextMode, props.modelValue, props.fallback)
      );
    };
    const updateColor = (event) => {
      emit(
        'update:modelValue',
        toNativeColorInputValue(event?.target?.value, props.fallback)
      );
    };
    return { mode, pickerValue, updateMode, updateColor };
  },
  template: `
    <div class="grid grid-cols-[minmax(0,1fr)_2.25rem] gap-1 items-center">
      <select
        :value="mode"
        @change="updateMode"
        class="form-input form-input-compact min-w-0"
        :aria-label="\`\${ariaLabel} mode\`"
      >
        <option value="auto">Auto</option>
        <option v-if="allowNone" value="none">None</option>
        <option value="color">Color</option>
      </select>
      <input
        type="color"
        :value="pickerValue"
        @input="updateColor"
        :disabled="mode !== 'color'"
        class="h-8 w-full p-0 border rounded disabled:opacity-40"
        :aria-label="ariaLabel"
      >
    </div>
  `
};

export const DefaultColorControl = {
  props: {
    modelValue: { default: null },
    fallback: { type: String, default: '#000000' },
    allowNone: { type: Boolean, default: true },
    ariaLabel: { type: String, default: 'Default color' }
  },
  emits: ['update:modelValue', 'accept'],
  setup(props, { emit }) {
    const draftColor = ref(toNativeColorInputValue(props.modelValue, props.fallback));
    const mode = computed(() => colorValueMode(props.modelValue));
    watch(
      () => [props.modelValue, props.fallback],
      () => {
        draftColor.value = toNativeColorInputValue(props.modelValue, props.fallback);
      }
    );
    const acceptValue = (value) => {
      emit('update:modelValue', value);
      emit('accept', value);
    };
    const updateMode = (event) => {
      const nextMode = String(event?.target?.value || 'auto');
      acceptValue(colorValueForMode(nextMode, props.modelValue, props.fallback));
    };
    const stageColor = (event) => {
      draftColor.value = toNativeColorInputValue(event?.target?.value, props.fallback);
    };
    const acceptColor = (event) => {
      stageColor(event);
      acceptValue(draftColor.value);
    };
    return { draftColor, mode, updateMode, stageColor, acceptColor };
  },
  template: `
    <div class="default-color-control grid grid-cols-[minmax(0,1fr)_2.25rem] gap-1 items-center">
      <select
        :value="mode"
        @change="updateMode"
        class="form-input form-input-compact min-w-0"
        :aria-label="\`\${ariaLabel} mode\`"
      >
        <option value="auto">Auto</option>
        <option v-if="allowNone" value="none">None</option>
        <option value="color">Color</option>
      </select>
      <input
        type="color"
        :value="draftColor"
        @input="stageColor"
        @change="acceptColor"
        :disabled="mode !== 'color'"
        class="h-8 w-full p-0 border rounded disabled:opacity-40"
        :aria-label="ariaLabel"
      >
    </div>
  `
};

export const FeatureFillColorControl = {
  props: {
    viewModel: { type: Object, required: true },
    ariaLabel: { type: String, default: 'Feature fill color' }
  },
  emits: ['request'],
  setup(props, { emit }) {
    const effectiveColor = computed(() => (
      toNativeColorInputValue(props.viewModel?.effectiveColor, '#cccccc')
    ));
    const mode = computed(() => (
      String(props.viewModel?.effectiveColor || '').toLowerCase() === 'none' ? 'none' : 'color'
    ));
    const requestMode = (event) => {
      const next = String(event?.target?.value || 'color');
      if (next === 'inherit') emit('request', { kind: 'inherit' });
      if (next === 'none') emit('request', { kind: 'none' });
    };
    const requestColor = (event) => {
      emit('request', {
        kind: 'color',
        color: toNativeColorInputValue(event?.target?.value, effectiveColor.value)
      });
    };
    return { effectiveColor, mode, requestMode, requestColor };
  },
  template: `
    <div class="feature-fill-color-control space-y-1">
      <div class="grid grid-cols-[minmax(0,1fr)_2.25rem] gap-1 items-center">
        <select
          :value="mode"
          @change="requestMode"
          class="form-input form-input-compact min-w-0"
          :aria-label="ariaLabel + ' action'"
        >
          <option value="color">Color</option>
          <option v-if="viewModel.allowNone !== false" value="none">No fill</option>
          <option value="inherit" :disabled="!viewModel.canReset">Use inherited</option>
        </select>
        <input
          type="color"
          :value="effectiveColor"
          @change="requestColor"
          class="h-8 w-full p-0 border rounded"
          :aria-label="ariaLabel"
        >
      </div>
      <p class="text-[10px] leading-tight text-slate-400" data-feature-fill-origin>
        {{ viewModel.originLabel }}
      </p>
    </div>
  `
};

export const FeatureStrokeColorControl = {
  props: {
    viewModel: { type: Object, required: true },
    ariaLabel: { type: String, default: 'Feature stroke' }
  },
  emits: ['request'],
  setup(props, { emit }) {
    const effectiveColor = computed(() => (
      toNativeColorInputValue(props.viewModel?.effectiveColor, '#000000')
    ));
    const effectiveWidth = computed(() => {
      const value = Number(props.viewModel?.effectiveWidth);
      return Number.isFinite(value) && value >= 0 ? value : 0.5;
    });
    const requestMode = (event) => {
      if (String(event?.target?.value || '') === 'inherit') {
        emit('request', { kind: 'inherit' });
      }
    };
    const requestColor = (event) => {
      emit('request', {
        kind: 'stroke',
        strokeColor: toNativeColorInputValue(event?.target?.value, effectiveColor.value)
      });
    };
    const requestWidth = (event) => {
      const strokeWidth = Number(event?.target?.value);
      if (Number.isFinite(strokeWidth) && strokeWidth >= 0) {
        emit('request', { kind: 'stroke', strokeWidth });
      }
    };
    return { effectiveColor, effectiveWidth, requestMode, requestColor, requestWidth };
  },
  template: `
    <div class="feature-stroke-color-control space-y-1">
      <div class="grid grid-cols-[minmax(0,1fr)_2.25rem_4.5rem] gap-1 items-center">
        <select
          value="stroke"
          @change="requestMode"
          class="form-input form-input-compact min-w-0"
          :aria-label="ariaLabel + ' action'"
        >
          <option value="stroke">Stroke</option>
          <option value="inherit" :disabled="!viewModel.canReset">Use inherited</option>
        </select>
        <input
          type="color"
          :value="effectiveColor"
          @change="requestColor"
          class="h-8 w-full p-0 border rounded"
          :aria-label="ariaLabel + ' color'"
        >
        <input
          type="number"
          :value="effectiveWidth"
          @change="requestWidth"
          min="0"
          step="0.1"
          class="form-input form-input-compact min-w-0"
          :aria-label="ariaLabel + ' width'"
        >
      </div>
      <p class="text-[10px] leading-tight text-slate-400" data-feature-stroke-origin>
        {{ viewModel.originLabel }}
      </p>
    </div>
  `
};

export const FileUploader = {
  template: '#file-uploader-template',
  props: ['label', 'accept', 'modelValue', 'small', 'multiple', 'testId'],
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
