<script setup lang="ts">
import type { FormState, AdvancedSettings, DiagramMode } from '@/types'
import { HelpTip } from '@/components/common'

defineProps<{
  form: FormState
  adv: AdvancedSettings
  mode: DiagramMode
}>()

const emit = defineEmits<{
  'update:form': [value: Partial<FormState>]
  'update:adv': [value: Partial<AdvancedSettings>]
}>()

function updateForm<K extends keyof FormState>(key: K, value: FormState[K]) {
  emit('update:form', { [key]: value })
}

function updateAdv<K extends keyof AdvancedSettings>(key: K, value: AdvancedSettings[K]) {
  emit('update:adv', { [key]: value })
}
</script>

<template>
  <div class="space-y-4">
    <!-- Output Prefix & Legend -->
    <div class="grid grid-cols-2 gap-4">
      <div>
        <label class="input-label">
          Output Prefix
          <HelpTip text="Prefix for output files. Leave empty to use the record name (e.g., NC_000xxx)." />
        </label>
        <input
          type="text"
          :value="form.prefix"
          @input="updateForm('prefix', ($event.target as HTMLInputElement).value)"
          class="form-input"
          placeholder="Optional (Default: Record ID)"
        />
      </div>
      <div>
        <label class="input-label">
          Legend
          <HelpTip text="Position of the legend. 'None' hides the legend." />
        </label>
        <select
          :value="form.legend"
          @change="updateForm('legend', ($event.target as HTMLSelectElement).value as FormState['legend'])"
          class="form-input"
        >
          <option value="right">Right</option>
          <option value="left">Left</option>
          <option v-if="mode === 'linear'" value="top">Top</option>
          <option v-if="mode === 'linear'" value="bottom">Bottom</option>
          <option v-if="mode === 'circular'" value="upper_left">Upper Left</option>
          <option v-if="mode === 'circular'" value="upper_right">Upper Right</option>
          <option value="none">None</option>
        </select>
      </div>
    </div>

    <!-- Circular: Track Layout -->
    <div v-if="mode === 'circular'">
      <label class="input-label">
        Track Layout
        <HelpTip text="Choose how features are displayed. 'Tuckin' is default and compact. 'Middle' aligns to center, 'Spreadout' separates them." />
      </label>
      <select
        :value="form.track_type"
        @change="updateForm('track_type', ($event.target as HTMLSelectElement).value as FormState['track_type'])"
        class="form-input"
      >
        <option value="tuckin">Tuckin (Compact)</option>
        <option value="middle">Middle</option>
        <option value="spreadout">Spreadout</option>
      </select>
    </div>

    <!-- Linear: Scale Style -->
    <div v-if="mode === 'linear'">
      <label class="input-label">
        Scale Style
        <HelpTip text="'Bar' draws a simple line, 'Ruler' draws a ruler with ticks and labels." />
      </label>
      <select
        :value="form.scale_style"
        @change="updateForm('scale_style', ($event.target as HTMLSelectElement).value as FormState['scale_style'])"
        class="form-input"
      >
        <option value="bar">Bar (Simple)</option>
        <option value="ruler">Ruler (Ticks)</option>
      </select>
    </div>

    <!-- Checkboxes -->
    <div class="grid grid-cols-2 gap-y-2 gap-x-1 p-2 bg-slate-50 rounded-lg">
      <!-- Linear: Show Labels dropdown -->
      <div v-if="mode === 'linear'" class="col-span-2 mb-1">
        <label class="input-label text-[10px] mb-1">
          Show Labels
          <HelpTip text="Display feature labels on the map. 'First' shows labels only for the top track." />
        </label>
        <select
          :value="form.show_labels_linear"
          @change="updateForm('show_labels_linear', ($event.target as HTMLSelectElement).value as FormState['show_labels_linear'])"
          class="form-input text-xs py-1"
        >
          <option value="none">None</option>
          <option value="all">All Records</option>
          <option value="first">First Record Only</option>
        </select>
      </div>

      <!-- Circular: Show Labels checkbox -->
      <label v-else class="flex items-center gap-2 text-xs text-slate-700 cursor-pointer select-none">
        <input
          type="checkbox"
          :checked="form.show_labels"
          @change="updateForm('show_labels', ($event.target as HTMLInputElement).checked)"
          class="form-checkbox"
        />
        Show Labels
        <HelpTip text="Display feature labels on the map." />
      </label>

      <!-- Separate Strands -->
      <label class="flex items-center gap-2 text-xs text-slate-700 cursor-pointer select-none">
        <input
          type="checkbox"
          :checked="form.separate_strands"
          @change="updateForm('separate_strands', ($event.target as HTMLInputElement).checked)"
          class="form-checkbox"
        />
        Separate Strands
        <HelpTip text="Display features on separate strands for better distinction." />
      </label>

      <!-- Circular-specific options -->
      <template v-if="mode === 'circular'">
        <label class="flex items-center gap-2 text-xs text-slate-700 cursor-pointer select-none">
          <input
            type="checkbox"
            :checked="form.allow_inner_labels"
            @change="updateForm('allow_inner_labels', ($event.target as HTMLInputElement).checked)"
            class="form-checkbox"
          />
          Inner Labels
          <HelpTip text="Enable labels inside the circle. Automatically suppresses GC content/skew to avoid overlap." />
        </label>
        <label class="flex items-center gap-2 text-xs text-slate-700 cursor-pointer select-none">
          <input
            type="checkbox"
            :checked="form.suppress_gc"
            @change="updateForm('suppress_gc', ($event.target as HTMLInputElement).checked)"
            class="form-checkbox"
          />
          Hide GC
          <HelpTip text="Suppress the GC content track." />
        </label>
        <label class="flex items-center gap-2 text-xs text-slate-700 cursor-pointer select-none">
          <input
            type="checkbox"
            :checked="form.suppress_skew"
            @change="updateForm('suppress_skew', ($event.target as HTMLInputElement).checked)"
            class="form-checkbox"
          />
          Hide Skew
          <HelpTip text="Suppress the GC skew track." />
        </label>
      </template>

      <!-- Linear-specific options -->
      <template v-if="mode === 'linear'">
        <label class="flex items-center gap-2 text-xs text-slate-700 cursor-pointer select-none">
          <input
            type="checkbox"
            :checked="adv.resolve_overlaps"
            @change="updateAdv('resolve_overlaps', ($event.target as HTMLInputElement).checked)"
            class="form-checkbox"
          />
          Resolve Overlaps
          <HelpTip text="Shift features vertically to avoid overlap. (Experimental)" />
        </label>
        <label class="flex items-center gap-2 text-xs text-slate-700 cursor-pointer select-none">
          <input
            type="checkbox"
            :checked="form.align_center"
            @change="updateForm('align_center', ($event.target as HTMLInputElement).checked)"
            class="form-checkbox"
          />
          Align Center
          <HelpTip text="Align the linear map to the center of the page." />
        </label>
        <label class="flex items-center gap-2 text-xs text-slate-700 cursor-pointer select-none">
          <input
            type="checkbox"
            :checked="form.normalize_length"
            @change="updateForm('normalize_length', ($event.target as HTMLInputElement).checked)"
            class="form-checkbox"
          />
          Normalize Len
          <HelpTip text="Normalize sequence lengths to be equal. Suppresses length bar." />
        </label>
        <label class="flex items-center gap-2 text-xs text-slate-700 cursor-pointer select-none">
          <input
            type="checkbox"
            :checked="form.show_gc"
            @change="updateForm('show_gc', ($event.target as HTMLInputElement).checked)"
            class="form-checkbox"
          />
          Show GC
          <HelpTip text="Display GC content track." />
        </label>
        <label class="flex items-center gap-2 text-xs text-slate-700 cursor-pointer select-none">
          <input
            type="checkbox"
            :checked="form.show_skew"
            @change="updateForm('show_skew', ($event.target as HTMLInputElement).checked)"
            class="form-checkbox"
          />
          Show Skew
          <HelpTip text="Display GC skew track." />
        </label>
      </template>
    </div>
  </div>
</template>
