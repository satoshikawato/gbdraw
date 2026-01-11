<script setup lang="ts">
import { HelpTip } from '@/components/common'

defineProps<{
  disabled: boolean
  exporting: boolean
  dpi: number
}>()

const emit = defineEmits<{
  'download-svg': []
  'download-png': []
  'download-pdf': []
  'update:dpi': [value: number]
}>()
</script>

<template>
  <div class="flex items-center gap-2 p-3 bg-white border-t border-slate-200">
    <span class="text-xs font-bold text-slate-500">Download:</span>

    <button
      @click="emit('download-svg')"
      :disabled="disabled || exporting"
      class="btn btn-secondary btn-sm"
    >
      <i class="ph ph-file-svg"></i> SVG
    </button>

    <button
      @click="emit('download-png')"
      :disabled="disabled || exporting"
      class="btn btn-secondary btn-sm"
    >
      <i class="ph ph-file-png"></i> PNG
    </button>

    <button
      @click="emit('download-pdf')"
      :disabled="disabled || exporting"
      class="btn btn-secondary btn-sm"
    >
      <i class="ph ph-file-pdf"></i> PDF
    </button>

    <div class="ml-auto flex items-center gap-2">
      <label class="text-xs text-slate-500 flex items-center gap-1">
        DPI
        <HelpTip text="Resolution for PNG export. Higher DPI = larger file size but better quality for print." />
      </label>
      <select
        :value="dpi"
        @change="emit('update:dpi', Number(($event.target as HTMLSelectElement).value))"
        class="form-input text-xs py-1 w-20"
        :disabled="exporting"
      >
        <option :value="72">72 (Web)</option>
        <option :value="150">150</option>
        <option :value="300">300 (Print)</option>
        <option :value="600">600 (HQ)</option>
      </select>
    </div>

    <!-- Loading indicator -->
    <div v-if="exporting" class="flex items-center gap-2 text-xs text-blue-600">
      <div class="animate-spin w-4 h-4 border-2 border-blue-600 border-t-transparent rounded-full"></div>
      Exporting...
    </div>
  </div>
</template>
