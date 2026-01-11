<script setup lang="ts">
import { ref, computed } from 'vue'
import { sanitizeSvg } from '@/utils/sanitize'

const props = defineProps<{
  svgContent: string | null
  zoom: number
}>()

const emit = defineEmits<{
  'update:zoom': [value: number]
}>()

const containerRef = ref<HTMLDivElement | null>(null)

const sanitizedSvg = computed(() => {
  if (!props.svgContent) return null
  return sanitizeSvg(props.svgContent)
})

const zoomStyle = computed(() => ({
  transform: `scale(${props.zoom})`,
  transformOrigin: 'center center',
}))

function zoomIn() {
  emit('update:zoom', Math.min(props.zoom + 0.1, 3))
}

function zoomOut() {
  emit('update:zoom', Math.max(props.zoom - 0.1, 0.1))
}

function resetZoom() {
  emit('update:zoom', 1)
}
</script>

<template>
  <div class="flex-1 flex flex-col overflow-hidden">
    <!-- Zoom controls -->
    <div class="flex items-center justify-center gap-2 p-2 bg-white border-b border-slate-200">
      <button @click="zoomOut" class="btn btn-secondary btn-sm" :disabled="zoom <= 0.1">
        <i class="ph ph-minus"></i>
      </button>
      <span class="text-xs font-mono text-slate-500 w-16 text-center">{{ Math.round(zoom * 100) }}%</span>
      <button @click="zoomIn" class="btn btn-secondary btn-sm" :disabled="zoom >= 3">
        <i class="ph ph-plus"></i>
      </button>
      <button @click="resetZoom" class="btn btn-secondary btn-sm ml-2">
        <i class="ph ph-arrows-in"></i> Reset
      </button>
    </div>

    <!-- SVG container -->
    <div ref="containerRef" class="flex-1 overflow-auto bg-slate-100 flex items-center justify-center p-4">
      <div v-if="sanitizedSvg" :style="zoomStyle" class="transition-transform duration-200">
        <div v-html="sanitizedSvg" class="svg-container"></div>
      </div>
      <div v-else class="text-center text-slate-400">
        <i class="ph ph-image text-6xl mb-4"></i>
        <p class="text-lg font-medium">Preview Area</p>
        <p class="text-sm mt-2">Upload a genome file and click Generate to see your diagram.</p>
      </div>
    </div>
  </div>
</template>

<style scoped>
.svg-container :deep(svg) {
  max-width: 100%;
  height: auto;
  background: white;
  box-shadow: 0 4px 6px -1px rgb(0 0 0 / 0.1), 0 2px 4px -2px rgb(0 0 0 / 0.1);
  border-radius: 0.5rem;
}
</style>
