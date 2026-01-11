<script setup lang="ts">
import type { DiagramMode } from '@/types'

defineProps<{
  mode: DiagramMode
}>()

const emit = defineEmits<{
  'update:mode': [value: DiagramMode]
  'save-config': []
  'load-config': []
}>()

function setMode(newMode: DiagramMode) {
  emit('update:mode', newMode)
}
</script>

<template>
  <header class="bg-white border-b border-slate-200 sticky top-0 z-40 shadow-sm/50 backdrop-blur-md bg-white/90 h-16 shrink-0">
    <div class="max-w-[1600px] mx-auto px-4 h-full flex items-center justify-between">
      <!-- Logo and title -->
      <div class="flex items-center gap-3">
        <div class="bg-blue-600 text-white w-9 h-9 rounded-lg flex items-center justify-center font-bold text-xl shadow-md">
          &#129516;
        </div>
        <h1 class="text-xl font-bold text-slate-800 tracking-tight flex items-baseline">
          gbdraw
          <span class="text-slate-400 text-xs font-normal ml-3">
            A genome diagram generator for microbes and organelles
          </span>
        </h1>
      </div>

      <div class="flex items-center gap-3">
        <!-- Config buttons -->
        <div class="flex gap-1 mr-4">
          <button
            @click="emit('save-config')"
            class="text-xs font-bold text-slate-500 hover:text-blue-600 px-2 py-1 rounded border border-slate-200 bg-white flex items-center gap-1"
          >
            <i class="ph ph-export"></i> Save Config
          </button>
          <button
            @click="emit('load-config')"
            class="text-xs font-bold text-slate-500 hover:text-blue-600 px-2 py-1 rounded border border-slate-200 bg-white flex items-center gap-1"
          >
            <i class="ph ph-download-simple"></i> Load Config
          </button>
        </div>

        <!-- Mode toggle -->
        <div class="flex gap-1 bg-slate-100 p-1 rounded-lg">
          <button
            @click="setMode('circular')"
            :class="[
              'px-3 py-1 rounded text-sm font-medium transition-all',
              mode === 'circular' ? 'bg-white shadow text-blue-600' : 'text-slate-500 hover:text-slate-700',
            ]"
          >
            Circular
          </button>
          <button
            @click="setMode('linear')"
            :class="[
              'px-3 py-1 rounded text-sm font-medium transition-all',
              mode === 'linear' ? 'bg-white shadow text-blue-600' : 'text-slate-500 hover:text-slate-700',
            ]"
          >
            Linear
          </button>
        </div>
      </div>
    </div>
  </header>
</template>
