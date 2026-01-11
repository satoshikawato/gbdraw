<script setup lang="ts">
import { ref } from 'vue'

defineProps<{
  label: string
  accept?: string
  modelValue: File | null
  small?: boolean
}>()

const emit = defineEmits<{
  'update:modelValue': [value: File | null]
}>()

const input = ref<HTMLInputElement | null>(null)

function handleFile(e: Event) {
  const target = e.target as HTMLInputElement
  if (target.files?.[0]) {
    emit('update:modelValue', target.files[0])
  }
  target.value = '' // Reset input to allow re-selecting same file
}

function clearFile() {
  emit('update:modelValue', null)
}

function triggerInput() {
  input.value?.click()
}
</script>

<template>
  <div class="w-full">
    <div class="flex justify-between items-end mb-1">
      <label class="input-label mb-0">{{ label }}</label>
      <button
        v-if="modelValue"
        @click.stop="clearFile"
        class="text-[10px] text-red-500 hover:text-red-700 flex items-center gap-1"
      >
        <i class="ph ph-trash"></i> Remove
      </button>
    </div>
    <div
      class="upload-zone cursor-pointer"
      :class="{ ready: modelValue, 'py-1 min-h-[36px]': small }"
      @click="triggerInput"
    >
      <input
        type="file"
        ref="input"
        @change="handleFile"
        :accept="accept"
        class="hidden"
      />
      <div
        v-if="!modelValue"
        class="flex items-center justify-center gap-2 text-slate-400 group-hover:text-blue-600 transition-colors"
      >
        <i class="ph ph-upload-simple" :class="small ? 'text-base' : 'text-xl'"></i>
        <span :class="small ? 'text-[10px]' : 'text-xs font-bold'">Click to Browse</span>
      </div>
      <div
        v-else
        class="text-green-700 font-bold truncate w-full px-2 flex items-center justify-center gap-2"
        :class="small ? 'text-[10px]' : 'text-xs'"
      >
        <i class="ph ph-check-circle text-lg"></i>
        <span>{{ modelValue.name }}</span>
      </div>
    </div>
  </div>
</template>
