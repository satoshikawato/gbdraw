<script setup lang="ts">
import { ref, reactive } from 'vue'

defineProps<{
  text: string
}>()

const visible = ref(false)
const style = reactive({
  top: '0px',
  left: '0px',
  transform: 'translate(-50%, -100%)',
})
const trigger = ref<HTMLElement | null>(null)

function show() {
  if (!trigger.value) return

  const rect = trigger.value.getBoundingClientRect()

  // Get viewport width and tooltip max width
  const viewportWidth = window.innerWidth
  const tooltipMaxWidth = 260 // max-w-250px + extra margin
  const halfWidth = tooltipMaxWidth / 2
  const gap = 12 // gap between icon and tooltip

  // Centered on the icon
  let left = rect.left + rect.width / 2
  let top = rect.top - gap
  let transform = 'translate(-50%, -100%)' // above the icon

  // Prevent overflow on left/right (keep at least 10px from viewport edges)
  if (left < halfWidth + 10) {
    left = halfWidth + 10
  } else if (left > viewportWidth - halfWidth - 10) {
    left = viewportWidth - halfWidth - 10
  }

  // Prevent overflow at the top (if too close to top edge, show below)
  // Considering header and browser frame, if y < 60px, show below
  if (rect.top < 60) {
    top = rect.bottom + gap
    transform = 'translate(-50%, 0)' // transform for below display
  }

  style.top = top + 'px'
  style.left = left + 'px'
  style.transform = transform

  visible.value = true
}

function hide() {
  visible.value = false
}
</script>

<template>
  <div class="inline-block ml-1 align-text-bottom">
    <i
      class="ph ph-question text-slate-400 hover:text-blue-500 cursor-help"
      @mouseenter="show"
      @mouseleave="hide"
      ref="trigger"
    ></i>
    <Teleport to="body">
      <div
        v-if="visible"
        :style="style"
        class="fixed z-[9999] px-3 py-2 bg-slate-800 text-white text-xs rounded-md shadow-lg max-w-[250px] pointer-events-none transition-opacity text-left leading-relaxed"
      >
        {{ text }}
      </div>
    </Teleport>
  </div>
</template>
