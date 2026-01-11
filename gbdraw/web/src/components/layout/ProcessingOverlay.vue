<script setup lang="ts">
defineProps<{
  message?: string
}>()
</script>

<template>
  <div class="fixed inset-0 bg-white/90 backdrop-blur-sm z-50 flex flex-col items-center justify-center">
    <div class="animate-spin w-16 h-16 border-4 border-blue-600 border-t-transparent rounded-full mb-6"></div>
    <h3 class="text-xl font-bold text-slate-800 animate-pulse">
      {{ message || 'Generating Diagram...' }}
    </h3>
  </div>
</template>
