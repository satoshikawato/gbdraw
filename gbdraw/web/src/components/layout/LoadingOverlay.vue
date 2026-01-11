<script setup lang="ts">
defineProps<{
  status: string
}>()
</script>

<template>
  <div class="fixed inset-0 bg-white z-[100] flex flex-col items-center justify-center p-4 text-center">
    <div class="animate-spin text-5xl mb-6">&#129516;</div>
    <h3 class="text-2xl font-bold text-slate-800">Initializing gbdraw Serverless...</h3>
    <p class="text-slate-500 mt-2">Loading Python engine (Pyodide) and dependencies.</p>
    <p class="text-xs text-slate-400 mt-4 font-mono">{{ status }}</p>
  </div>
</template>
