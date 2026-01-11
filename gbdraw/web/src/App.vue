<script setup lang="ts">
import { ref, computed, onMounted, watch } from 'vue'
import { usePyodide, useFormState, useExport } from '@/composables'
import {
  LoadingOverlay,
  ProcessingOverlay,
  AppHeader,
} from '@/components/layout'
import { InputGenomes, BasicSettings, ColorsFilters } from '@/components/settings'
import { SvgPreview, DownloadButtons } from '@/components/preview'
import type { GbdrawResult } from '@/types'

// Pyodide composable
const {
  ready: pyodideReady,
  processing,
  loadingStatus,
  errorLog,
  paletteNames,
  initialize,
  runAnalysis,
  getColorPalette,
  writeFile,
  writeRawFile,
} = usePyodide()

// Form state composable
const formState = useFormState()

// Export composable
const { downloadDpi, exporting, downloadSvg, downloadPng, downloadPdf, exportConfig, importConfig } = useExport()

// Results state
const results = ref<GbdrawResult[]>([])
const selectedResultIndex = ref(0)
const zoom = ref(1.0)

const svgContent = computed(() => {
  if (results.value.length > 0 && results.value[selectedResultIndex.value]) {
    return results.value[selectedResultIndex.value].svg || null
  }
  return null
})

const currentFilename = computed(() => {
  if (results.value.length > 0 && results.value[selectedResultIndex.value]) {
    return results.value[selectedResultIndex.value].name || 'diagram.svg'
  }
  return 'diagram.svg'
})

// Initialize Pyodide on mount
onMounted(async () => {
  await initialize()
  // Load default palette
  formState.setColors(getColorPalette('default'))
})

// Watch palette changes
watch(
  () => formState.selectedPalette.value,
  (newPalette) => {
    formState.setColors(getColorPalette(newPalette))
  }
)

// Helper: Convert colors to TOML format
function colorsToToml(colors: Record<string, string>): string {
  const lines = ['[default]']
  for (const [key, value] of Object.entries(colors)) {
    lines.push(`${key} = "${value}"`)
  }
  return lines.join('\n')
}

// Helper: Convert specific rules to TOML format
function specificRulesToToml(rules: typeof formState.specificRules): string {
  if (rules.length === 0) return ''
  const lines: string[] = []
  for (const rule of rules) {
    lines.push(`[[specific]]`)
    lines.push(`feature = "${rule.feat}"`)
    lines.push(`qualifier = "${rule.qual}"`)
    lines.push(`value = "${rule.val}"`)
    lines.push(`color = "${rule.color}"`)
    lines.push('')
  }
  return lines.join('\n')
}

// Run analysis
async function generate() {
  const args: string[] = []

  if (formState.mode.value === 'circular') {
    if (formState.cInputType.value === 'gb' && formState.files.c_gb) {
      await writeFile(formState.files.c_gb, '/input.gb')
      args.push('--gbk', '/input.gb')
    } else if (
      formState.cInputType.value === 'gff' &&
      formState.files.c_gff &&
      formState.files.c_fasta
    ) {
      await writeFile(formState.files.c_gff, '/input.gff')
      await writeFile(formState.files.c_fasta, '/input.fasta')
      args.push('--gff', '/input.gff', '--fasta', '/input.fasta')
    }
  } else {
    // Linear mode
    for (let i = 0; i < formState.linearSeqs.length; i++) {
      const seq = formState.linearSeqs[i]
      if (formState.lInputType.value === 'gb' && seq.gb) {
        await writeFile(seq.gb, `/seq_${i}.gb`)
        args.push('--gbk', `/seq_${i}.gb`)
      } else if (formState.lInputType.value === 'gff' && seq.gff && seq.fasta) {
        await writeFile(seq.gff, `/seq_${i}.gff`)
        await writeFile(seq.fasta, `/seq_${i}.fasta`)
        args.push('--gff', `/seq_${i}.gff`, '--fasta', `/seq_${i}.fasta`)
      }
      if (seq.blast && i < formState.linearSeqs.length - 1) {
        await writeFile(seq.blast, `/blast_${i}.txt`)
        args.push('-b', `/blast_${i}.txt`)
      }
    }
  }

  // Add common options
  if (formState.form.prefix) args.push('-p', formState.form.prefix)
  if (formState.form.legend !== 'right') args.push('-l', formState.form.legend)
  if (formState.form.track_type !== 'tuckin') args.push('--track_type', formState.form.track_type)
  if (formState.form.show_labels) args.push('--show_labels')
  if (!formState.form.separate_strands) args.push('--no_separate_strands')
  if (formState.form.suppress_gc) args.push('--suppress_gc')
  if (formState.form.suppress_skew) args.push('--suppress_skew')

  // Write color table as TOML and pass with -d
  const colorToml = colorsToToml(formState.currentColors.value)
  const colorData = new TextEncoder().encode(colorToml)
  writeRawFile('/colors.toml', colorData)
  args.push('-d', '/colors.toml')

  // Write specific rules as TOML and pass with -t (if any)
  if (formState.specificRules.length > 0) {
    const rulesText = specificRulesToToml(formState.specificRules)
    const rulesData = new TextEncoder().encode(rulesText)
    writeRawFile('/specific.toml', rulesData)
    args.push('-t', '/specific.toml')
  }

  const result = await runAnalysis(formState.mode.value, args)
  if (result) {
    results.value = result
    selectedResultIndex.value = 0
  }
}

// Config handlers
function handleSaveConfig() {
  const config = formState.exportConfig()
  exportConfig(config)
}

async function handleLoadConfig() {
  const input = document.createElement('input')
  input.type = 'file'
  input.accept = '.json'
  input.onchange = async (e) => {
    const file = (e.target as HTMLInputElement).files?.[0]
    if (file) {
      try {
        const config = await importConfig(file)
        formState.importConfig(config)
      } catch (err) {
        console.error('Failed to load config:', err)
        alert('Failed to load configuration file')
      }
    }
  }
  input.click()
}

// Export handlers
function handleDownloadSvg() {
  if (svgContent.value) {
    downloadSvg(svgContent.value, currentFilename.value)
  }
}

async function handleDownloadPng() {
  if (svgContent.value) {
    await downloadPng(svgContent.value, currentFilename.value)
  }
}

async function handleDownloadPdf() {
  if (svgContent.value) {
    await downloadPdf(svgContent.value, currentFilename.value)
  }
}
</script>

<template>
  <!-- Loading overlay -->
  <LoadingOverlay v-if="!pyodideReady" :status="loadingStatus" />

  <!-- Processing overlay -->
  <ProcessingOverlay v-if="processing" />

  <!-- Main app -->
  <div v-if="pyodideReady" class="h-screen flex flex-col relative">
    <!-- Header -->
    <AppHeader
      :mode="formState.mode.value"
      @update:mode="formState.mode.value = $event"
      @save-config="handleSaveConfig"
      @load-config="handleLoadConfig"
    />

    <!-- Main content -->
    <main class="flex-1 flex overflow-hidden">
      <!-- Left panel (settings) -->
      <div class="w-96 border-r border-slate-200 bg-white overflow-y-auto p-4 space-y-4">
        <!-- Input Genomes -->
        <InputGenomes
          :mode="formState.mode.value"
          :cInputType="formState.cInputType.value"
          :lInputType="formState.lInputType.value"
          :files="formState.files"
          :linearSeqs="formState.linearSeqs"
          @update:cInputType="formState.cInputType.value = $event"
          @update:lInputType="formState.lInputType.value = $event"
          @update:files="formState.updateFiles"
          @update:linearSeq="formState.updateLinearSeq"
          @add-linear-seq="formState.addLinearSeq"
          @remove-linear-seq="formState.removeLinearSeq"
        />

        <!-- Basic Settings -->
        <div class="card">
          <div class="card-header">
            <i class="ph ph-sliders"></i> Basic Settings
          </div>
          <BasicSettings
            :form="formState.form"
            :adv="formState.adv"
            :mode="formState.mode.value"
            @update:form="formState.updateForm"
            @update:adv="formState.updateAdv"
          />
        </div>

        <!-- Colors & Filters -->
        <ColorsFilters
          :paletteNames="paletteNames"
          :selectedPalette="formState.selectedPalette.value"
          :currentColors="formState.currentColors.value"
          :specificRules="formState.specificRules"
          :newSpecRule="formState.newSpecRule"
          :filterMode="formState.filterMode.value"
          :manualBlacklist="formState.manualBlacklist.value"
          :manualWhitelist="formState.manualWhitelist"
          :priorityRules="formState.priorityRules"
          :newPriorityRule="formState.newPriorityRule"
          :customFeature="formState.customFeature"
          @update:selectedPalette="formState.selectedPalette.value = $event; formState.setColors(getColorPalette($event))"
          @update:currentColors="formState.updateColor"
          @update:filterMode="formState.filterMode.value = $event"
          @update:manualBlacklist="formState.manualBlacklist.value = $event"
          @update:newSpecRule="(formState.newSpecRule as Record<string, string>)[$event[0]] = $event[1]"
          @update:newPriorityRule="(formState.newPriorityRule as Record<string, string>)[$event[0]] = $event[1]"
          @update:customFeature="(formState.customFeature as Record<string, string>)[$event[0]] = $event[1]"
          @reset-colors="formState.setColors(getColorPalette(formState.selectedPalette.value))"
          @add-custom-feature="formState.addCustomFeature"
          @add-specific-rule="formState.addSpecificRule"
          @remove-specific-rule="formState.removeSpecificRule"
          @add-priority-rule="formState.addPriorityRule"
          @remove-priority-rule="formState.removePriorityRule"
          @add-whitelist-item="formState.addWhitelistItem"
          @remove-whitelist-item="formState.removeWhitelistItem"
        />

        <!-- Generate Button -->
        <button
          @click="generate"
          :disabled="processing || !formState.hasInput.value"
          class="btn btn-generate w-full text-lg py-3"
        >
          <i class="ph ph-play-circle text-xl"></i>
          Generate Diagram
        </button>

        <!-- Error display -->
        <div v-if="errorLog" class="card bg-red-50 border-red-200">
          <div class="card-header text-red-600">
            <i class="ph ph-warning"></i> Error
          </div>
          <pre class="text-xs text-red-700 whitespace-pre-wrap overflow-auto max-h-48">{{ errorLog }}</pre>
        </div>
      </div>

      <!-- Right panel (preview) -->
      <div class="flex-1 flex flex-col overflow-hidden">
        <SvgPreview :svgContent="svgContent" :zoom="zoom" @update:zoom="zoom = $event" />
        <DownloadButtons
          :disabled="!svgContent"
          :exporting="exporting"
          :dpi="downloadDpi"
          @update:dpi="downloadDpi = $event"
          @download-svg="handleDownloadSvg"
          @download-png="handleDownloadPng"
          @download-pdf="handleDownloadPdf"
        />
      </div>
    </main>
  </div>
</template>
