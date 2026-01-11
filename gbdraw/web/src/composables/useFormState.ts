/**
 * Vue composable for form state management
 *
 * Centralizes all form-related state including basic settings,
 * advanced options, colors, and filters.
 */

import { reactive, ref, computed } from 'vue'
import type {
  DiagramMode,
  InputType,
  FormState,
  AdvancedSettings,
  FileState,
  LinearSeq,
  SpecificRule,
  PriorityRule,
  FilterMode,
} from '@/types'
import { createDefaultFormState, createDefaultAdvancedSettings } from '@/types'

export function useFormState() {
  // Mode state
  const mode = ref<DiagramMode>('circular')
  const cInputType = ref<InputType>('gb')
  const lInputType = ref<InputType>('gb')

  // Form state
  const form = reactive<FormState>(createDefaultFormState())
  const adv = reactive<AdvancedSettings>(createDefaultAdvancedSettings())

  // File state
  const files = reactive<FileState>({
    c_gb: null,
    c_gff: null,
    c_fasta: null,
    d_color: null,
    t_color: null,
    blacklist: null,
    whitelist: null,
    qualifier_priority: null,
  })

  const linearSeqs = reactive<LinearSeq[]>([
    { gb: null, gff: null, fasta: null, blast: null },
  ])

  // Color state
  const selectedPalette = ref('default')
  const currentColors = ref<Record<string, string>>({})
  const specificRules = reactive<SpecificRule[]>([])
  const newSpecRule = reactive<SpecificRule>({
    feat: 'CDS',
    qual: 'product',
    val: '',
    color: '#ff0000',
    cap: '',
  })

  // Filter state
  const filterMode = ref<FilterMode>('None')
  const manualBlacklist = ref('hypothetical, uncharacterized, putative, unknown')
  const manualWhitelist = reactive<string[]>([])
  const priorityRules = reactive<PriorityRule[]>([])
  const newPriorityRule = reactive<PriorityRule>({
    feat: 'CDS',
    order: 'product,gene,locus_tag',
  })

  // Custom feature colors
  const customFeature = reactive({ name: '', color: '#808080' })

  // Computed: check if input files are ready
  const hasCircularInput = computed(() => {
    if (cInputType.value === 'gb') {
      return files.c_gb !== null
    }
    return files.c_gff !== null && files.c_fasta !== null
  })

  const hasLinearInput = computed(() => {
    return linearSeqs.some((seq) => {
      if (lInputType.value === 'gb') {
        return seq.gb !== null
      }
      return seq.gff !== null && seq.fasta !== null
    })
  })

  const hasInput = computed(() => {
    return mode.value === 'circular' ? hasCircularInput.value : hasLinearInput.value
  })

  // Methods
  function updateForm(updates: Partial<FormState>) {
    Object.assign(form, updates)
  }

  function updateAdv(updates: Partial<AdvancedSettings>) {
    Object.assign(adv, updates)
  }

  function updateFiles(key: keyof FileState, value: File | null) {
    files[key] = value
  }

  function updateLinearSeq(index: number, key: keyof LinearSeq, value: File | null) {
    if (linearSeqs[index]) {
      linearSeqs[index][key] = value
    }
  }

  function addLinearSeq() {
    linearSeqs.push({ gb: null, gff: null, fasta: null, blast: null })
  }

  function removeLinearSeq() {
    if (linearSeqs.length > 1) {
      linearSeqs.pop()
    }
  }

  function setColors(colors: Record<string, string>) {
    currentColors.value = { ...colors }
  }

  function updateColor(feat: string, color: string) {
    currentColors.value[feat] = color
  }

  function addCustomFeature() {
    if (customFeature.name && !currentColors.value[customFeature.name]) {
      currentColors.value[customFeature.name] = customFeature.color
      customFeature.name = ''
      customFeature.color = '#808080'
    }
  }

  function addSpecificRule() {
    if (!newSpecRule.val) return

    // Warn about long regex patterns
    if (newSpecRule.val.length > 50) {
      if (!confirm('Regular expression is quite long (>50 chars). This might impact performance. Continue?')) {
        return
      }
    }

    // Check for potentially dangerous patterns (ReDoS)
    if (/\(.+[+*]\)[+*]/.test(newSpecRule.val) || /\(.*\)\+/.test(newSpecRule.val)) {
      if (!confirm('This regular expression contains patterns that may freeze the browser (ReDoS risk). Are you sure you want to add it?')) {
        return
      }
    }

    try {
      new RegExp(newSpecRule.val)
    } catch {
      alert('Invalid regular expression')
      return
    }

    specificRules.push({ ...newSpecRule })
    newSpecRule.val = ''
    newSpecRule.cap = ''
  }

  function removeSpecificRule(index: number) {
    specificRules.splice(index, 1)
  }

  function addPriorityRule() {
    if (!newPriorityRule.order) return

    const idx = priorityRules.findIndex((r) => r.feat === newPriorityRule.feat)
    if (idx >= 0) {
      priorityRules[idx].order = newPriorityRule.order
    } else {
      priorityRules.push({ ...newPriorityRule })
    }
  }

  function removePriorityRule(index: number) {
    priorityRules.splice(index, 1)
  }

  function addWhitelistItem(item: string) {
    if (item && !manualWhitelist.includes(item)) {
      manualWhitelist.push(item)
    }
  }

  function removeWhitelistItem(index: number) {
    manualWhitelist.splice(index, 1)
  }

  function resetState() {
    Object.assign(form, createDefaultFormState())
    Object.assign(adv, createDefaultAdvancedSettings())
    specificRules.length = 0
    priorityRules.length = 0
    manualWhitelist.length = 0
    filterMode.value = 'None'
    manualBlacklist.value = 'hypothetical, uncharacterized, putative, unknown'
  }

  // Export config as JSON
  function exportConfig(): object {
    return {
      mode: mode.value,
      cInputType: cInputType.value,
      lInputType: lInputType.value,
      form: { ...form },
      adv: { ...adv },
      selectedPalette: selectedPalette.value,
      currentColors: { ...currentColors.value },
      specificRules: [...specificRules],
      filterMode: filterMode.value,
      manualBlacklist: manualBlacklist.value,
      manualWhitelist: [...manualWhitelist],
      priorityRules: [...priorityRules],
    }
  }

  // Import config from JSON
  function importConfig(config: Record<string, unknown>) {
    if (config.mode) mode.value = config.mode as DiagramMode
    if (config.cInputType) cInputType.value = config.cInputType as InputType
    if (config.lInputType) lInputType.value = config.lInputType as InputType
    if (config.form) Object.assign(form, config.form)
    if (config.adv) Object.assign(adv, config.adv)
    if (config.selectedPalette) selectedPalette.value = config.selectedPalette as string
    if (config.currentColors) currentColors.value = config.currentColors as Record<string, string>
    if (config.specificRules) {
      specificRules.length = 0
      specificRules.push(...(config.specificRules as SpecificRule[]))
    }
    if (config.filterMode) filterMode.value = config.filterMode as FilterMode
    if (config.manualBlacklist) manualBlacklist.value = config.manualBlacklist as string
    if (config.manualWhitelist) {
      manualWhitelist.length = 0
      manualWhitelist.push(...(config.manualWhitelist as string[]))
    }
    if (config.priorityRules) {
      priorityRules.length = 0
      priorityRules.push(...(config.priorityRules as PriorityRule[]))
    }
  }

  return {
    // State
    mode,
    cInputType,
    lInputType,
    form,
    adv,
    files,
    linearSeqs,
    selectedPalette,
    currentColors,
    specificRules,
    newSpecRule,
    filterMode,
    manualBlacklist,
    manualWhitelist,
    priorityRules,
    newPriorityRule,
    customFeature,

    // Computed
    hasInput,
    hasCircularInput,
    hasLinearInput,

    // Methods
    updateForm,
    updateAdv,
    updateFiles,
    updateLinearSeq,
    addLinearSeq,
    removeLinearSeq,
    setColors,
    updateColor,
    addCustomFeature,
    addSpecificRule,
    removeSpecificRule,
    addPriorityRule,
    removePriorityRule,
    addWhitelistItem,
    removeWhitelistItem,
    resetState,
    exportConfig,
    importConfig,
  }
}
