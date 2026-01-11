/**
 * Vue composable for Pyodide integration
 *
 * Provides reactive state and methods for interacting with the Pyodide runtime.
 */

import { ref, readonly, computed } from 'vue'
import { pyodideService, type ColorPalettes, type LegendEntrySvg, type ExtractedFeaturesResult } from '@/services/pyodideService'
import type { GbdrawResult } from '@/types'

export function usePyodide() {
  // Reactive state
  const ready = ref(false)
  const processing = ref(false)
  const loadingStatus = ref('Initializing...')
  const errorLog = ref<string | null>(null)

  // Color palettes
  const colorPalettes = ref<ColorPalettes>({})
  const paletteNames = computed(() =>
    Object.keys(colorPalettes.value).filter((k) => k !== 'title').sort()
  )

  /**
   * Initialize Pyodide and load gbdraw
   */
  async function initialize(): Promise<void> {
    if (ready.value) return

    try {
      await pyodideService.initialize((status) => {
        loadingStatus.value = status
      })

      // Load color palettes
      colorPalettes.value = pyodideService.getColorPalettes()

      ready.value = true
    } catch (e) {
      const message = e instanceof Error ? e.message : String(e)
      loadingStatus.value = `Startup Error: ${message}`
      console.error('Pyodide initialization error:', e)
    }
  }

  /**
   * Run gbdraw analysis
   */
  async function runAnalysis(
    mode: 'circular' | 'linear',
    args: string[]
  ): Promise<GbdrawResult[] | null> {
    if (!ready.value) return null

    processing.value = true
    errorLog.value = null

    try {
      const result = pyodideService.runGbdraw(mode, args)

      if ('error' in result) {
        errorLog.value = result.error
        return null
      }

      return result as GbdrawResult[]
    } catch (e) {
      const message = e instanceof Error ? e.message : String(e)
      errorLog.value = message
      return null
    } finally {
      processing.value = false
    }
  }

  /**
   * Get color palette by name
   */
  function getColorPalette(name: string): Record<string, string> {
    return colorPalettes.value[name] || {}
  }

  /**
   * Generate legend entry SVG elements
   */
  function generateLegendEntry(
    caption: string,
    color: string,
    yOffset: number,
    rectSize?: number,
    fontSize?: number,
    fontFamily?: string,
    xOffset?: number
  ): LegendEntrySvg {
    return pyodideService.generateLegendEntry(
      caption,
      color,
      yOffset,
      rectSize,
      fontSize,
      fontFamily,
      xOffset
    )
  }

  /**
   * Extract features from GenBank file
   */
  function extractFeatures(gbPath: string): ExtractedFeaturesResult {
    return pyodideService.extractFeatures(gbPath)
  }

  /**
   * Write file to virtual filesystem
   */
  async function writeFile(file: File, path: string): Promise<boolean> {
    return pyodideService.writeFileFromBlob(file, path)
  }

  /**
   * Write raw data to virtual filesystem
   */
  function writeRawFile(path: string, data: Uint8Array): void {
    pyodideService.writeFile(path, data)
  }

  /**
   * Delete file from virtual filesystem
   */
  function deleteFile(path: string): void {
    pyodideService.deleteFile(path)
  }

  return {
    // State (readonly)
    ready: readonly(ready),
    processing: readonly(processing),
    loadingStatus: readonly(loadingStatus),
    errorLog: readonly(errorLog),
    paletteNames,
    colorPalettes: readonly(colorPalettes),

    // Methods
    initialize,
    runAnalysis,
    getColorPalette,
    generateLegendEntry,
    extractFeatures,
    writeFile,
    writeRawFile,
    deleteFile,
  }
}
