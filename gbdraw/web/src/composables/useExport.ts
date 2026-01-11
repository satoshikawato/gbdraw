/**
 * Vue composable for export functionality
 *
 * Handles exporting diagrams to SVG, PNG, and PDF formats.
 */

import { ref } from 'vue'
import { setDpiInPng } from '@/utils/png'

export function useExport() {
  const downloadDpi = ref(300)
  const exporting = ref(false)

  /**
   * Download SVG file
   */
  function downloadSvg(svgContent: string, filename: string) {
    const blob = new Blob([svgContent], { type: 'image/svg+xml' })
    const url = URL.createObjectURL(blob)
    const link = document.createElement('a')
    link.href = url
    link.download = filename.endsWith('.svg') ? filename : `${filename}.svg`
    link.click()
    URL.revokeObjectURL(url)
  }

  /**
   * Download PNG file with DPI metadata
   */
  async function downloadPng(svgContent: string, filename: string, dpi?: number) {
    exporting.value = true
    const targetDpi = dpi ?? downloadDpi.value

    try {
      // Create canvas from SVG
      const parser = new DOMParser()
      const svgDoc = parser.parseFromString(svgContent, 'image/svg+xml')
      const svgElement = svgDoc.documentElement

      // Get SVG dimensions
      let width = parseFloat(svgElement.getAttribute('width') || '800')
      let height = parseFloat(svgElement.getAttribute('height') || '600')

      // Check viewBox for dimensions if width/height are not set
      const viewBox = svgElement.getAttribute('viewBox')
      if (viewBox) {
        const parts = viewBox.split(/[\s,]+/).map(Number)
        if (parts.length === 4) {
          width = parts[2]
          height = parts[3]
        }
      }

      // Scale for DPI
      const scale = targetDpi / 72
      const canvasWidth = Math.round(width * scale)
      const canvasHeight = Math.round(height * scale)

      // Create canvas
      const canvas = document.createElement('canvas')
      canvas.width = canvasWidth
      canvas.height = canvasHeight
      const ctx = canvas.getContext('2d')

      if (!ctx) {
        throw new Error('Could not get canvas context')
      }

      // Create image from SVG
      const img = new Image()
      const svgBlob = new Blob([svgContent], { type: 'image/svg+xml;charset=utf-8' })
      const svgUrl = URL.createObjectURL(svgBlob)

      await new Promise<void>((resolve, reject) => {
        img.onload = () => {
          ctx.fillStyle = 'white'
          ctx.fillRect(0, 0, canvasWidth, canvasHeight)
          ctx.drawImage(img, 0, 0, canvasWidth, canvasHeight)
          URL.revokeObjectURL(svgUrl)
          resolve()
        }
        img.onerror = () => {
          URL.revokeObjectURL(svgUrl)
          reject(new Error('Failed to load SVG image'))
        }
        img.src = svgUrl
      })

      // Convert to blob
      const blob = await new Promise<Blob>((resolve, reject) => {
        canvas.toBlob(
          (b) => {
            if (b) resolve(b)
            else reject(new Error('Failed to create PNG blob'))
          },
          'image/png',
          1.0
        )
      })

      // Add DPI metadata
      const fixedBlob = await setDpiInPng(blob, targetDpi)

      // Download
      const url = URL.createObjectURL(fixedBlob)
      const link = document.createElement('a')
      link.href = url
      link.download = filename.replace(/\.svg$/i, '.png')
      link.click()
      URL.revokeObjectURL(url)
    } finally {
      exporting.value = false
    }
  }

  /**
   * Download PDF file
   */
  async function downloadPdf(svgContent: string, filename: string) {
    exporting.value = true

    try {
      // Dynamic import of jspdf and svg2pdf
      const [{ jsPDF }, { svg2pdf }] = await Promise.all([
        import('jspdf'),
        import('svg2pdf.js'),
      ])

      // Parse SVG
      const parser = new DOMParser()
      const svgDoc = parser.parseFromString(svgContent, 'image/svg+xml')
      const svgElement = svgDoc.documentElement as unknown as SVGSVGElement

      // Get SVG dimensions
      let width = parseFloat(svgElement.getAttribute('width') || '800')
      let height = parseFloat(svgElement.getAttribute('height') || '600')

      const viewBox = svgElement.getAttribute('viewBox')
      if (viewBox) {
        const parts = viewBox.split(/[\s,]+/).map(Number)
        if (parts.length === 4) {
          width = parts[2]
          height = parts[3]
        }
      }

      // Create PDF with SVG dimensions (in points, 72 dpi)
      const orientation = width > height ? 'landscape' : 'portrait'
      const pdf = new jsPDF({
        orientation,
        unit: 'pt',
        format: [width, height],
      })

      // Convert SVG to PDF
      await svg2pdf(svgElement, pdf, {
        x: 0,
        y: 0,
        width,
        height,
      })

      // Save PDF
      pdf.save(filename.replace(/\.svg$/i, '.pdf'))
    } finally {
      exporting.value = false
    }
  }

  /**
   * Export configuration as JSON file
   */
  function exportConfig(config: object, filename = 'gbdraw-config.json') {
    const json = JSON.stringify(config, null, 2)
    const blob = new Blob([json], { type: 'application/json' })
    const url = URL.createObjectURL(blob)
    const link = document.createElement('a')
    link.href = url
    link.download = filename
    link.click()
    URL.revokeObjectURL(url)
  }

  /**
   * Import configuration from JSON file
   */
  async function importConfig(file: File): Promise<Record<string, unknown>> {
    const text = await file.text()
    const config = JSON.parse(text)

    // Security: Check for prototype pollution
    if ('__proto__' in config || 'constructor' in config || 'prototype' in config) {
      throw new Error('Invalid configuration file: contains forbidden keys')
    }

    return config
  }

  return {
    downloadDpi,
    exporting,
    downloadSvg,
    downloadPng,
    downloadPdf,
    exportConfig,
    importConfig,
  }
}
