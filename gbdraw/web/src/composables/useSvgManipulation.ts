/**
 * Vue composable for SVG manipulation
 *
 * Handles SVG DOM manipulation including color application,
 * feature selection, and style updates.
 */

import { ref } from 'vue'
import type { SpecificRule } from '@/types'
import { parseViewBox, createViewBox } from '@/utils/svg'
import { interpolateColor, estimateColorFactor } from '@/utils/color'

export function useSvgManipulation() {
  const svgContainer = ref<HTMLDivElement | null>(null)

  /**
   * Get SVG element from container
   */
  function getSvg(): SVGSVGElement | null {
    return svgContainer.value?.querySelector('svg') ?? null
  }

  /**
   * Apply color palette to SVG features
   */
  function applyPaletteToSvg(colors: Record<string, string>) {
    const svg = getSvg()
    if (!svg) return

    // Find all feature paths and apply colors based on class
    Object.entries(colors).forEach(([featureType, color]) => {
      const elements = svg.querySelectorAll(`[class*="${featureType}"], [data-feature="${featureType}"]`)
      elements.forEach((el) => {
        if (el instanceof SVGElement) {
          el.setAttribute('fill', color)
        }
      })
    })
  }

  /**
   * Apply specific rules to SVG
   */
  function applySpecificRulesToSvg(rules: SpecificRule[]) {
    const svg = getSvg()
    if (!svg) return

    rules.forEach((rule) => {
      try {
        const regex = new RegExp(rule.val, 'i')
        const elements = svg.querySelectorAll('text')

        elements.forEach((textEl) => {
          const content = textEl.textContent || ''
          if (regex.test(content)) {
            // Find associated path element
            const parent = textEl.parentElement
            if (parent) {
              const paths = parent.querySelectorAll('path')
              paths.forEach((path) => {
                path.setAttribute('fill', rule.color)
              })
            }
          }
        })
      } catch {
        // Invalid regex, skip
      }
    })
  }

  /**
   * Apply stroke styles to SVG
   */
  function applyStylesToSvg(styles: {
    blockStrokeWidth?: number | null
    blockStrokeColor?: string | null
    lineStrokeWidth?: number | null
    lineStrokeColor?: string | null
    axisStrokeWidth?: number | null
    axisStrokeColor?: string | null
  }) {
    const svg = getSvg()
    if (!svg) return

    // Apply block stroke styles
    if (styles.blockStrokeWidth != null || styles.blockStrokeColor != null) {
      const blocks = svg.querySelectorAll('[class*="feature"], [class*="block"]')
      blocks.forEach((el) => {
        if (styles.blockStrokeWidth != null) {
          el.setAttribute('stroke-width', String(styles.blockStrokeWidth))
        }
        if (styles.blockStrokeColor != null) {
          el.setAttribute('stroke', styles.blockStrokeColor)
        }
      })
    }

    // Apply line stroke styles
    if (styles.lineStrokeWidth != null || styles.lineStrokeColor != null) {
      const lines = svg.querySelectorAll('line, [class*="line"]')
      lines.forEach((el) => {
        if (styles.lineStrokeWidth != null) {
          el.setAttribute('stroke-width', String(styles.lineStrokeWidth))
        }
        if (styles.lineStrokeColor != null) {
          el.setAttribute('stroke', styles.lineStrokeColor)
        }
      })
    }

    // Apply axis stroke styles
    if (styles.axisStrokeWidth != null || styles.axisStrokeColor != null) {
      const axis = svg.querySelectorAll('[id*="axis"], [class*="axis"]')
      axis.forEach((el) => {
        if (styles.axisStrokeWidth != null) {
          el.setAttribute('stroke-width', String(styles.axisStrokeWidth))
        }
        if (styles.axisStrokeColor != null) {
          el.setAttribute('stroke', styles.axisStrokeColor)
        }
      })
    }
  }

  /**
   * Apply track visibility
   */
  function applyTrackVisibility(options: {
    showGc?: boolean
    showSkew?: boolean
    suppressGc?: boolean
    suppressSkew?: boolean
  }) {
    const svg = getSvg()
    if (!svg) return

    // GC content track
    const gcTrack = svg.querySelector('#gc_content, [id*="gc_content"]')
    if (gcTrack) {
      const show = options.showGc === true || options.suppressGc === false
      gcTrack.setAttribute('visibility', show ? 'visible' : 'hidden')
    }

    // GC skew track
    const skewTrack = svg.querySelector('#gc_skew, [id*="gc_skew"]')
    if (skewTrack) {
      const show = options.showSkew === true || options.suppressSkew === false
      skewTrack.setAttribute('visibility', show ? 'visible' : 'hidden')
    }
  }

  /**
   * Set feature color by ID
   */
  function setFeatureColor(featureId: string, color: string) {
    const svg = getSvg()
    if (!svg) return false

    const element = svg.querySelector(`#${CSS.escape(featureId)}`)
    if (element) {
      element.setAttribute('fill', color)
      return true
    }
    return false
  }

  /**
   * Get feature color by ID
   */
  function getFeatureColor(featureId: string): string | null {
    const svg = getSvg()
    if (!svg) return null

    const element = svg.querySelector(`#${CSS.escape(featureId)}`)
    return element?.getAttribute('fill') ?? null
  }

  /**
   * Apply pairwise match colors (for linear comparison)
   */
  function applyPairwiseMatchColors(
    minColor: string,
    maxColor: string,
    factors: Record<string, number>
  ) {
    const svg = getSvg()
    if (!svg) return

    Object.entries(factors).forEach(([pathId, factor]) => {
      const element = svg.querySelector(`#${CSS.escape(pathId)}`)
      if (element) {
        const newColor = interpolateColor(minColor, maxColor, factor)
        element.setAttribute('fill', newColor)
      }
    })
  }

  /**
   * Extract pairwise match factors from current SVG colors
   */
  function extractPairwiseMatchFactors(
    minColor: string,
    maxColor: string
  ): Record<string, number> {
    const svg = getSvg()
    if (!svg) return {}

    const factors: Record<string, number> = {}
    const matchPaths = svg.querySelectorAll('[id^="match_"], [class*="pairwise"]')

    matchPaths.forEach((el) => {
      const id = el.getAttribute('id')
      const fill = el.getAttribute('fill')
      if (id && fill) {
        factors[id] = estimateColorFactor(fill, minColor, maxColor)
      }
    })

    return factors
  }

  /**
   * Apply canvas padding
   */
  function applyCanvasPadding(padding: {
    top: number
    right: number
    bottom: number
    left: number
  }) {
    const svg = getSvg()
    if (!svg) return

    const viewBox = parseViewBox(svg.getAttribute('viewBox'))
    if (!viewBox) return

    const newViewBox = createViewBox(
      viewBox.x - padding.left,
      viewBox.y - padding.top,
      viewBox.w + padding.left + padding.right,
      viewBox.h + padding.top + padding.bottom
    )

    svg.setAttribute('viewBox', newViewBox)
  }

  /**
   * Reset canvas padding to original
   */
  function resetCanvasPadding(originalViewBox: string) {
    const svg = getSvg()
    if (!svg) return

    svg.setAttribute('viewBox', originalViewBox)
  }

  /**
   * Get all legend groups
   */
  function getLegendGroups(): Element[] {
    const svg = getSvg()
    if (!svg) return []

    const groups: Element[] = []

    const hLegend = svg.getElementById('horizontal_feature_legend')
    if (hLegend) groups.push(hLegend)

    const vLegend = svg.getElementById('vertical_feature_legend')
    if (vLegend) groups.push(vLegend)

    if (groups.length === 0) {
      const legend = svg.getElementById('legend')
      if (legend) groups.push(legend)
    }

    return groups
  }

  /**
   * Get SVG content as string
   */
  function getSvgContent(): string | null {
    const svg = getSvg()
    if (!svg) return null

    const serializer = new XMLSerializer()
    return serializer.serializeToString(svg)
  }

  /**
   * Update SVG content
   */
  function updateSvgContent(content: string) {
    if (!svgContainer.value) return

    svgContainer.value.innerHTML = content
  }

  return {
    svgContainer,
    getSvg,
    applyPaletteToSvg,
    applySpecificRulesToSvg,
    applyStylesToSvg,
    applyTrackVisibility,
    setFeatureColor,
    getFeatureColor,
    applyPairwiseMatchColors,
    extractPairwiseMatchFactors,
    applyCanvasPadding,
    resetCanvasPadding,
    getLegendGroups,
    getSvgContent,
    updateSvgContent,
  }
}
