/**
 * SVG manipulation utilities
 */

export interface TransformXY {
  x: number
  y: number
}

/**
 * Parse translate transform from SVG element
 * Handles both "translate(x, y)" and "translate(x,y)" formats
 */
export function parseTransform(transformStr: string | null): TransformXY {
  if (!transformStr) return { x: 0, y: 0 }
  const match = transformStr.match(/translate\(\s*([-\d.]+)\s*,?\s*([-\d.]+)?\s*\)/)
  if (match) {
    return { x: parseFloat(match[1]) || 0, y: parseFloat(match[2]) || 0 }
  }
  return { x: 0, y: 0 }
}

/**
 * Parse X and Y from translate transform
 * Alias for parseTransform for compatibility
 */
export function parseTransformXY(transform: string | null): TransformXY {
  if (!transform) return { x: 0, y: 0 }
  const match = transform.match(/translate\(\s*([\d.-]+)\s*,\s*([\d.-]+)\s*\)/)
  return match ? { x: parseFloat(match[1]), y: parseFloat(match[2]) } : { x: 0, y: 0 }
}

/**
 * Parse only Y value from translate transform
 */
export function parseTransformY(transform: string | null): number | null {
  if (!transform) return null
  const match = transform.match(/translate\(\s*[\d.-]+\s*,\s*([\d.-]+)\s*\)/)
  return match ? parseFloat(match[1]) : null
}

/**
 * Create a translate transform string
 */
export function createTranslate(x: number, y: number): string {
  return `translate(${x}, ${y})`
}

/**
 * Parse viewBox from SVG element
 */
export function parseViewBox(viewBoxStr: string | null): { x: number; y: number; w: number; h: number } | null {
  if (!viewBoxStr) return null
  const parts = viewBoxStr.split(/[\s,]+/).map(Number)
  if (parts.length !== 4) return null
  return { x: parts[0], y: parts[1], w: parts[2], h: parts[3] }
}

/**
 * Create a viewBox string
 */
export function createViewBox(x: number, y: number, w: number, h: number): string {
  return `${x} ${y} ${w} ${h}`
}

/**
 * Get element by ID from SVG, handling special characters in IDs
 */
export function getSvgElementById(svg: SVGSVGElement, id: string): Element | null {
  // Escape special characters for CSS selector
  const escapedId = CSS.escape(id)
  return svg.querySelector(`#${escapedId}`)
}

/**
 * Get all legend groups from SVG
 * Returns both horizontal and vertical feature legend groups
 */
export function getAllFeatureLegendGroups(svg: SVGSVGElement): Element[] {
  const groups: Element[] = []

  // Check for horizontal legends (linear mode)
  const hLegend = svg.getElementById('horizontal_feature_legend')
  if (hLegend) groups.push(hLegend)

  // Check for vertical legends (circular mode or linear mode)
  const vLegend = svg.getElementById('vertical_feature_legend')
  if (vLegend) groups.push(vLegend)

  // Fallback to #legend for single legend mode
  if (groups.length === 0) {
    const legend = svg.getElementById('legend')
    if (legend) groups.push(legend)
  }

  return groups
}

/**
 * Get diagram group IDs that should move together
 */
export function getDiagramElementIds(mode: 'circular' | 'linear'): string[] {
  if (mode === 'circular') {
    return ['tick', 'labels', 'axis', 'definition']
  } else {
    // Linear mode: find all record groups
    return ['tick', 'labels', 'axis']
  }
}

/**
 * Get SVG dimensions from viewBox
 */
export function getSvgDimensions(svg: SVGSVGElement): { width: number; height: number } | null {
  const viewBox = parseViewBox(svg.getAttribute('viewBox'))
  if (viewBox) {
    return { width: viewBox.w, height: viewBox.h }
  }
  // Fallback to width/height attributes
  const width = parseFloat(svg.getAttribute('width') || '0')
  const height = parseFloat(svg.getAttribute('height') || '0')
  if (width && height) {
    return { width, height }
  }
  return null
}
