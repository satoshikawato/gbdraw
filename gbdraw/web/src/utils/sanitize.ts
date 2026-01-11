/**
 * SVG sanitization utilities using DOMPurify
 */

import DOMPurify from 'dompurify'

/**
 * DOMPurify configuration for SVG sanitization
 * Tailored for svgwrite output from gbdraw
 */
const SVG_SANITIZE_CONFIG = {
  USE_PROFILES: { svg: true },
  // Only allow main tags used by svgwrite
  ADD_TAGS: [
    'use',
    'g',
    'defs',
    'linearGradient',
    'radialGradient',
    'stop',
    'path',
    'rect',
    'circle',
    'line',
    'polyline',
    'polygon',
    'text',
    'tspan',
  ],
  // Allowed attributes list
  ADD_ATTR: [
    'xlink:href',
    'href',
    'id',
    'class',
    'fill',
    'fill-opacity',
    'stroke',
    'stroke-width',
    'stroke-opacity',
    'stroke-dasharray',
    'stroke-linecap',
    'stroke-linejoin',
    'd',
    'x',
    'y',
    'x1',
    'y1',
    'x2',
    'y2',
    'cx',
    'cy',
    'r',
    'rx',
    'ry',
    'width',
    'height',
    'transform',
    'viewBox',
    'preserveAspectRatio',
    'font-family',
    'font-size',
    'font-weight',
    'text-anchor',
    'dominant-baseline',
    'writing-mode',
    'letter-spacing',
  ],
  // Forbid style tags (to prevent breaking parent page CSS) and script-related tags
  FORBID_TAGS: [
    'style',
    'script',
    'foreignObject',
    'iframe',
    'embed',
    'object',
    'animate',
    'set',
    'animateTransform',
    'image',
  ],
  FORBID_ATTR: ['name', 'onload', 'onclick', 'onmouseover', 'onfocus', 'onerror'],
}

/**
 * Sanitize SVG content for safe rendering
 */
export function sanitizeSvg(svg: string): string {
  return DOMPurify.sanitize(svg, SVG_SANITIZE_CONFIG)
}

/**
 * Check if content contains potentially dangerous patterns
 */
export function containsDangerousPatterns(content: string): boolean {
  // Check for script tags
  if (/<script/i.test(content)) return true
  // Check for event handlers
  if (/\bon\w+\s*=/i.test(content)) return true
  // Check for javascript: URLs
  if (/javascript:/i.test(content)) return true
  // Check for data: URLs with scripts
  if (/data:.*script/i.test(content)) return true
  return false
}
