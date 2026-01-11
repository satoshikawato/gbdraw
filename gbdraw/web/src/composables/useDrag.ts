/**
 * Vue composable for drag functionality
 *
 * Handles drag and drop for legend and diagram elements in SVG.
 */

import { ref, reactive, onUnmounted } from 'vue'
import { parseTransform } from '@/utils/svg'

export interface DragState {
  dragging: boolean
  startX: number
  startY: number
  currentOffsetX: number
  currentOffsetY: number
}

export function useLegendDrag(getSvg: () => SVGSVGElement | null) {
  const dragging = ref(false)
  const dragStart = reactive({ x: 0, y: 0 })
  const originalTransform = ref({ x: 0, y: 0 })
  const initialTransform = ref({ x: 0, y: 0 })
  const currentOffset = reactive({ x: 0, y: 0 })

  let legendGroup: SVGGElement | null = null

  function startDrag(e: MouseEvent) {
    const svg = getSvg()
    if (!svg) return

    legendGroup = svg.getElementById('legend') as SVGGElement | null
    if (!legendGroup) return

    e.preventDefault()
    e.stopPropagation()

    dragging.value = true
    dragStart.x = e.clientX
    dragStart.y = e.clientY

    const currentTransform = parseTransform(legendGroup.getAttribute('transform'))
    originalTransform.value = { ...currentTransform }

    legendGroup.style.cursor = 'grabbing'
    legendGroup.style.opacity = '0.8'

    document.addEventListener('mousemove', onDrag)
    document.addEventListener('mouseup', endDrag)
  }

  function onDrag(e: MouseEvent) {
    if (!dragging.value || !legendGroup) return

    const svg = getSvg()
    if (!svg) return

    const dx = e.clientX - dragStart.x
    const dy = e.clientY - dragStart.y

    // Get SVG scale factor
    const svgRect = svg.getBoundingClientRect()
    const viewBox = svg.viewBox.baseVal
    const scaleX = viewBox.width / svgRect.width
    const scaleY = viewBox.height / svgRect.height

    const newX = originalTransform.value.x + dx * scaleX
    const newY = originalTransform.value.y + dy * scaleY

    legendGroup.setAttribute('transform', `translate(${newX}, ${newY})`)

    currentOffset.x = newX - initialTransform.value.x
    currentOffset.y = newY - initialTransform.value.y
  }

  function endDrag() {
    if (!dragging.value) return

    dragging.value = false

    if (legendGroup) {
      legendGroup.style.cursor = 'grab'
      legendGroup.style.opacity = '1'
    }

    document.removeEventListener('mousemove', onDrag)
    document.removeEventListener('mouseup', endDrag)
  }

  function setupDrag() {
    const svg = getSvg()
    if (!svg) return

    legendGroup = svg.getElementById('legend') as SVGGElement | null
    if (!legendGroup) return

    const transform = parseTransform(legendGroup.getAttribute('transform'))
    initialTransform.value = { ...transform }

    currentOffset.x = 0
    currentOffset.y = 0

    legendGroup.style.cursor = 'grab'
    legendGroup.removeEventListener('mousedown', startDrag)
    legendGroup.addEventListener('mousedown', startDrag)
  }

  function resetPosition() {
    const svg = getSvg()
    if (!svg) return

    legendGroup = svg.getElementById('legend') as SVGGElement | null
    if (!legendGroup) return

    legendGroup.setAttribute(
      'transform',
      `translate(${initialTransform.value.x}, ${initialTransform.value.y})`
    )
    currentOffset.x = 0
    currentOffset.y = 0
  }

  onUnmounted(() => {
    document.removeEventListener('mousemove', onDrag)
    document.removeEventListener('mouseup', endDrag)
  })

  return {
    dragging,
    currentOffset,
    setupDrag,
    resetPosition,
  }
}

export function useDiagramDrag(getSvg: () => SVGSVGElement | null, mode: () => 'circular' | 'linear') {
  const dragging = ref(false)
  const dragStart = reactive({ x: 0, y: 0 })
  const offset = reactive({ x: 0, y: 0 })
  const originalTransforms = ref<Map<Element, { x: number; y: number }>>(new Map())

  function getElementsToMove(): Element[] {
    const svg = getSvg()
    if (!svg) return []

    const elements: Element[] = []
    const ids = mode() === 'circular'
      ? ['tick', 'labels', 'axis', 'definition']
      : ['tick', 'labels', 'axis']

    ids.forEach((id) => {
      const el = svg.getElementById(id)
      if (el) elements.push(el)
    })

    // For linear mode, also get record groups
    if (mode() === 'linear') {
      const records = svg.querySelectorAll('[id^="record_"]')
      records.forEach((el) => elements.push(el))
    }

    return elements
  }

  function startDrag(e: MouseEvent) {
    const svg = getSvg()
    if (!svg) return

    // Check if clicking on a diagram element
    const target = e.target as Element
    const elements = getElementsToMove()
    const isOnElement = elements.some((el) => el.contains(target) || el === target)

    if (!isOnElement) return

    e.preventDefault()
    e.stopPropagation()

    dragging.value = true
    dragStart.x = e.clientX
    dragStart.y = e.clientY

    // Store original transforms
    originalTransforms.value.clear()
    elements.forEach((el) => {
      const transform = parseTransform(el.getAttribute('transform'))
      originalTransforms.value.set(el, transform)
    })

    elements.forEach((el) => {
      (el as SVGElement).style.cursor = 'grabbing'
    })

    document.addEventListener('mousemove', onDrag)
    document.addEventListener('mouseup', endDrag)
  }

  function onDrag(e: MouseEvent) {
    if (!dragging.value) return

    const svg = getSvg()
    if (!svg) return

    const dx = e.clientX - dragStart.x
    const dy = e.clientY - dragStart.y

    const svgRect = svg.getBoundingClientRect()
    const viewBox = svg.viewBox.baseVal
    const scaleX = viewBox.width / svgRect.width
    const scaleY = viewBox.height / svgRect.height

    originalTransforms.value.forEach((origTransform, el) => {
      const newX = origTransform.x + dx * scaleX
      const newY = origTransform.y + dy * scaleY
      el.setAttribute('transform', `translate(${newX}, ${newY})`)
    })

    offset.x = dx * scaleX
    offset.y = dy * scaleY
  }

  function endDrag() {
    if (!dragging.value) return

    dragging.value = false

    originalTransforms.value.forEach((_, el) => {
      (el as SVGElement).style.cursor = 'grab'
    })

    document.removeEventListener('mousemove', onDrag)
    document.removeEventListener('mouseup', endDrag)
  }

  function setupDrag() {
    const svg = getSvg()
    if (!svg) return

    const elements = getElementsToMove()

    // Store original transforms
    originalTransforms.value.clear()
    elements.forEach((el) => {
      const transform = parseTransform(el.getAttribute('transform'))
      originalTransforms.value.set(el, transform)
      ;(el as SVGElement).style.cursor = 'grab'
    })

    offset.x = 0
    offset.y = 0

    svg.removeEventListener('mousedown', startDrag)
    svg.addEventListener('mousedown', startDrag)
  }

  function resetPositions() {
    originalTransforms.value.forEach((origTransform, el) => {
      el.setAttribute('transform', `translate(${origTransform.x}, ${origTransform.y})`)
    })
    offset.x = 0
    offset.y = 0
  }

  onUnmounted(() => {
    document.removeEventListener('mousemove', onDrag)
    document.removeEventListener('mouseup', endDrag)
  })

  return {
    dragging,
    offset,
    setupDrag,
    resetPositions,
  }
}
