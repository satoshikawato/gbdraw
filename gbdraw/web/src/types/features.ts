/**
 * Type definitions for feature and legend state
 */

export interface ExtractedFeature {
  id: string
  svg_id: string
  label: string
  color: string
  feat: string
  qualifier?: Record<string, string>
}

export interface ClickedFeature {
  id: string
  svg_id: string
  label: string
  color: string
  feat: string
}

export interface LegendEntry {
  caption: string
  color: string
  yPos: number
  id?: string
}

export interface CanvasPadding {
  top: number
  right: number
  bottom: number
  left: number
}

export interface DragPosition {
  x: number
  y: number
}

export interface CircularBaseConfig {
  viewBoxWidth: number
  viewBoxHeight: number
  diagramCenterX: number
  diagramCenterY: number
  legendWidth: number
  legendHeight: number
  generatedPosition: string
}

export interface LinearViewBox {
  x: number
  y: number
  w: number
  h: number
}

export interface LinearBaseConfig {
  verticalViewBox: LinearViewBox
  horizontalViewBox: LinearViewBox
  diagramBaseTransforms: Map<string, string>
  legendBaseTransform: string
  legendWidth: number
  legendHeight: number
  generatedPosition: string
}

export interface GbdrawResult {
  success: boolean
  svg?: string
  error?: string
  name?: string
}
