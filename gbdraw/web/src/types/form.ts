/**
 * Type definitions for form state
 */

export type DiagramMode = 'circular' | 'linear'
export type InputType = 'gb' | 'gff'
export type TrackType = 'tuckin' | 'middle' | 'spreadout'
export type LegendPosition = 'right' | 'left' | 'top' | 'bottom' | 'upper_left' | 'upper_right' | 'none'
export type ScaleStyle = 'bar' | 'ruler'
export type LabelDisplay = 'none' | 'all' | 'first'
export type FilterMode = 'None' | 'Blacklist' | 'Whitelist'

export interface FormState {
  prefix: string
  species: string
  strain: string
  track_type: TrackType
  legend: LegendPosition
  scale_style: ScaleStyle
  show_labels: boolean
  show_labels_linear: LabelDisplay
  separate_strands: boolean
  allow_inner_labels: boolean
  suppress_gc: boolean
  suppress_skew: boolean
  align_center: boolean
  show_gc: boolean
  show_skew: boolean
  normalize_length: boolean
}

export interface AdvancedSettings {
  // Feature types
  features: string[]

  // GC settings
  window_size: number | null
  step_size: number | null
  nt: string

  // Font settings
  def_font_size: number
  label_font_size: number | null

  // Stroke styles
  block_stroke_width: number | null
  block_stroke_color: string | null
  line_stroke_width: number | null
  line_stroke_color: string | null
  axis_stroke_width: number | null
  axis_stroke_color: string | null

  // Legend styles
  legend_box_size: number | null
  legend_font_size: number | null

  // Linear specific
  resolve_overlaps: boolean
  feature_height: number | null
  gc_height: number | null
  comparison_height: number | null
  min_bitscore: number
  evalue: string
  identity: number
  scale_interval: number | null
  scale_font_size: number | null
  scale_stroke_width: number | null
  scale_stroke_color: string | null

  // Circular specific
  outer_label_x_offset: number | null
  outer_label_y_offset: number | null
  inner_label_x_offset: number | null
  inner_label_y_offset: number | null
}

export interface FileState {
  c_gb: File | null
  c_gff: File | null
  c_fasta: File | null
  d_color: File | null
  t_color: File | null
  blacklist: File | null
  whitelist: File | null
  qualifier_priority: File | null
}

export interface LinearSeq {
  gb: File | null
  gff: File | null
  fasta: File | null
  blast: File | null
}

export interface SpecificRule {
  feat: string
  qual: string
  val: string
  color: string
  cap: string
}

export interface PriorityRule {
  feat: string
  order: string
}

export interface ColorState {
  paletteNames: string[]
  selectedPalette: string
  currentColors: Record<string, string>
  specificRules: SpecificRule[]
}

export interface FilterState {
  mode: FilterMode
  manualBlacklist: string
  manualWhitelist: string[]
  priorityRules: PriorityRule[]
}

export function createDefaultFormState(): FormState {
  return {
    prefix: '',
    species: '',
    strain: '',
    track_type: 'tuckin',
    legend: 'right',
    scale_style: 'bar',
    show_labels: false,
    show_labels_linear: 'none',
    separate_strands: true,
    allow_inner_labels: false,
    suppress_gc: false,
    suppress_skew: false,
    align_center: false,
    show_gc: false,
    show_skew: false,
    normalize_length: false,
  }
}

export function createDefaultAdvancedSettings(): AdvancedSettings {
  return {
    features: ['CDS', 'rRNA', 'tRNA', 'tmRNA', 'ncRNA', 'repeat_region'],
    window_size: null,
    step_size: null,
    nt: 'GC',
    def_font_size: 18,
    label_font_size: null,
    block_stroke_width: null,
    block_stroke_color: null,
    line_stroke_width: null,
    line_stroke_color: null,
    axis_stroke_width: null,
    axis_stroke_color: null,
    legend_box_size: null,
    legend_font_size: null,
    resolve_overlaps: false,
    feature_height: null,
    gc_height: null,
    comparison_height: null,
    min_bitscore: 50,
    evalue: '1e-2',
    identity: 0,
    scale_interval: null,
    scale_font_size: null,
    scale_stroke_width: null,
    scale_stroke_color: null,
    outer_label_x_offset: null,
    outer_label_y_offset: null,
    inner_label_x_offset: null,
    inner_label_y_offset: null,
  }
}
