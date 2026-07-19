import { createAnnotationSet, createDefaultAnnotationStyle, normalizeAnnotationSets } from './state.js';
import {
  annotationRecordSelectorValue,
  coordinateTarget,
  featureTarget,
  parseAnnotationRecordSelectorValue
} from './target-actions.js';

const REQUIRED = ['set_id', 'id', 'mark'];
const COLUMNS = [
  'set_id', 'id', 'mark', 'record', 'start', 'end', 'coordinate_space', 'wraps_origin', 'out_of_bounds',
  'feature_selector', 'envelope', 'circular_path', 'label', 'lane', 'legend_label', 'stroke', 'stroke_width',
  'stroke_dasharray', 'line_cap', 'fill', 'fill_opacity', 'hatch_angle', 'hatch_spacing', 'hatch_color',
  'hatch_width', 'hatch_cross', 'label_color', 'label_font_size', 'label_orientation', 'label_position', 'label_offset'
];
const splitTsv = (line) => line.split('\t').map((value) => value.trim());
const boolValue = (value, fallback = false) => {
  if (value == null || value === '') return fallback;
  const text = String(value).toLowerCase();
  if (['true', '1', 'yes', 'on'].includes(text)) return true;
  if (['false', '0', 'no', 'off'].includes(text)) return false;
  throw new Error(`Invalid boolean value: ${value}`);
};
const styleFromRow = (row) => {
  const hasStyle = COLUMNS.slice(15).some((name) => row[name] !== '');
  if (!hasStyle) return null;
  const style = createDefaultAnnotationStyle();
  if (row.stroke) style.stroke = row.stroke;
  if (row.stroke_width) style.strokeWidth = Number(row.stroke_width);
  if (row.stroke_dasharray) style.strokeDasharray = row.stroke_dasharray.split(/[ ,]+/).filter(Boolean).map(Number);
  if (row.line_cap) style.lineCap = row.line_cap;
  if (row.fill) style.fill = row.fill;
  if (row.fill_opacity) style.fillOpacity = Number(row.fill_opacity);
  if (row.label_color) style.labelColor = row.label_color;
  if (row.label_font_size) style.labelFontSize = Number(row.label_font_size);
  if (row.label_orientation) style.labelOrientation = row.label_orientation;
  if (row.label_position) style.labelPosition = row.label_position;
  if (row.label_offset) style.labelOffset = Number(row.label_offset);
  if (['hatch_angle', 'hatch_spacing', 'hatch_color', 'hatch_width', 'hatch_cross'].some((name) => row[name])) {
    style.hatch = {
      angle: Number(row.hatch_angle || 45), spacing: Number(row.hatch_spacing || 5), color: row.hatch_color || '#666666',
      width: Number(row.hatch_width || 1), cross: boolValue(row.hatch_cross)
    };
  }
  return style;
};

export const parseAnnotationTable = (text) => {
  const lines = String(text || '').split(/\r?\n/).filter((line) => line.trim() && !line.trimStart().startsWith('#'));
  if (lines.length === 0) return [];
  const header = splitTsv(lines[0]);
  REQUIRED.forEach((name) => { if (!header.includes(name)) throw new Error(`Annotation table is missing required column: ${name}`); });
  const sets = new Map();
  lines.slice(1).forEach((line, index) => {
    const values = splitTsv(line);
    const row = Object.fromEntries(header.map((name, column) => [name, values[column] ?? '']));
    const rowNumber = index + 2;
    if (!row.set_id || !row.id) throw new Error(`Annotation table row ${rowNumber}: set_id and id are required.`);
    const parsedRecord = parseAnnotationRecordSelectorValue(row.record);
    if (parsedRecord.error) {
      throw new Error(`Annotation table row ${rowNumber}, column 'record': ${parsedRecord.error}`);
    }
    const recordSelector = parsedRecord.selector;
    const record = recordSelector?.kind === 'recordIndex'
      ? { recordId: null, recordIndex: recordSelector.index }
      : { recordId: recordSelector?.value || null, recordIndex: null };
    const hasCoordinates = Boolean(row.start || row.end);
    const hasFeature = Boolean(row.feature_selector);
    if (hasCoordinates === hasFeature) throw new Error(`Annotation table row ${rowNumber}: provide exactly one coordinate or feature target.`);
    const target = hasCoordinates
      ? coordinateTarget({ ...record, start: row.start, end: row.end, coordinateSpace: row.coordinate_space })
      : featureTarget({ ...record, selector: row.feature_selector, extent: row.envelope, circularPath: row.circular_path });
    if (hasCoordinates) {
      target.wrapsOrigin = boolValue(row.wraps_origin, Number(row.start) > Number(row.end));
      target.outOfBounds = row.out_of_bounds || 'clip';
    }
    if (!sets.has(row.set_id)) sets.set(row.set_id, createAnnotationSet({ id: row.set_id }));
    sets.get(row.set_id).annotations.push({
      id: row.id, target, label: row.label || '', mark: row.mark || 'bracket',
      lane: row.lane === '' || row.lane == null ? null : Number(row.lane), style: styleFromRow(row),
      legendLabel: row.legend_label || null, metadata: {}
    });
  });
  return normalizeAnnotationSets(Array.from(sets.values()));
};

export const encodeAnnotationTable = (sets) => {
  const rows = [COLUMNS.join('\t')];
  normalizeAnnotationSets(sets).forEach((set) => set.annotations.forEach((annotation) => {
    const target = annotation.target || {};
    const coordinate = target.kind === 'coordinateSpan';
    const selectors = Array.isArray(target.selectors) ? target.selectors.map((item) => item.key ? `${item.key}=${item.value}` : item.value).join(';') : '';
    const style = annotation.style || set.defaultStyle || {};
    const hatch = style.hatch || {};
    const row = {
      set_id: set.id, id: annotation.id, mark: annotation.mark, record: annotationRecordSelectorValue(target.record),
      start: coordinate ? target.start : '', end: coordinate ? target.end : '', coordinate_space: coordinate ? target.coordinateSpace : '',
      wraps_origin: coordinate && (target.wrapsOrigin || Number(target.start) > Number(target.end)) ? 'true' : '', out_of_bounds: coordinate ? target.outOfBounds : '',
      feature_selector: coordinate ? '' : selectors, envelope: coordinate ? '' : target.envelope, circular_path: coordinate ? '' : target.circularPath,
      label: annotation.label, lane: annotation.lane ?? '', legend_label: annotation.legendLabel ?? set.legendLabel ?? '', stroke: style.stroke ?? '',
      stroke_width: style.strokeWidth ?? '', stroke_dasharray: (style.strokeDasharray || []).join(','), line_cap: style.lineCap ?? '',
      fill: style.fill ?? '', fill_opacity: style.fillOpacity ?? '', hatch_angle: hatch.angle ?? '', hatch_spacing: hatch.spacing ?? '',
      hatch_color: hatch.color ?? '', hatch_width: hatch.width ?? '', hatch_cross: hatch.cross ? 'true' : '', label_color: style.labelColor ?? '',
      label_font_size: style.labelFontSize ?? '', label_orientation: style.labelOrientation ?? '', label_position: style.labelPosition ?? '', label_offset: style.labelOffset ?? ''
    };
    rows.push(COLUMNS.map((name) => String(row[name] ?? '').replace(/[\t\r\n]/g, ' ')).join('\t'));
  }));
  return `${rows.join('\n')}\n`;
};
