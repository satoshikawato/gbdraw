/**
 * Pyodide service for running Python code in the browser
 *
 * This service manages the Pyodide runtime, including initialization,
 * package installation, and execution of gbdraw Python code.
 */

import type { GbdrawResult } from '@/types'

// Python helper code that will be executed after Pyodide initialization
const PYTHON_HELPERS = `
import json
import glob
import os
import traceback
import tomllib
from importlib import resources
from gbdraw.circular import circular_main
from gbdraw.linear import linear_main

def get_palettes_json():
    try:
        with resources.files("gbdraw.data").joinpath("color_palettes.toml").open("rb") as fh:
            return json.dumps(tomllib.load(fh))
    except:
        return "{}"

def run_gbdraw_wrapper(mode, args):
    for f in glob.glob("*.svg"):
        try:
            os.remove(f)
        except:
            pass

    full_args = args + ["-f", "svg"]
    try:
        if mode == 'circular':
            circular_main(full_args)
        else:
            linear_main(full_args)
    except SystemExit as e:
        if e.code != 0:
            return json.dumps({"error": f"SystemExit: {e}"})
    except Exception:
        return json.dumps({"error": traceback.format_exc()})

    files = glob.glob("*.svg")
    if not files:
        return json.dumps({"error": "No output files generated."})
    results = []
    for fname in sorted(files):
        with open(fname, "r") as f:
            results.append({"name": fname, "svg": f.read(), "success": True})
    return json.dumps(results)

def generate_legend_entry_svg(caption, color, y_offset, rect_size=14, font_size=14, font_family="Arial", x_offset=0):
    """Generate SVG elements for a single legend entry"""
    half = rect_size / 2
    rect_d = f"M 0,{-half} L {rect_size},{-half} L {rect_size},{half} L 0,{half} z"
    rect_svg = f'<path d="{rect_d}" fill="{color}" stroke="none" transform="translate({x_offset}, {y_offset})"/>'

    x_margin = (22 / 14) * rect_size
    text_svg = f'<text font-size="{font_size}" font-family="{font_family}" dominant-baseline="central" text-anchor="start" transform="translate({x_offset + x_margin}, {y_offset})">{caption}</text>'

    return json.dumps({"rect": rect_svg, "text": text_svg})

def extract_features_from_genbank(gb_path):
    """Extract feature info from GenBank file for UI display"""
    import hashlib
    from Bio import SeqIO
    features = []
    record_ids = []
    idx = 0
    try:
        for rec_idx, record in enumerate(SeqIO.parse(gb_path, "genbank")):
            record_id = record.id or f"Record_{rec_idx}"
            record_ids.append(record_id)
            for feat in record.features:
                if feat.type in ['CDS', 'tRNA', 'rRNA', 'ncRNA', 'misc_RNA', 'tmRNA',
                                 'repeat_region', 'misc_feature', 'mobile_element',
                                 'regulatory', 'gene', 'mRNA', 'exon', 'intron']:
                    start = int(feat.location.start)
                    end = int(feat.location.end)
                    strand_raw = feat.location.strand

                    if hasattr(feat.location, 'parts') and feat.location.parts:
                        first_part = feat.location.parts[0]
                        hash_start = int(first_part.start)
                        hash_end = int(first_part.end)
                        hash_strand = first_part.strand
                    else:
                        hash_start = start
                        hash_end = end
                        hash_strand = strand_raw

                    key = f"{feat.type}:{hash_start}:{hash_end}:{hash_strand}"
                    svg_id = "f" + hashlib.md5(key.encode()).hexdigest()[:8]
                    features.append({
                        "id": f"f{idx}",
                        "svg_id": svg_id,
                        "record_id": record_id,
                        "record_idx": rec_idx,
                        "type": feat.type,
                        "start": start,
                        "end": end,
                        "strand": "+" if strand_raw == 1 else "-",
                        "locus_tag": feat.qualifiers.get("locus_tag", [""])[0],
                        "gene": feat.qualifiers.get("gene", [""])[0],
                        "product": feat.qualifiers.get("product", [""])[0],
                        "note": feat.qualifiers.get("note", [""])[0][:50] if feat.qualifiers.get("note") else "",
                    })
                    idx += 1
    except Exception as e:
        return json.dumps({"error": str(e)})
    return json.dumps({"features": features, "record_ids": record_ids})
`

export interface ColorPalettes {
  [key: string]: Record<string, string>
}

export interface LegendEntrySvg {
  rect: string
  text: string
}

export interface ExtractedFeaturesResult {
  features?: Array<{
    id: string
    svg_id: string
    record_id: string
    record_idx: number
    type: string
    start: number
    end: number
    strand: string
    locus_tag: string
    gene: string
    product: string
    note: string
  }>
  record_ids?: string[]
  error?: string
}

class PyodideService {
  private pyodide: PyodideInterface | null = null
  private initialized = false

  /**
   * Check if Pyodide is initialized
   */
  get isReady(): boolean {
    return this.initialized
  }

  /**
   * Initialize Pyodide and install required packages
   */
  async initialize(onProgress: (status: string) => void): Promise<void> {
    if (this.initialized) return

    // Load Pyodide script dynamically
    onProgress('Loading Pyodide runtime...')
    const script = document.createElement('script')
    script.src = 'https://cdn.jsdelivr.net/pyodide/v0.29.0/full/pyodide.js'
    script.crossOrigin = 'anonymous'
    script.integrity = 'sha384-l95tshxQlbjf4kdyWZf10uUL5Dw8/iN9q16SQ+ttOEWA8SN0cLG6BGDGY17GxToh'

    await new Promise<void>((resolve, reject) => {
      script.onload = () => resolve()
      script.onerror = () => reject(new Error('Failed to load Pyodide'))
      document.head.appendChild(script)
    })

    onProgress('Initializing Python environment...')
    this.pyodide = await window.loadPyodide()

    onProgress('Installing micropip...')
    await this.pyodide.loadPackage('micropip')
    const micropip = this.pyodide.pyimport('micropip')

    onProgress('Installing BioPython...')
    await micropip.install('biopython')

    onProgress('Installing svgwrite...')
    await micropip.install('svgwrite')

    onProgress('Installing pandas...')
    await micropip.install('pandas')

    onProgress('Installing fonttools...')
    await micropip.install('fonttools')

    onProgress('Installing bcbio-gff...')
    await micropip.install('bcbio-gff')

    onProgress('Installing gbdraw...')
    await micropip.install('gbdraw-0.8.3-py3-none-any.whl')

    // Initialize helper functions
    onProgress('Initializing helper functions...')
    await this.pyodide.runPythonAsync(PYTHON_HELPERS)

    this.initialized = true
    onProgress('Ready!')
  }

  /**
   * Run gbdraw to generate diagram
   */
  runGbdraw(mode: 'circular' | 'linear', args: string[]): GbdrawResult[] | { error: string } {
    if (!this.pyodide) throw new Error('Pyodide not initialized')

    const runWrapper = this.pyodide.globals.get('run_gbdraw_wrapper') as PyProxy
    const result = runWrapper(mode, this.pyodide.toPy(args.map(String)))
    return JSON.parse(result as string)
  }

  /**
   * Get available color palettes
   */
  getColorPalettes(): ColorPalettes {
    if (!this.pyodide) throw new Error('Pyodide not initialized')

    const result = this.pyodide.runPython('get_palettes_json()')
    return JSON.parse(result as string)
  }

  /**
   * Generate SVG elements for a legend entry
   */
  generateLegendEntry(
    caption: string,
    color: string,
    yOffset: number,
    rectSize = 14,
    fontSize = 14,
    fontFamily = 'Arial',
    xOffset = 0
  ): LegendEntrySvg {
    if (!this.pyodide) throw new Error('Pyodide not initialized')

    // Escape caption and fontFamily for Python string
    const escapedCaption = caption.replace(/\\/g, '\\\\').replace(/"/g, '\\"')
    const escapedFontFamily = fontFamily.replace(/\\/g, '\\\\').replace(/"/g, '\\"')

    const result = this.pyodide.runPython(
      `generate_legend_entry_svg("${escapedCaption}", "${color}", ${yOffset}, ${rectSize}, ${fontSize}, "${escapedFontFamily}", ${xOffset})`
    )
    return JSON.parse(result as string)
  }

  /**
   * Extract features from a GenBank file
   */
  extractFeatures(gbPath: string): ExtractedFeaturesResult {
    if (!this.pyodide) throw new Error('Pyodide not initialized')

    const extractFunc = this.pyodide.globals.get('extract_features_from_genbank') as PyProxy
    const result = extractFunc(gbPath)
    return JSON.parse(result as string)
  }

  /**
   * Write a file to the virtual filesystem
   */
  writeFile(path: string, data: Uint8Array): void {
    if (!this.pyodide) throw new Error('Pyodide not initialized')
    this.pyodide.FS.writeFile(path, data)
  }

  /**
   * Write a File object to the virtual filesystem
   */
  async writeFileFromBlob(file: File, path: string): Promise<boolean> {
    if (!this.pyodide || !file) return false
    const buffer = await file.arrayBuffer()
    this.pyodide.FS.writeFile(path, new Uint8Array(buffer))
    return true
  }

  /**
   * Read a file from the virtual filesystem
   */
  readFile(path: string): Uint8Array {
    if (!this.pyodide) throw new Error('Pyodide not initialized')
    return this.pyodide.FS.readFile(path) as Uint8Array
  }

  /**
   * Delete a file from the virtual filesystem
   */
  deleteFile(path: string): void {
    if (!this.pyodide) throw new Error('Pyodide not initialized')
    try {
      this.pyodide.FS.unlink(path)
    } catch {
      // File might not exist, ignore
    }
  }

  /**
   * Run arbitrary Python code
   */
  runPython(code: string): unknown {
    if (!this.pyodide) throw new Error('Pyodide not initialized')
    return this.pyodide.runPython(code)
  }
}

// Export singleton instance
export const pyodideService = new PyodideService()
