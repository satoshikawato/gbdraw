/**
 * PNG DPI metadata injection utility
 *
 * This utility injects pHYs (physical dimensions) chunk into PNG files
 * to set the DPI metadata. This is used for high-resolution PNG exports.
 */

/**
 * Calculate CRC32 for PNG chunk data
 */
function createCrc32Table(): number[] {
  const table: number[] = []
  for (let n = 0; n < 256; n++) {
    let c = n
    for (let k = 0; k < 8; k++) {
      c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1
    }
    table[n] = c
  }
  return table
}

function crc32(buf: Uint8Array, table: number[]): number {
  let crc = 0 ^ -1
  for (let i = 0; i < buf.length; i++) {
    crc = (crc >>> 8) ^ table[(crc ^ buf[i]) & 0xff]
  }
  return (crc ^ -1) >>> 0
}

/**
 * Set DPI metadata in a PNG blob
 *
 * Injects or replaces the pHYs chunk in a PNG file to set the
 * pixels per meter value, which determines the DPI of the image.
 *
 * @param blob - The PNG blob to modify
 * @param dpi - The desired DPI value (e.g., 300 for print quality)
 * @returns A new blob with the DPI metadata set
 */
export async function setDpiInPng(blob: Blob, dpi: number): Promise<Blob> {
  const pixelsPerMeter = Math.round(dpi / 0.0254)
  const buffer = await blob.arrayBuffer()
  const view = new DataView(buffer)
  const uint8 = new Uint8Array(buffer)

  // Find pHYs or IDAT chunk
  let offset = 8 // Skip PNG signature
  let physChunk: number | null = null
  let idatOffset: number | null = null

  while (offset < view.byteLength) {
    const length = view.getUint32(offset)
    const type = String.fromCharCode(uint8[offset + 4], uint8[offset + 5], uint8[offset + 6], uint8[offset + 7])

    if (type === 'pHYs') {
      physChunk = offset
      break
    }
    if (type === 'IDAT') {
      idatOffset = offset
      break
    }
    offset += 12 + length
  }

  // Create pHYs chunk data (9 bytes)
  const ppm = pixelsPerMeter
  const physData = new Uint8Array(9)
  const dv = new DataView(physData.buffer)
  dv.setUint32(0, ppm) // X pixels per unit
  dv.setUint32(4, ppm) // Y pixels per unit
  physData[8] = 1 // Unit is meter

  // Calculate CRC32
  const crcTable = createCrc32Table()
  const crcInput = new Uint8Array(13)
  crcInput.set([112, 72, 89, 115], 0) // "pHYs" in bytes
  crcInput.set(physData, 4)
  const checksum = crc32(crcInput, crcTable)

  // Create new chunk (21 bytes: 4 length + 4 type + 9 data + 4 CRC)
  const newChunk = new Uint8Array(21)
  const newDv = new DataView(newChunk.buffer)
  newDv.setUint32(0, 9) // Length
  newChunk.set([112, 72, 89, 115], 4) // "pHYs"
  newChunk.set(physData, 8)
  newDv.setUint32(17, checksum)

  if (physChunk !== null) {
    // Replace existing pHYs chunk
    const oldLen = view.getUint32(physChunk)
    const newBuffer = new Uint8Array(buffer.byteLength)
    newBuffer.set(uint8.slice(0, physChunk), 0)
    newBuffer.set(newChunk, physChunk)
    newBuffer.set(uint8.slice(physChunk + 12 + oldLen), physChunk + 21)
    return new Blob([newBuffer], { type: 'image/png' })
  } else if (idatOffset !== null) {
    // Insert new pHYs chunk before IDAT
    const newBuffer = new Uint8Array(buffer.byteLength + 21)
    newBuffer.set(uint8.slice(0, idatOffset), 0)
    newBuffer.set(newChunk, idatOffset)
    newBuffer.set(uint8.slice(idatOffset), idatOffset + 21)
    return new Blob([newBuffer], { type: 'image/png' })
  }

  // If we can't find a place to insert, return original blob
  return blob
}
