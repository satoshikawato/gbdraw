export const setDpiInPng = async (blob, dpi) => {
  const pixelsPerMeter = Math.round(dpi / 0.0254);
  const buffer = await blob.arrayBuffer();
  const view = new DataView(buffer);
  const uint8 = new Uint8Array(buffer);
  let offset = 8;
  let physChunk = null;
  let idatOffset = null;
  while (offset < view.byteLength) {
    const length = view.getUint32(offset);
    const type = String.fromCharCode(
      uint8[offset + 4],
      uint8[offset + 5],
      uint8[offset + 6],
      uint8[offset + 7]
    );
    if (type === 'pHYs') {
      physChunk = offset;
      break;
    }
    if (type === 'IDAT') {
      idatOffset = offset;
      break;
    }
    offset += 12 + length;
  }
  const ppm = pixelsPerMeter;
  const physData = new Uint8Array(9);
  const dv = new DataView(physData.buffer);
  dv.setUint32(0, ppm);
  dv.setUint32(4, ppm);
  physData[8] = 1;
  const crcTable = [];
  for (let n = 0; n < 256; n++) {
    let c = n;
    for (let k = 0; k < 8; k++) {
      c = (c & 1) ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
    }
    crcTable[n] = c;
  }
  const crc32 = (buf) => {
    let crc = 0 ^ (-1);
    for (let i = 0; i < buf.length; i++) {
      crc = (crc >>> 8) ^ crcTable[(crc ^ buf[i]) & 0xff];
    }
    return (crc ^ (-1)) >>> 0;
  };
  const newChunk = new Uint8Array(21);
  const newDv = new DataView(newChunk.buffer);
  newDv.setUint32(0, 9);
  newChunk.set([112, 72, 89, 115], 4);
  newChunk.set(physData, 8);
  const crcInput = new Uint8Array(13);
  crcInput.set([112, 72, 89, 115], 0);
  crcInput.set(physData, 4);
  newDv.setUint32(17, crc32(crcInput));
  if (physChunk !== null) {
    const newBuffer = new Uint8Array(buffer.byteLength);
    newBuffer.set(uint8.slice(0, physChunk), 0);
    newBuffer.set(newChunk, physChunk);
    const oldLen = view.getUint32(physChunk);
    newBuffer.set(uint8.slice(physChunk + 12 + oldLen), physChunk + 21);
    return new Blob([newBuffer], { type: 'image/png' });
  }
  const newBuffer = new Uint8Array(buffer.byteLength + 21);
  newBuffer.set(uint8.slice(0, idatOffset), 0);
  newBuffer.set(newChunk, idatOffset);
  newBuffer.set(uint8.slice(idatOffset), idatOffset + 21);
  return new Blob([newBuffer], { type: 'image/png' });
};
