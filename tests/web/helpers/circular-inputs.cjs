// The file input that chooses a Circular input file of a slot ('gb', 'gff' or
// 'fasta') the way a user does: File 1's uploader of that slot once there is a
// card, else the Add files zone, whose pick makes the first cards (one file a
// card; GFF3 + FASTA pairs by file-name stem). A spec that chooses File 1's
// file in one step whether or not the card exists yet uses it.
const circularFileInput = (page, slot = 'gb') => page.locator([
  `[data-circular-source-slot="${slot}"] input[type=file]`,
  '[data-setting=inputs]:not(:has([data-circular-source-card])) [data-circular-file-add] input[type=file]'
].join(', ')).first();

module.exports = { circularFileInput };
