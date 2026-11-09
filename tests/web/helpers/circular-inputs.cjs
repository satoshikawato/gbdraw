// The file input that chooses a Circular input file of a slot ('gb', 'gff' or
// 'fasta') the way a user does, by the control's accessible label: File 1's
// uploader of that slot once there is a card, else the Add files zone (a
// button named by its visible label), whose pick makes the first cards (one
// file a card; GFF3 + FASTA pairs by file-name stem). A spec that chooses File
// 1's file in one step whether or not the card exists yet uses it.
// docs/capture/flows/web_capture.py `circular_file_input` is the capture twin.
const SLOT_LABELS = { gb: 'GenBank/DDBJ File', gff: 'GFF3 File', fasta: 'FASTA File' };
const ZONE_NAMES = {
  gb: ['Add GenBank/DDBJ files'],
  gff: ['Add GFF3 and FASTA files'],
  fasta: ['Add GFF3 and FASTA files', 'Add FASTA files']
};

const circularFileInput = (page, slot = 'gb') => {
  const noCard = page.locator('body').filter({
    hasNot: page.getByRole('region', { name: 'Circular input file 1', exact: true, includeHidden: true })
  });
  const zone = ZONE_NAMES[slot]
    .map((name) => noCard.getByRole('button', { name, exact: true, includeHidden: true }))
    .reduce((either, named) => either.or(named));
  return page.getByLabel(SLOT_LABELS[slot], { exact: true }).or(zone.locator('input[type=file]')).first();
};

module.exports = { circularFileInput };
