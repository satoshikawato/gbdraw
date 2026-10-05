// One value rule for every TSV table the web app writes as text for Python: a tab or
// line break would add a column or a row, so each run becomes one space, and the cell
// is trimmed (gbdraw/features/source.py trims the same way when it reads).
export const normalizeTsvCell = (value) => String(value ?? '').replace(/[\t\r\n]+/g, ' ').trim();
