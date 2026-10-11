// @ts-check
export const WEB_UX_PROFILE_VERSION = 1;

export const WEB_UX_PROFILE = Object.freeze({
  separateStrands: true,
  circular: Object.freeze({
    singleRecordGrouping: 'single',
    multiRecordGrouping: 'batch',
    gridByDefault: true,
    legend: 'left',
    plotTitlePosition: 'none'
  }),
  // A file with more records opens its record list after an upload (D-04).
  recordList: Object.freeze({ autoOpenAbove: 20 }),
  linear: Object.freeze({
    arrangeInRowsByDefault: true,
    legend: 'bottom',
    plotTitlePosition: 'bottom'
  })
});
