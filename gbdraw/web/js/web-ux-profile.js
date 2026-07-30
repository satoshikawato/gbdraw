export const WEB_UX_PROFILE_VERSION = 1;

export const WEB_UX_PROFILE = Object.freeze({
  separateStrands: true,
  circular: Object.freeze({
    singleRecordGrouping: 'single',
    multiRecordGrouping: 'batch',
    gridByDefault: false,
    legend: 'left',
    plotTitlePosition: 'none'
  }),
  linear: Object.freeze({
    legend: 'bottom',
    plotTitlePosition: 'bottom'
  })
});
