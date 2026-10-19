// @ts-check
import {
  ChoiceDialog,
  dialogFocus,
  OperationError,
  RecordDisplayControl,
  SpecificRuleRows,
  ActionMenu,
  AutoValueField,
  ColorValueControl,
  DefaultColorScopeLine,
  HelpTip,
  FileUploader,
  GbdrawLoader,
  LinearRecordNames,
  MarkupPreview,
  PreviewZoomPercent,
  previewMarkup,
  previewPanning,
  previewTransform,
  RecordName,
  SettingsSectionIndex,
  SettingNote,
  SettingsSectionTitle,
  settingsSection
} from './components.js';
import { CircularMeasureInput } from './app/circular-track-slots/measure-input.js';
import { DepthTrackIndexInput } from './app/circular-track-slots/track-index-input.js';
import { createAppSetup } from './app/app-setup.js';
import { formatFeatureLocation } from './services/feature-utils.js';
import { t, tn, tx } from './utils/i18n.js';
import { settingApplicationNote, settingApplicationSentence, settingApplicationText } from './services/setting-application.js';

const { createApp } = window.Vue;

const app = createApp({
  components: { ActionMenu, ChoiceDialog, OperationError, SpecificRuleRows, CircularMeasureInput, DepthTrackIndexInput, RecordDisplayControl, AutoValueField, ColorValueControl, DefaultColorScopeLine, FileUploader, GbdrawLoader, HelpTip, LinearRecordNames, MarkupPreview, PreviewZoomPercent, RecordName, SettingNote, SettingsSectionIndex, SettingsSectionTitle },
  methods: { formatFeatureLocation, settingApplicationNote, settingApplicationSentence, settingApplicationText },
  setup: createAppSetup
});
// Every template (the root, the x-templates, the inline component templates)
// shows UI text through these (utils/i18n.js).
Object.assign(app.config.globalProperties, { t, tn, tx });
app.directive('dialog-focus', dialogFocus);
app.directive('settings-section', settingsSection);
app.directive('preview-transform', previewTransform);
app.directive('preview-panning', previewPanning);
app.directive('preview-markup', previewMarkup);
// Vue's production build only logs an error raised in a render, watcher,
// event handler or hook; throw it, as the development build did, so it stays
// an uncaught page error.
app.config.throwUnhandledErrorInProduction = true;

const mountedApp = app.mount('#app');
window.__GBDRAW_APP__ = mountedApp;
