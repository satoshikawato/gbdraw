// @ts-check
import {
  ChoiceDialog,
  dialogFocus,
  OperationError,
  RecordDisplayControl,
  SpecificRuleRows,
  ActionMenu,
  AutoValueField,
  ColorSwatch,
  ColorValueControl,
  DefaultColorScopeLine,
  HelpTip,
  I18nText,
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
import { chooseUiLanguage, showStoredUiLanguage } from './app/ui-language.js';
import { formatFeatureLocation } from './services/feature-utils.js';
import { UI_LOCALES, joinSentences, sentenceSeparator, t, tm, tn, tx, uiLocale } from './utils/i18n.js';
import { errorText } from './utils/error-text.js';
import { settingApplicationNote, settingApplicationSentence, settingApplicationText } from './services/setting-application.js';

const { createApp } = window.Vue;

const app = createApp({
  components: { ActionMenu, ChoiceDialog, OperationError, SpecificRuleRows, CircularMeasureInput, DepthTrackIndexInput, RecordDisplayControl, AutoValueField, ColorSwatch, ColorValueControl, DefaultColorScopeLine, FileUploader, GbdrawLoader, HelpTip, LinearRecordNames, MarkupPreview, PreviewZoomPercent, RecordName, SettingNote, SettingsSectionIndex, SettingsSectionTitle },
  methods: { chooseUiLanguage, formatFeatureLocation, settingApplicationNote, settingApplicationSentence, settingApplicationText },
  setup: createAppSetup
});
// Every template (the root, the x-templates, the inline component templates)
// shows UI text through these (utils/i18n.js). `uiLocale` is for a v-memo whose
// body shows UI text: it lists `uiLocale.value` in its dependencies. An error
// text is shown from its parts with `errorText` (utils/error-text.js). Two
// sentences shown together are separated by `sentenceSeparator()` (none in
// Japanese and Chinese); `UI_LOCALES` lists the header Language menu.
Object.assign(app.config.globalProperties, { t, tn, tx, tm, uiLocale, errorText, UI_LOCALES, joinSentences, sentenceSeparator });
// A sentence with elements inside it: <i18n-text :text="t('…{name}…')"> with a #name slot.
app.component('I18nText', I18nText);
app.directive('dialog-focus', dialogFocus);
app.directive('settings-section', settingsSection);
app.directive('preview-transform', previewTransform);
app.directive('preview-panning', previewPanning);
app.directive('preview-markup', previewMarkup);
// Vue's production build only logs an error raised in a render, watcher,
// event handler or hook; throw it, as the development build did, so it stays
// an uncaught page error.
app.config.throwUnhandledErrorInProduction = true;

// The first render shows the language this browser chose; a catalog that
// fails to load leaves English.
await showStoredUiLanguage();
const mountedApp = app.mount('#app');
window.__GBDRAW_APP__ = mountedApp;
