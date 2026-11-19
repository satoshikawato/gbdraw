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
  refuseStrayFileDrop,
  GbdrawLoader,
  LinearRecordNames,
  MarkupPreview,
  NoticeLines,
  PreviewZoomPercent,
  previewMarkup,
  previewPanning,
  previewTransform,
  RecordName,
  resetLayoutState,
  RunBlockedReason,
  SequenceReadLine,
  SettingsSectionIndex,
  SettingNote,
  SettingsSectionTitle,
  settingsSection
} from './components.js';
import { CircularMeasureInput } from './app/circular-track-slots/measure-input.js';
import { DepthTrackIndexInput } from './app/circular-track-slots/track-index-input.js';
import { FileCard } from './app/file-card.js';
import { createAppSetup } from './app/app-setup.js';
import { chooseUiLanguage, showStoredUiLanguage, translationFeedbackUrl } from './app/ui-language.js';
import { headerFit } from './app/header-fit.js';
import { committedValue, settleCommittedField } from './app/committed-value.js';
import { blankNull } from './app/blank-number.js';
import { keepsListFocus } from './app/list-focus.js';
import { titleCutHint } from './app/hint-title.js';
import { formatFeatureLocation } from './services/feature-utils.js';
import { UI_LOCALES, joinSentences, numberText, sentenceSeparator, t, thenArrow, thenText, tm, tn, tx, uiLocale } from './utils/i18n.js';
import { itemListText } from './utils/item-list.js';
import { SELECTOR_TEXT } from './utils/selector-syntax.js';
import { errorText } from './utils/error-text.js';
import { readOtherDigitsInNumberFields } from './utils/number-field-digits.js';
import { INTERACTIVE_GALLERY_URL } from './utils/gallery-routes.js';
import { DEFAULT_RECORD_GAP_PX, ROW_GAPS_MAX_SHARE } from './contracts/schema.generated.js';
import { settingApplicationHelp, settingApplicationSentence } from './services/setting-application.js';
import { comparisonPairLabel } from './services/linear-comparisons.js';

const { createApp } = window.Vue;

const app = createApp({
  components: { ActionMenu, ChoiceDialog, OperationError, SequenceReadLine, SpecificRuleRows, CircularMeasureInput, DepthTrackIndexInput, FileCard, RecordDisplayControl, AutoValueField, ColorSwatch, ColorValueControl, DefaultColorScopeLine, FileUploader, GbdrawLoader, HelpTip, LinearRecordNames, MarkupPreview, NoticeLines, PreviewZoomPercent, RecordName, RunBlockedReason, SettingNote, SettingsSectionIndex, SettingsSectionTitle },
  methods: { chooseUiLanguage, formatFeatureLocation, settingApplicationHelp, settingApplicationSentence, translationFeedbackUrl },
  setup: createAppSetup
});
// Every template (the root, the x-templates, the inline component templates)
// shows UI text through these (utils/i18n.js). `uiLocale` is for a v-memo whose
// body shows UI text: it lists `uiLocale.value` in its dependencies. An error
// text is shown from its parts with `errorText` (utils/error-text.js). Two
// sentences shown together are separated by `sentenceSeparator()` (none in
// Japanese and Chinese); `UI_LOCALES` lists the header Language menu;
// `INTERACTIVE_GALLERY_URL` is the Gallery the header and the footer open.
// `DEFAULT_RECORD_GAP_PX` and `ROW_GAPS_MAX_SHARE` are Python's Auto Linear
// record gap, which the Record gap help states (OV-3572).
// `settleCommittedField` hands an asynchronous commit of a `v-committed-value`
// field to that directive (app/committed-value.js). `thenText` writes "a, then b"
// (`thenArrow` the arrow alone between two elements) and `comparisonPairLabel` a
// comparison pair ("#1 → #2"), in the reading direction of the UI language.
Object.assign(app.config.globalProperties, { t, tn, tx, tm, numberText, itemListText, SELECTOR_TEXT, thenArrow, thenText, comparisonPairLabel, uiLocale, errorText, UI_LOCALES, joinSentences, sentenceSeparator, INTERACTIVE_GALLERY_URL, DEFAULT_RECORD_GAP_PX, ROW_GAPS_MAX_SHARE, settleCommittedField });
// A sentence with elements inside it: <i18n-text :text="t('…{name}…')"> with a #name slot.
app.component('I18nText', I18nText);
app.directive('dialog-focus', dialogFocus);
app.directive('settings-section', settingsSection);
app.directive('preview-transform', previewTransform);
app.directive('preview-panning', previewPanning);
app.directive('reset-layout-state', resetLayoutState);
app.directive('preview-markup', previewMarkup);
app.directive('header-fit', headerFit);
// A field the user types into shows its model's value through `v-committed-value`, not `:value`.
app.directive('committed-value', committedValue);
// A number field bound with `v-model` writes null when it is blank, never '' (app/blank-number.js).
app.directive('blank-null', blankNull);
// A keyed list whose rows move keeps the focus in the row it was in (app/list-focus.js).
app.directive('keeps-list-focus', keepsListFocus);
// Vue's production build only logs an error raised in a render, watcher,
// event handler or hook; throw it, as the development build did, so it stays
// an uncaught page error.
app.config.throwUnhandledErrorInProduction = true;

// The first render shows the language this browser chose; a catalog that
// fails to load leaves English.
await showStoredUiLanguage();
const mountedApp = app.mount('#app');
window.addEventListener('dragover', refuseStrayFileDrop);
window.addEventListener('drop', refuseStrayFileDrop);
// A field's hint cut with an ellipsis shows whole as its title (app/hint-title.js).
window.addEventListener('pointerover', titleCutHint);
window.addEventListener('focusin', titleCutHint);
// A number field reads Persian, Arabic-Indic, Bengali and Devanagari digits typed or pasted
// into it as 0-9 (utils/number-field-digits.js).
window.addEventListener('beforeinput', readOtherDigitsInNumberFields, true);
window.__GBDRAW_APP__ = mountedApp;
