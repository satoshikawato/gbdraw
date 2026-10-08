// @ts-check
import {
  ChoiceDialog,
  dialogFocus,
  OperationError,
  RecordDisplayControl,
  AutoValueField,
  ColorValueControl,
  HelpTip,
  FileUploader
} from './components.js';
import { CircularMeasureInput } from './app/circular-track-slots/measure-input.js';
import { DepthTrackIndexInput } from './app/circular-track-slots/track-index-input.js';
import { createAppSetup } from './app/app-setup.js';
import { formatFeatureLocation } from './services/feature-utils.js';

const { createApp } = window.Vue;

const app = createApp({
  components: { ChoiceDialog, OperationError, CircularMeasureInput, DepthTrackIndexInput, RecordDisplayControl, AutoValueField, ColorValueControl, FileUploader, HelpTip },
  methods: { formatFeatureLocation },
  setup: createAppSetup
});
app.directive('dialog-focus', dialogFocus);

const mountedApp = app.mount('#app');
window.__GBDRAW_APP__ = mountedApp;
