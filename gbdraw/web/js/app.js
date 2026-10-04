import {
  OperationError,
  RecordDisplayControl,
  AutoValueField,
  ColorValueControl,
  HelpTip,
  FileUploader
} from './components.js';
import { CircularMeasureInput } from './app/circular-track-slots/measure-input.js';
import { createAppSetup } from './app/app-setup.js';
import { formatFeatureLocation } from './app/feature-utils.js';

const { createApp } = window.Vue;

const app = createApp({
  components: { OperationError, CircularMeasureInput, RecordDisplayControl, AutoValueField, ColorValueControl, FileUploader, HelpTip },
  methods: { formatFeatureLocation },
  setup: createAppSetup
});

const mountedApp = app.mount('#app');
window.__GBDRAW_APP__ = mountedApp;
