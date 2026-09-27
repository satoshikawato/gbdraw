import {
  RecordDisplayControl,
  AutoValueField,
  ColorValueControl,
  HelpTip,
  FileUploader
} from './components.js';
import { CircularMeasureInput } from './app/circular-track-slots/measure-input.js';
import { createAppSetup } from './app/app-setup.js';

const { createApp } = window.Vue;

const app = createApp({
  components: { CircularMeasureInput, RecordDisplayControl, AutoValueField, ColorValueControl, FileUploader, HelpTip },
  setup: createAppSetup
});

const mountedApp = app.mount('#app');
window.__GBDRAW_APP__ = mountedApp;
