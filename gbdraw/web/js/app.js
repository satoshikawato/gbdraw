import { HelpTip, FileUploader } from './components.js';
import { createAppSetup } from './app/app-setup.js';

const { createApp } = window.Vue;

const app = createApp({
  components: { FileUploader, HelpTip },
  setup: createAppSetup
});

const mountedApp = app.mount('#app');
window.__GBDRAW_APP__ = mountedApp;
