import { HelpTip, FileUploader } from './components.js';
import { createAppSetup } from './app/app-setup.js';

const { createApp } = window.Vue;

createApp({
  components: { FileUploader, HelpTip },
  setup: createAppSetup
}).mount('#app');
