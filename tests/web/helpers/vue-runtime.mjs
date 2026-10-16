// The Vue build that gbdraw/web/index.html loads, evaluated in a fresh node:vm
// context. Node tests that need Vue's reactivity or template compiler load it
// from here, so the page stays the one owner of which build that is.
import { readFileSync } from 'node:fs';
import vm from 'node:vm';

const repositoryFile = (path) => readFileSync(new URL(`../../../${path}`, import.meta.url), 'utf8');

const pageScript = repositoryFile('gbdraw/web/index.html').match(/<script src="\.\/(vendor\/vue\/[^"]+)"><\/script>/);
if (!pageScript) throw new Error('gbdraw/web/index.html loads no vendor/vue build');
export const PAGE_VUE_PATH = `gbdraw/web/${pageScript[1]}`;
const source = repositoryFile(PAGE_VUE_PATH);

export const loadPageVue = () => {
  const context = vm.createContext({ console });
  vm.runInContext(source, context);
  return context.Vue;
};
