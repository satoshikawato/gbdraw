/// <reference types="vite/client" />

declare module '*.vue' {
  import type { DefineComponent } from 'vue'
  const component: DefineComponent<object, object, unknown>
  export default component
}

// Pyodide types
interface PyodideInterface {
  loadPackage(packages: string | string[]): Promise<void>
  pyimport(module: string): PyProxy
  runPython(code: string): unknown
  runPythonAsync(code: string): Promise<unknown>
  globals: PyProxy
  FS: {
    writeFile(path: string, data: Uint8Array | string): void
    readFile(path: string, opts?: { encoding?: string }): Uint8Array | string
    unlink(path: string): void
    readdir(path: string): string[]
  }
  toPy<T>(obj: T): PyProxy
}

interface PyProxy {
  get(key: string): unknown
  toJs(): unknown
  destroy(): void
  install(pkg: string): Promise<void>
  // Callable interface for Python functions
  (this: void, ...args: unknown[]): unknown
}

declare function loadPyodide(options?: {
  indexURL?: string
  fullStdLib?: boolean
}): Promise<PyodideInterface>

declare global {
  interface Window {
    loadPyodide: typeof loadPyodide
  }
}
