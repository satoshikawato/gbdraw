/** @type {import('tailwindcss').Config} */
module.exports = {
  content: [
    './index.html',
    './js/**/*.js'
  ],
  safelist: [
    // Dynamic Vue class bindings in index.html (:class / ternary / object forms)
    'bg-white',
    'text-blue-600',
    'shadow',
    'ring-1',
    'ring-black/5',
    'text-slate-500',
    'hover:text-slate-700',
    'bg-blue-600',
    'text-white',
    'text-slate-600',
    'border',
    'bg-blue-50',
    'border-blue-300',
    'translate-x-0',
    '-translate-x-80',
    'translate-x-full',
    'w-7',
    'h-5',
    'rounded',
    'shrink-0',
    'cursor-pointer',
    'cursor-not-allowed',
    'opacity-50',
    'mb-0',
    'py-0.5',
    'min-h-[30px]',
    'text-base',
    'text-xl',
    'text-[10px]',
    'text-xs',
    'font-bold'
  ]
};
