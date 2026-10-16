# Tailwind CSS Play CDN bundle

`tailwindcss-play.js` is the compiler that `tools/generate_web_tailwind_css.mjs`
runs to write `gbdraw/web/css/tailwind.css`. It is a maintenance input: the Web
app does not load it, and it is not deployed or included in the Python package.
`tests/web/tailwind-stylesheet.test.mjs` runs it in CI to check that the tracked
stylesheet is current.

- Source: a snapshot of the unversioned endpoint <https://cdn.tailwindcss.com>
  (the same bytes the Web app loaded from `gbdraw/web/vendor/tailwindcss/`
  before the stylesheet was precompiled).
- Version: the bundle reports `tailwindcss v3.4.17`.
- sha256: `176e894661aa9cdc9a5cba6c720044cbbf7b8bd80d1c9a142a7c24b1b6c50d15`
- A different snapshot can change the generated CSS. Replace the file only
  together with a regenerated `gbdraw/web/css/tailwind.css`.

## License

Tailwind CSS is distributed under the MIT License. The bundle also contains
MIT-licensed runtime dependencies.

MIT License

Copyright (c) Tailwind Labs, Inc.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
