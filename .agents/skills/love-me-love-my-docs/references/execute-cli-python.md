# Executable CLI and Python documentation

Use the published command or code as the executable authority.

## Clean-directory contract

For each scenario:

1. create a new temporary directory;
2. acquire original sequence inputs from the documented authoritative database
   by accession; offline automation may copy only checksum-verified mirrors of
   those exact downloads;
3. create only the files whose complete contents appear in the executable
   source named by the scenario;
4. run the exact fenced CLI block or literal Python program without editing;
5. capture exit status and safe stdout/stderr;
6. reject undeclared inputs and unexpected outputs;
7. validate every promised filename, format, and semantic result;
8. compare the generated figure with the published artifact at readable scale.

Do not make repository-bundled sequence files or finished sessions part of the
reader workflow. Do not rely on the repository root, shell aliases, previous
outputs, undeclared network access, or test-only helper code that the reader
cannot see.

## Evidence by surface

- For CLI, execute the exact commands as subprocesses. Test expected failures
  and overwrite behavior explicitly. Prefer the generated figure and parsed
  artifact evidence to a terminal screenshot.
- For Python, extract or import the literal code block. Assert public return
  types as well as saved output. A parallel sample in a test file does not
  prove the published snippet.
- Parse SVG/XML and structured outputs for stable biological or domain
  semantics. Validate raster/vector signatures and dimensions where relevant.
- Keep generated reference artifacts owned by the recipe. Regenerate and
  visually inspect them; never patch their bytes by hand.

When the scenario supports a public procedure, document the empty-directory
preparation, direct input links, exact save names, working-tree layout, run
command, expected console message, generated files, visible checks, and common
failures. Internal regression evidence may retain only the literal inputs and
executable block needed by its runner.
