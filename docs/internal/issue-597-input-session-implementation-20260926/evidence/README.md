# 基準動作と今後の測定

既存観察のbaseは `d457b7189b137185a8dec800819a312c30b969fa`。
[current-browser.json](./current-browser.json) は1440×1000/Chromium149の限定的observation。
[probe_current_behavior.py](./probe_current_behavior.py) はnative upload→saved preview→native replacementを再実行する。
scriptをserverと同じcheckout/inputsへ向け、local serverはtask専用portを使う。

```bash
python3 -m http.server 4597 --bind 127.0.0.1
```

別terminalから:

```bash
python docs/internal/issue-597-input-session-implementation-20260926/evidence/probe_current_behavior.py \
  --url http://127.0.0.1:4597 --output /tmp/issue597-current-browser.json
```

このprobeはGenerate/refreshを直接呼ばない。raw dataはnative2 records/spinbuttons2/Workers0、
saved preview records0/status loading/Workers0、replacement ready1、page errors0。
S03以降はstatus仕様が変更されるため、基準観察を改ざんせず新contractsで比較する。
過去の30 Node checksはsession-file/session-request/session-active-files/record-display-options/record-metadata-inference。
real large Save/Load、transport/heap、GFF/DDBJ variantsはS01以降で実測する。
このevidence directoryを新しいauthorityまたは第二のruntime fixture ownerにしない。
