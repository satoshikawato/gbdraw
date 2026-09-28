import json,os,socket,subprocess
with socket.socket() as server:
    server.bind(('127.0.0.1',0))
    port=server.getsockname()[1]
env=dict(os.environ,GBDRAW_WEB_TEST_PORT=str(port))
command=['npx','playwright','test','--config=playwright.config.js','--project=chromium','--workers=1','--retries=0','--reporter=list','--output=/tmp/issue597-S05-trusted-guard-integration-20260927-evidence/browser-results','tests/web/session-import-worker.playwright.spec.js','tests/web/session-operation-consistency.playwright.spec.js','--grep=real import Worker|unsafe keys, malformed|teardown cancels|occupies all semantic owners']
print(json.dumps({'command':command,'environment':{'GBDRAW_WEB_TEST_PORT':str(port)}},indent=2),flush=True)
raise SystemExit(subprocess.call(command,env=env))
