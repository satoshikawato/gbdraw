import hashlib,json,pathlib,subprocess,tarfile,tempfile
root=pathlib.Path('/tmp/issue597-S05-trusted-guard-integration-20260927-evidence')
overlay=pathlib.Path(tempfile.mkdtemp(prefix='issue597-integrated-guard-negative-',dir='/tmp'))
tree=subprocess.check_output(['git','write-tree'],text=True).strip()
paths=['.github','tools','tests/web','docs/internal','gbdraw/web/js','gbdraw/web/index.html','package.json','playwright.config.js','playwright.functional.config.js','playwright.pr-smoke.config.js']
archive=subprocess.Popen(['git','archive',tree,*paths],stdout=subprocess.PIPE)
with tarfile.open(fileobj=archive.stdout,mode='r|') as tar: tar.extractall(overlay,filter='data')
assert archive.wait()==0
client=overlay/'gbdraw/web/js/services/session-import-client.js'
original=client.read_bytes()
assert original.count(b'new Worker(')==1
pattern='Worker construction and the diagram-generation client have explicit owners|shared privileged detectors preserve the characterized current-source facts'
results=[]
mutations=[('extra-constructor',original+b'\nnew Worker("./extra.js");\n',None),('missing-constructor',original.replace(b'new Worker(',b'createWorker('),None),('extra-owner',original,overlay/'gbdraw/web/js/services/session-import-client-extra.js')]
for name,body,extra in mutations:
    client.write_bytes(body)
    if extra: extra.write_text('new Worker("./extra.js");\n')
    command=['node','--test','--test-name-pattern='+pattern,'tests/web/architecture-contracts.test.mjs']
    raw=root/('negative-corrected-'+name+'.log')
    assert not raw.exists()
    with raw.open('wb') as stream: result=subprocess.run(command,cwd=overlay,stdout=stream,stderr=subprocess.STDOUT)
    output=raw.read_text()
    assert result.returncode==1, (name,result.returncode)
    assert 'fail 2' in output and 'pass 0' in output, (name,output)
    assert 'Worker construction and the diagram-generation client have explicit owners' in output
    assert 'shared privileged detectors preserve the characterized current-source facts' in output
    results.append({'mutation':name,'command':command,'cwd':str(overlay),'exit':result.returncode,'expectedExit':1,'summary':{'pass':0,'fail':2},'log':str(raw),'sha256':hashlib.sha256(raw.read_bytes()).hexdigest()})
    if extra: extra.unlink()
    client.write_bytes(original)
print(json.dumps({'sourceTree':tree,'overlay':str(overlay),'clientOriginalSha256':hashlib.sha256(original).hexdigest(),'acceptanceUnchanged':True,'invocations':results},indent=2))
