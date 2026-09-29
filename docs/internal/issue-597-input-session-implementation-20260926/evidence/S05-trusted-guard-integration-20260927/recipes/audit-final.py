import hashlib,json,pathlib,subprocess
root=pathlib.Path('/tmp/issue597-S05-trusted-guard-integration-20260927-evidence')
git=lambda *args:subprocess.check_output(['git',*args])
sha=lambda b:hashlib.sha256(b).hexdigest()
oldBase='494091aa68ca59ffa27ecaa6c3df19da4fbf5090'
newBase='98c21116f6439d5721e7ea62ae1a21a7cf2d4319'
first='e6f3227aba4ee3f7c911cfe08fbc3528c5d6dca8'
head=git('rev-parse','HEAD').decode().strip()
assert git('rev-parse','origin/dev').decode().strip()==newBase
assert git('rev-list','--parents','-n','1',head).decode().split()[1:]==[first,newBase]
subprocess.run(['git','merge-base','--is-ancestor',oldBase,newBase],check=True)
assert not git('status','--porcelain=v1','--untracked-files=all')
paths=git('diff','--name-only',oldBase,newBase).decode().splitlines()
assert paths==['docs/internal/issue-601-pr-smoke-inventory-20260927.md','tests/ci/playwright-inventory.test.mjs']
assert git('diff',first,head)==git('diff',oldBase,newBase)
assert not git('diff','--name-only',first,head,'--','gbdraw','tools','tests/web','.github')
initial=json.loads((root/'fingerprints.json').read_text())
assert not (root/'fingerprints-initial.json').exists()
(root/'fingerprints-initial.json').write_bytes((root/'fingerprints.json').read_bytes())
initial['initialTrustedDev']=oldBase
initial['trustedDev']=newBase
initial['initialIntegrationTree']=initial['integratedTree']
initial['integratedTree']=git('rev-parse',head+'^{tree}').decode().strip()
initial['integrationCommits']=[first,head]
for p,values in initial['protectedFilesEqualTrustedDev'].items():
    values['trustedDev']=sha(git('show',newBase+':'+p))
    values['integrated']=sha(git('show',head+':'+p))
    assert values['trustedDev']==values['integrated']
for p,expected in initial['source']['files'].items(): assert sha(git('show',head+':'+p))==expected
for p,expected in initial['pinnedEvidenceUnchanged'].items(): assert sha(git('show',head+':'+p))==expected
initial['lateDevTransition']={'from':oldBase,'to':newBase,'paths':paths,'secondIntegrationCommit':head,'exactUpstreamDiffImported':True,'runtimeCheckerGuardMappedSourceUnchanged':True,'testOnlyInventoryUpperBound':'PR #625 changes the CI-only PR-smoke case ceiling from 13 to 19. Imported unchanged from trusted dev, independently authored before this session; no S05 timeout or performance budget changes.','verificationReuse':'139 architecture PASS (including all four named contracts), 790 Node PASS, five Chromium PASS and three negative inventory mutations retain identical source/tests/checker/input/environment/acceptance after this two-path transition. Latest-base policy rerun and affected CI contracts rerun.'}
(root/'fingerprints.json').write_text(json.dumps(initial,indent=2)+'\n')
print(json.dumps({'finalIntegrationCommit':head,'parents':[first,newBase],'lateDevTransition':initial['lateDevTransition'],'protectedFilesIdentical':len(initial['protectedFilesEqualTrustedDev']),'sourceManifestSha256':initial['source']['manifestSha256'],'pinnedEvidenceUnchanged':len(initial['pinnedEvidenceUnchanged'])},indent=2))
