import hashlib,json,pathlib,re,subprocess
START='c3ae518b914024d1299c1f3f6487640c8a2fbaaf'
BASE='494091aa68ca59ffa27ecaa6c3df19da4fbf5090'
AUTO_TREE='5782d74f991d5d7efa971981a00b6952e975ce4e'
ROOT=pathlib.Path('/tmp/issue597-S05-trusted-guard-integration-20260927-evidence')
PREFIX='docs/internal/issue-597-input-session-implementation-20260926/'
git=lambda *args:subprocess.check_output(['git',*args])
sha=lambda b:hashlib.sha256(b).hexdigest()
blob=lambda ref,path:git('show',ref+':'+path)
def tracked(ref,*paths):return git('ls-tree','-r','--name-only',ref,*paths).decode().splitlines()
tree=git('write-tree').decode().strip()
assert git('rev-parse','origin/dev').decode().strip()==BASE
assert git('diff','--name-only',AUTO_TREE,tree).decode().splitlines()==['gbdraw/web/index.html']
assert not git('diff','--name-only')
protected=set(json.loads(pathlib.Path('docs/internal/issue-597-session-import-guard-20260927/fingerprints.json').read_text())['unchangedCheckerPolicyAuthorityWorkflowFiles'])
protected.update(json.loads(pathlib.Path(PREFIX+'evidence/S05-resume-fingerprints.json').read_text())['protectedActiveFilesEqualTrustedDev'])
protected.add('tests/web/architecture-contracts.test.mjs')
identities={p:{'trustedDev':sha(blob(BASE,p)),'integrated':sha(blob(tree,p))} for p in sorted(protected)}
assert all(v['trustedDev']==v['integrated'] for v in identities.values())
arch='tests/web/architecture-contracts.test.mjs'
names=['production rendering crosses the canonical request and Worker boundary','History intent and SVG admission have one production ownership path']
def body(text,name):
    start=text.index("test('"+name+"'")
    end=text.find('\ntest(',start+1)
    return text[start:end if end!=-1 else len(text)].encode()
bodies={name:{ref:sha(body(blob(rev,arch).decode(),name)) for ref,rev in [('start',START),('trustedDev',BASE),('integrated',tree)]} for name in names}
assert all(len(set(v.values()))==1 for v in bodies.values())
map_path='tools/web-product-impact-map.json'
maps={ref:json.loads(blob(rev,map_path)) for ref,rev in [('start',START),('trustedDev',BASE),('integrated',tree)]}
assert maps['integrated']==maps['trustedDev']==maps['start']
refs=[]
for concern in maps['trustedDev']['concerns']:
 for contract in concern['contracts']:
  if any(contract['ref'].endswith('::'+name) for name in names): refs.append({'concern':concern['key'],**contract})
assert len(refs)==2
client='gbdraw/web/js/services/session-import-client.js'
worker='gbdraw/web/js/workers/session-import-worker.js'
unchanged_transport={p:{'start':sha(blob(START,p)),'integrated':sha(blob(tree,p))} for p in [client,worker]}
assert all(v['start']==v['integrated'] for v in unchanged_transport.values())
constructors={p.removeprefix('gbdraw/web/js/'):len(re.findall(rb'\bnew\s+Worker\s*\(',blob(tree,p))) for p in tracked(tree,'gbdraw/web/js') if p.endswith('.js')}
constructors={p:n for p,n in constructors.items() if n}
assert constructors=={'services/diagram-generation.js':1,'services/losat.js':2,'services/session-import-client.js':1,'workers/losat-threaded-worker.js':2}
pinned={p:sha(blob(START,p)) for p in tracked(START,PREFIX)}
assert all(sha(blob(tree,p))==h for p,h in pinned.items())
contract='docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md'
def section(ref,name):
 text=blob(ref,contract).decode()
 start=text.index('### '+name+':')
 end=text.find('\n### ',start+1)
 if end==-1:end=text.index('\n## Acceptance contract catalog',start)
 return sha(text[start:end].encode())
receipts={name:{label:section(ref,name) for label,ref in [('start',START),('trustedDev',BASE),('integrated',tree)]} for name in ['PD-OI-044','PD-OI-045']}
assert all(len(set(v.values()))==1 for v in receipts.values())
source={p:sha(blob(tree,p)) for p in tracked(tree,'gbdraw')}
old=json.loads(pathlib.Path(PREFIX+'evidence/S05-resume-fingerprints.json').read_text())
python_inputs={p:h for p,h in old['runtime']['sha256'].items() if not p.startswith('gbdraw/web/')}
python_unchanged=all(source.get(p)==h for p,h in python_inputs.items())
assert python_unchanged
inputs={i['path']:i['sha256'] for i in old['fixture']['sources']}
assert all(sha(pathlib.Path(p).read_bytes())==h for p,h in inputs.items())
fixtures=['tests/web/architecture-ratchet-fixtures.test.mjs','tests/web/product-impact-ratchet-fixtures.test.mjs','tests/web/promotion-readiness.test.mjs','tests/web/web-promotion-context.test.mjs']
fixtures.extend(p for p in protected if p.startswith('tools/'))
fixture_match={p:sha(blob(tree,p)) for p in sorted(set(fixtures))}
prior_guard=json.loads(pathlib.Path('docs/internal/issue-597-session-import-guard-20260927/fingerprints.json').read_text())
assert all(fixture_match[p]==sha(blob(prior_guard['base']['revision'],p)) for p in fixture_match)
prior_log=pathlib.Path('docs/internal/issue-597-session-import-guard-20260927/logs/guard-fixtures.log')
prior_validation=json.loads(pathlib.Path('docs/internal/issue-597-session-import-guard-20260927/validation.json').read_text())
fixture_invocation=next(i for i in prior_validation['invocations'] if i['name']=='guard-fixtures')
assert sha(prior_log.read_bytes())==fixture_invocation['sha256']
assert fixture_invocation['summary']=={'tests':78,'pass':78,'fail':0}
diff_reviews={category:git('diff','--name-status',START,tree,'--',*paths).decode().splitlines() for category,paths in {
 'production':['gbdraw'], 'tests':['tests'], 'docsEvidence':['docs'],
 'generated':['dist','gbdraw.egg-info','tests/reference_outputs','gbdraw/web/gallery','examples/gbdraw_social_preview.png']}.items()}
assert not diff_reviews['generated']
for category,paths in {'production':['gbdraw'],'tests':['tests'],'docs-evidence':['docs'],'generated':['dist','gbdraw.egg-info','tests/reference_outputs','gbdraw/web/gallery','examples/gbdraw_social_preview.png']}.items():
 (ROOT/('review-'+category+'.diff')).write_bytes(git('diff',START,tree,'--',*paths))
result={'startHead':START,'trustedDev':BASE,'integratedTree':tree,
 'automaticMergeTree':AUTO_TREE,'manualResolutionPaths':['gbdraw/web/index.html'],
 'noIndependentRuntimeDesign':True,'protectedFilesEqualTrustedDev':identities,
 'mappedBodies':bodies,'mappedAuthorityReferences':refs,'completeProductMapIdentical':True,
 'unchangedImportTransport':unchanged_transport,'workerConstructorInventory':constructors,
 'issue597DecisionSections':receipts,'contractRevision':24,'pinnedEvidenceUnchanged':pinned,
 'source':{'files':source,'manifestSha256':sha(json.dumps(source,sort_keys=True,separators=(',',':')).encode())},
 'unchangedPythonSourceAndData':{'count':len(python_inputs),'identical':python_unchanged},
 'unchangedRealFixtureInputs':inputs,
 'reusedGuardFixtures':{'invocation':fixture_invocation,'sourceFingerprints':fixture_match,'conditions':'Same exact four fixture test bodies, unchanged checker/detector/evaluator/parsers/policy/rules/map/decisions, same Node v26.8.2 and OS; fixture-owned input and acceptance unchanged. No production source execution in these pure/mechanical fixtures. 78 PASS remains scoped to these mechanics.'},
 'priorRuntimeMeasurements':'Historical only: integration changes 14 Web runtime files and the aggregate runtime fingerprint. No old browser/performance PASS is claimed for the new integrated runtime. Old FAIL/null invocations remain unchanged.',
 'reviewedDiffCategories':diff_reviews}
(ROOT/'fingerprints.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['source','pinnedEvidenceUnchanged','protectedFilesEqualTrustedDev','reviewedDiffCategories','reusedGuardFixtures']},indent=2))
print('Protected identity count:',len(identities),'Pinned evidence count:',len(pinned),'Source count:',len(source),'Reused fixture count: 78')
