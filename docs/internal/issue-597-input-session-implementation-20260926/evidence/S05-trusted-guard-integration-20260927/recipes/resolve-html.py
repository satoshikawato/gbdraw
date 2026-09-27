import hashlib,json,pathlib,re,subprocess
path=pathlib.Path('gbdraw/web/index.html')
s=path.read_text()
pattern=re.compile(r'^<<<<<<< HEAD\n(.*?)^=======\n(.*?)^>>>>>>> origin/dev\n',re.M|re.S)
blocks=list(pattern.finditer(s))
assert len(blocks)==6
resolved=[]
records=[]
for number,b in enumerate(blocks,1):
    ours,theirs=b[1],b[2]
    if number==1:
        assert 'saveSessionWithTitle' in ours and 'session-save-help' in theirs
        result=theirs.replace('<button @click=', '<button data-history-ignore @click=',1).replace(':disabled="sessionSavePending"', ':disabled="!sessionSaveAvailable"')
        reason='Preserve S05 History exclusion and canonical Save availability; import dev Save explanation and help reference.'
    elif number==2:
        assert not ours and 'refreshCircularRecordOrder' in theirs
        result=''.join(line for line in theirs.splitlines(keepends=True) if '<button' not in line)
        reason='Import dev rotation timing explanation; preserve S03 explicit inspection by keeping removed lazy rotation-load button absent.'
    elif number==3:
        assert 'runAnalysis' in ours and 'generate-help' in theirs
        result=theirs.replace(':disabled="processing"', ':disabled="!semanticMutationAvailable || (processing)"')
        reason='Preserve S04 canonical operation availability; import dev generation feedback, live announcement and help.'
    elif number==4:
        assert 'keep_definition_left_aligned' in ours and 'linear-definition-lock-help' in theirs
        result=theirs.replace('aria-describedby="linear-definition-lock-help">','aria-describedby="linear-definition-lock-help" :disabled="!semanticMutationAvailable">')
        reason='Preserve S04 mutation availability; import dev responsive label, accessible name and help.'
    elif number==5:
        assert 'row.visibilityType' in ours and 'linear-label-visibility-' in theirs
        result=theirs.replace(':aria-label="row.label + \' visibility\'">', ':aria-label="row.label + \' visibility\'" :disabled="!semanticMutationAvailable">')
        reason='Preserve S04 mutation availability; import dev label/control association.'
    else:
        assert 'selectedPalette' in ours and 'data-palette-application-help' in theirs
        result=theirs.replace('class="form-input form-input-compact mb-2">','class="form-input form-input-compact mb-2" :disabled="!semanticMutationAvailable">')
        reason='Preserve S04 mutation availability; import dev palette application explanation.'
    assert result!=theirs or number==2
    resolved.append(result)
    records.append({'block':number,'ours':ours,'theirs':theirs,'resolved':result,'reason':reason})
iterator=iter(resolved)
new=pattern.sub(lambda _:next(iterator),s)
assert not re.search(r'^(<<<<<<< |=======\s*$|>>>>>>> )',new,re.M)
assert not subprocess.check_output(['git','diff','--name-only','--diff-filter=U']).decode().strip().splitlines() != ['gbdraw/web/index.html']
path.write_text(new)
print(json.dumps({'path':str(path),'conflictSourceSha256':hashlib.sha256(s.encode()).hexdigest(),'resolvedSha256':hashlib.sha256(new.encode()).hexdigest(),'resolutions':records},indent=2))
