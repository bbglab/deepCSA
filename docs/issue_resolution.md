# Issue Resolution Documentation

## Introduction

This document provides an overview of the issues identified during the development process and the corresponding resolutions.

## Table of Contents

1. [BedTools merge: 0-len interval wrongly merged](#issue-1)
2. [Issue 2: Description of the issue](#issue-2)
3. [Issue 3: Description of the issue](#issue-3)
4. ...

## 1. BedTools merge: 0-len interval wrongly merged  <a name="issue-1"></a>

The `BedTools merge` command adds a -1 and a +1 to 0-len interval when merged. This does not happen if the line is not merged to anything.


- Example:
```sh
>> cat test.bed 
chr1       119915306       119915306
..    ...    ....
chr1       119916711       119916711

>> bedtools merge -i test.bed
chr1    119915305       119916712  # expected 119915306 119916711
```

This is already known by the devs of bedtools [arq5x/bedtools2#359](https://github.com/arq5x/bedtools2/issues/359#issuecomment-2066711835)

### Processes affected: 
- `CREATEPANELS:*`


### GitHub tracking
- **GitHub Issue:** [bbglab/deepCSA#82](https://github.com/bbglab/deepCSA/issues/82)
- **GitHub PR (Resolution):** [bbglab/deepCSA#88](https://github.com/bbglab/deepCSA/pull/88)

### Resolution <a name="resolution"></a>

A dirty and quick solution was decided: Adding a line of awk that transform the output of bedtool merge into the expected coordinates.

```sh
awk 'BEGIN { FS = "\\t"; OFS = "\\t" } { if (\$2 != \$3) { \$2 += 1 ; \$3 -= 1 } else { \$2 = \$2 ; \$3 = \$3 }; print }' > [filename].bed
```


## 2. Issue 2: Description of the issue <a name="issue-2"></a>

Description of the issue, including its impact and context with example.

### Processes affected: 
- `<PROCESSNAME>`


### GitHub tracking

- **GitHub Issue:** [Link to GitHub Issue](link_to_issue)
- **GitHub PR (Resolution):** [Link to GitHub PR](link_to_pr)

### Resolution <a name="resolution"></a>

## 3. Issue 3: Description of the issue <a name="issue-3"></a>

Description of the issue, including its impact and context.

### Processes affected: 
- `<PROCESSNAME>`


### GitHub tracking

- **GitHub Issue:** [Link to GitHub Issue](link_to_issue)
- **GitHub PR (Resolution):** [Link to GitHub PR](link_to_pr)

### Resolution <a name="resolution"></a>


## Conclusion

In this document, we have documented the issues identified during development and their corresponding resolutions.
