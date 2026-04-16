# PSF Privacy Hooks

Git hooks that enforce the public/private split between the Codeberg public edition and private editions of this repository.

## What this protects

Public branch (`main`) is the Codeberg-facing edition. Commits here must not contain:

1. **Gitignored paths** — i.e. paths listed in the private section of `.gitignore`, which may not be committed to the public branch even via `git add -f`.
2. **Blocked keywords** — a small set of public program acronyms and partner-specific phrases whose *presence* in a public commit would reveal the private fork's focus. Names of individuals and private companies are maintained in a *local-only* keyword list that is never committed.

Private branches (`private`, `private/*`, `bitbucket`, `bitbucket/*`) are exempt.

## What this does NOT protect

- **Secrets** — API tokens, passwords, private keys. Use a dedicated secret scanner (e.g. `git-secrets`, `trufflehog`) for that.
- **License violations** — unrelated concern.
- **Humans who intentionally bypass** — the emergency escape hatch (`PSF_SKIP_PRIVACY_HOOK=1`) is logged but not blocked. Review the log periodically.

## Install

From the repo root:

```bash
scripts/git-hooks/install.sh
```

This sets `git config core.hooksPath scripts/git-hooks`. All hook scripts in this directory become active.

Run once per clone. The setting is stored in `.git/config` (per clone, not committed), so a fresh clone requires running it again.

## Uninstall

```bash
scripts/git-hooks/uninstall.sh
```

Clears `core.hooksPath` but leaves the hook scripts in place.

## Files

| File | Purpose | Committed |
|---|---|---|
| `pre-commit` | Runs on every `git commit` on non-private branches | yes |
| `pre-push` | Runs on every `git push`; scans entire push range | yes |
| `install.sh` | Sets `core.hooksPath`, seeds local keyword file | yes |
| `uninstall.sh` | Unsets `core.hooksPath` | yes |
| `blocked-keywords.txt` | Public-safe keyword patterns (acronyms, program names) | yes |
| `blocked-keywords.local.txt.example` | Template for local additions | yes |
| `blocked-keywords.local.txt` | **User's real sensitive patterns** (names, partner companies) | **no, gitignored** |
| `.skipped.log` | Audit log of emergency bypasses | no, gitignored |

## Behavior matrix

| Branch | Remote | pre-commit | pre-push |
|---|---|---|---|
| `private` | any | **skipped** | skipped if remote is private |
| `main` | any | enforces | enforces if remote is public |
| detached HEAD | any | enforces | enforces |
| anything else | codeberg/github/gitlab | enforces | enforces |
| anything else | bitbucket | enforces | skipped |
| anything else | unclassified | enforces | enforces (warn) |

Private remotes are identified by host in `pre-push`: `bitbucket.org`, `topmistbpu`. Add more by editing `PRIVATE_REMOTE_PATTERNS` in the `pre-push` script.

## How rejection looks

```
✗ privacy hook rejected this commit on branch 'main'.

✗ Gitignored paths were force-added (use of 'git add -f'?):
    docs/roadmap-v07-bpu-private.md

How to resolve
  - Unstage the offending files: git reset HEAD <file>
  - If the content belongs in the private edition, switch branches:
      git checkout private && git cherry-pick <commit>
  - If this is a genuine false positive, edit the keyword list:
      scripts/git-hooks/blocked-keywords.local.txt  (gitignored)
  - For an emergency bypass (logged):
      PSF_SKIP_PRIVACY_HOOK=1 git commit ...
```

## Escape hatches

Three ways to bypass:

1. **`PSF_SKIP_PRIVACY_HOOK=1`** — documented bypass, append-only log at `.skipped.log`. Use sparingly.
2. **`git commit --no-verify`** — skips all hooks. Git's escape hatch, not ours. Leaves no log.
3. **Edit `blocked-keywords.local.txt`** — correct fix for a genuine false positive. Commit history of the local file is still local-only.

## Extending

- **New private path:** add to `.gitignore`. The path-check is automatic.
- **New public-safe keyword:** edit `blocked-keywords.txt` and commit.
- **New sensitive name or partner:** edit `blocked-keywords.local.txt` (do not commit).
- **New private branch name pattern:** edit `PRIVATE_BRANCH_PATTERNS` in `pre-commit`.
- **New private remote host:** edit `PRIVATE_REMOTE_PATTERNS` in `pre-push`.

## Testing the hook locally

Dry-run test (stages a file, runs hook, does not commit):

```bash
# On the public branch:
git checkout main
echo "This mentions O-Circuit and BPU testing" > /tmp/test-violation.md
git add -f /tmp/test-violation.md
git commit -m "test" --dry-run   # hook runs, rejects
git reset HEAD /tmp/test-violation.md
```

Expected: rejection message listing the blocked keyword matches.

## Performance

- `pre-commit`: scans staged content only. Fast (< 1 s for typical commits; large binary files are skipped by extension).
- `pre-push`: scans every commit in the push range. Can be slow for large pushes (new branch with many commits). Files > 512 KB are skipped by size.

To tune: edit `MAX_SCAN_BYTES` and `SKIP_CONTENT_SCAN_EXTS` in the hook scripts.

## Philosophy

This hook is a **safety net, not a security boundary**. A motivated insider can bypass it (two ways shown above). What it prevents is the routine case: a hurried commit that accidentally includes a private file, or a copy-paste that pulls in consortium-partner names. For that, it is sufficient.

The real security is:

- Separate branches (`main` public, `private` private)
- `.gitignore` listing every private path explicitly
- Separate remotes (`codeberg` public, `bitbucket` private) with independent access control
- Periodic review of `.skipped.log`

The hook automates the easy failure modes so attention can focus on the hard ones.
