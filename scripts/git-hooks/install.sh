#!/usr/bin/env bash
# Install PSF privacy hooks by pointing git at this directory.
#
# Uses `core.hooksPath` (single git config change) rather than symlinking
# into .git/hooks/. This means the hooks are source-controlled and travel
# with the repo; no per-clone setup beyond running this script once.
#
# Idempotent. Safe to re-run.

set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel)"
HOOKS_DIR="${REPO_ROOT}/scripts/git-hooks"

if [[ ! -d "${HOOKS_DIR}" ]]; then
    echo "error: ${HOOKS_DIR} does not exist" >&2
    exit 1
fi

# Make all hook scripts executable.
for hook in pre-commit pre-push; do
    if [[ -f "${HOOKS_DIR}/${hook}" ]]; then
        chmod +x "${HOOKS_DIR}/${hook}"
    fi
done

chmod +x "${HOOKS_DIR}/install.sh"
[[ -f "${HOOKS_DIR}/uninstall.sh" ]] && chmod +x "${HOOKS_DIR}/uninstall.sh"

# Point git at our hooks directory.
current="$(git config --get core.hooksPath || true)"
target="scripts/git-hooks"

if [[ "${current}" == "${target}" ]]; then
    echo "core.hooksPath already set to ${target} — nothing to do"
else
    if [[ -n "${current}" ]]; then
        echo "note: previous core.hooksPath was '${current}' — overwriting"
    fi
    git config core.hooksPath "${target}"
    echo "set core.hooksPath = ${target}"
fi

# Seed a local keyword file if the user does not have one yet.
local_kw="${HOOKS_DIR}/blocked-keywords.local.txt"
example_kw="${HOOKS_DIR}/blocked-keywords.local.txt.example"

if [[ ! -f "${local_kw}" && -f "${example_kw}" ]]; then
    cp "${example_kw}" "${local_kw}"
    echo "created ${local_kw} from example (edit to add sensitive names)"
fi

# Quick sanity check: warn if the user is on the public branch and already
# has a working-tree mismatch that would fail on commit.
branch="$(git symbolic-ref --short HEAD 2>/dev/null || echo '')"
if [[ "${branch}" == "main" ]]; then
    echo "reminder: you are on '${branch}' — the PUBLIC branch. Hook will enforce."
fi

echo "done."
echo
echo "Test:  echo '# test' > /tmp/psf_test && git add /tmp/psf_test &&  \\"
echo "       git commit -m test --dry-run     # hook runs, nothing committed"
echo
echo "Uninstall:  git config --unset core.hooksPath"
