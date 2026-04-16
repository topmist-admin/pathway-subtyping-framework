#!/usr/bin/env bash
# Remove PSF privacy hooks by clearing core.hooksPath.
# Does NOT delete the hook scripts themselves.

set -euo pipefail

current="$(git config --get core.hooksPath || true)"
target="scripts/git-hooks"

if [[ -z "${current}" ]]; then
    echo "core.hooksPath is not set — hooks are not installed"
    exit 0
fi

if [[ "${current}" != "${target}" ]]; then
    echo "core.hooksPath is '${current}' (not PSF's '${target}') — leaving alone"
    exit 0
fi

git config --unset core.hooksPath
echo "unset core.hooksPath (was ${current})"
echo "hook scripts remain at ${target}/ — delete them if you want."
