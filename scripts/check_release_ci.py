#!/usr/bin/env python3
"""REL-13: refuse publication unless the latest CI run passed for this exact commit."""

import argparse
import json
import os
import re
from urllib.parse import urlencode
from urllib.request import Request, urlopen


def successful_ci_run(runs: list[dict], sha: str, repository: str) -> dict:
    matching = [
        run
        for run in runs
        if run.get("head_sha") == sha
        and run.get("path") == ".github/workflows/ci.yml"
        and run.get("event") == "push"
        and run.get("head_branch") == "main"
        and (run.get("head_repository") or {}).get("full_name") == repository
    ]
    if not matching:
        raise ValueError(f"No main-branch CI run exists for exact commit {sha}.")
    latest = max(matching, key=lambda run: (int(run["run_number"]), int(run["run_attempt"])))
    if latest.get("status") != "completed" or latest.get("conclusion") != "success":
        raise ValueError(
            f"Latest CI run {latest.get('id')} for {sha} is "
            f"{latest.get('status')}/{latest.get('conclusion')}; publication refused."
        )
    return latest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sha", required=True)
    parser.add_argument("--repository", default=os.environ.get("GITHUB_REPOSITORY", ""))
    args = parser.parse_args()
    if not re.fullmatch(r"[0-9a-f]{40}", args.sha):
        parser.error("--sha must be a complete commit SHA")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", args.repository):
        parser.error("--repository must be owner/repository")
    token = os.environ.get("GH_TOKEN")
    if not token:
        parser.error("GH_TOKEN is required to read workflow runs")
    query = urlencode({"head_sha": args.sha, "event": "push", "per_page": 100})
    request = Request(
        f"https://api.github.com/repos/{args.repository}/actions/workflows/ci.yml/runs?{query}",
        headers={
            "Accept": "application/vnd.github+json",
            "Authorization": f"Bearer {token}",
            "X-GitHub-Api-Version": "2026-03-10",
        },
    )
    with urlopen(request, timeout=30) as response:
        runs = json.load(response)["workflow_runs"]
    try:
        run = successful_ci_run(runs, args.sha, args.repository)
    except (ValueError, KeyError, TypeError) as exc:
        parser.error(str(exc))
    print(json.dumps({"source_sha": args.sha, "ci_run_id": run["id"], "ci_url": run["html_url"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
