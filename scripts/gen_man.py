#!/usr/bin/env python3
"""Generate man/*.Rd from the roxygen blocks in R/*.R.

`roxygen2` is the right tool and is not installable here -- `pkgload`, one of its
dependencies, fails to build in this environment. So this script does the narrow subset
that `R CMD check` actually requires:

  * one `.Rd` per exported object, named after it;
  * `\\name`, `\\alias`, `\\title`, `\\description`, `\\usage`, `\\keyword`.

`\\usage` is read from the function's own formals rather than written by hand, because
`R CMD check` compares the deparsed call in `\\usage` against the actual definition and
reports a mismatch as an ERROR. Deriving it means a signature change cannot silently
desynchronise the docs.

The prose comes from the `#'` block immediately above the definition, minus roxygen tags,
so the `.Rd` and the source comment cannot disagree. Anything the block does not cover
(arguments, details) falls back to pointing at the upstream CellChat documentation, which
is where the authoritative argument semantics live -- this package's contract is to be a
drop-in for those functions, not to restate them.

Run:  python3 scripts/gen_man.py
"""
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MAN = ROOT / "man"
R_DIR = ROOT / "R"


def exported_names() -> list[str]:
    ns = (ROOT / "NAMESPACE").read_text()
    return re.findall(r"^export\((.+)\)$", ns, re.M)


def r_script(body: str) -> str:
    """Run a snippet through R and return its stdout, so R itself does the parsing."""
    # `R_LIBS` has to point at the repo-local library. `scripts/reinstall.sh` installs into
    # `<root>/.rlib` precisely so nothing touches the user's site library, and a subprocess that
    # does not inherit that setting cannot `library(CellChat)` -- so every `\\usage` lookup fails
    # and the script reports "no roxygen block found" for exports whose docs are perfectly fine.
    # The failure is silent, which is the part worth fixing: R writes the load error to stderr.
    env = dict(os.environ)
    lib = ROOT / ".rlib"
    if lib.is_dir():
        existing = env.get("R_LIBS", "")
        env["R_LIBS"] = f"{lib}:{existing}" if existing else str(lib)
    out = subprocess.run(
        ["R", "--vanilla", "-q", "-e", body],
        capture_output=True, text=True, cwd=ROOT, env=env,
    )
    return out.stdout


def blocks_for(fname: str) -> dict[str, str]:
    """Map function name -> roxygen prose (tags stripped) for one R file."""
    lines = (R_DIR / fname).read_text().split("\n")
    result: dict[str, str] = {}
    pending: list[str] = []
    for i, line in enumerate(lines):
        stripped = line.lstrip()
        if stripped.startswith("#'"):
            pending.append(stripped[2:].lstrip())
            continue
        m = re.match(r"^([A-Za-z._][A-Za-z0-9._]*)\s*<-\s*function", line)
        if m and pending:
            prose = []
            for p in pending:
                t = p.strip()
                if t.startswith("@"):
                    continue
                prose.append(t)
            # Collapse runs of blank lines and drop leading/trailing blanks.
            while prose and not prose[0]:
                prose.pop(0)
            while prose and not prose[-1]:
                prose.pop()
            result[m.group(1)] = "\n".join(prose)
        if not stripped.startswith("#'"):
            pending = []
    return result


def usage_of(name: str) -> str | None:
    """`\\usage` from the function's formals, via R's own deparse."""
    src = (R_DIR / "modeling.R").read_text()
    esc = re.sub(r"([.\\+*?\[\]$^])", r"\\\1", name)
    body = (
        'suppressWarnings(suppressMessages(library(CellChat)));'
        f'f <- get("{name}", envir = asNamespace("CellChat"));'
        'cat(paste(deparse(args(f)), collapse = " "), sep = "\\n")'
    )
    out = r_script(body)
    # `args()` gives `function (...) NULL` for the `...` forwarders; the deparsed formals
    # are what we want, and R prints them one per line.
    lines = [l for l in out.split("\n") if l.strip() and not l.startswith(">")]
    if not lines:
        return None
    formals = " ".join(l.strip() for l in lines)
    formals = re.sub(r"\s+", " ", formals).strip()
    if not formals.startswith("function"):
        return None
    inner = formals[len("function"):].strip()
    # `args()` renders `function (a, b = 1) NULL` -- drop the trailing NULL and one layer of
    # parentheses, or `\usage` reads `f((a, b = 1))` and check reports a codoc mismatch.
    inner = re.sub(r"\s*NULL$", "", inner).strip()
    if inner.startswith("(") and inner.endswith(")"):
        inner = inner[1:-1].strip()
    if inner in ("", "..."):
        # A zero-argument function must be documented as `f()`, not `f(...)`: `R CMD check`
        # compares the two and reports a codoc mismatch otherwise.
        return (f"{name}()", []) if inner == "" else (f"{name}(...)", ["..."])
    names = []
    depth = 0
    cur = ""
    for ch in inner:
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        if ch == "," and depth == 0:
            names.append(cur.strip())
            cur = ""
        else:
            cur += ch
    if cur.strip():
        names.append(cur.strip())
    return f"{name}({inner})", [n.split("=")[0].strip() for n in names]


def arg_doc(a: str) -> str:
    """One line per argument.

    Deliberately a pointer rather than a paraphrase. These functions exist to be argument-
    compatible with their `CellChat` counterparts, and the authoritative semantics -- including
    the several that are subtler than they look, like `thresh.pc` being compared against a
    fraction in one place and a percentage in another -- live in CellChat's own
    documentation and in `docs/SEMANTICS.md`. Restating them here would create a second copy
    that can drift.
    """
    if a == "...":
        return ("Arguments passed through unchanged. See the corresponding function in the "
                "\\pkg{CellChat} package.")
    return (f"Argument `{a}` of the corresponding \\pkg{{CellChat}} function; see its "
            "documentation and \\code{docs/SEMANTICS.md} in this package.")


def title_of(prose: str, name: str) -> str:
    for line in prose.split("\n"):
        t = line.strip()
        if t:
            return t.rstrip(".")
    return name


def escape(s: str) -> str:
    return s.replace("\\", "\\\\").replace("%", "\\%").replace("{", "\\{").replace("}", "\\}")


def prune(names: list[str]) -> None:
    r"""Delete every stale *function* man/*.Rd, leaving class, data and method pages alone.

    Without this the directory only ever grows. `main()` iterates `exported_names()`, so a
    file for anything the package stopped exporting is simply never rewritten -- it is also
    never removed, and it survives every later run. `R CMD check` then reports it under
    "checking for code/documentation mismatches" as *usage in documentation object X but not
    in code*, for a file no amount of re-running the generator will ever touch.

    That is not hypothetical: 17 such files (`thresholdedMean`, `alpha_centrality`,
    `node_distance`, ...) survived here, documenting upstream internals the package does not
    define.

    This module owns function help pages only. Class, data and S4 method pages are written by
    a different route -- `R/CellChat_class.R` and the `data/` exports -- and `exported_names()`
    cannot see them, since it matches `^export(...)` and those objects are not exported that
    way. Pruning on "not written this run" alone therefore deleted eight legitimate pages
    (`CellChat-class.Rd`, `CellChatDB.*.Rd`, `PPI.*.Rd`, `show-CellChat-method.Rd`,
    `AnyMatrix-class.Rd`) along with the 17. `\docType` is what tells the two apart, and the
    name suffixes are checked as well so a page that omits `\docType` is still spared.

    The five dataset pages are a separate, deliberate absence and are the one case where
    sparing a `\docType{data}` page is wrong: `CellChatDB.{human,mouse,zebrafish}` and
    `PPI.{human,mouse}` are byte-identical duplicates of the copies in
    `inst/upstream/<sha>/data/`, and shipping both is most of what puts the source tarball over
    CRAN's 5 MB limit, so `.Rbuildignore` keeps `data/*.rda` out and the pages went with them.
    Documenting an object the package does not ship is the `code/documentation mismatches`
    WARNING, which is the same failure this function exists to prevent -- so do not re-add
    them. `tests/test-public-surface.R` still asserts those datasets are byte-identical to
    upstream, by loading them from `inst/upstream` directly rather than through `data()`.
    Nothing generates these pages, so nothing has to be taught to suppress them; this note is
    here because they are the obvious thing to restore and restoring them breaks the check.
    """
    keep = set(names)
    for path in sorted(MAN.glob("*.Rd")):
        if path.stem in keep:
            continue
        text = path.read_text(errors="replace")
        if re.search(r"\\docType\{(data|class|methods)\}", text):
            continue
        if path.stem.endswith("-class") or "-method" in path.stem:
            continue
        path.unlink()


def main() -> int:
    blocks: dict[str, str] = {}
    for f in sorted(p.name for p in R_DIR.glob("*.R")):
        blocks.update(blocks_for(f))
    names = exported_names()
    MAN.mkdir(exist_ok=True)
    missing = []
    for name in names:
        upstream_doc = ROOT / "inst/upstream/CellChat-75253cd0/man" / f"{name}.Rd"
        if upstream_doc.exists():
            (MAN / upstream_doc.name).write_text(upstream_doc.read_text())
            continue
        prose = blocks.get(name)
        if prose is None:
            missing.append(name)
            continue
        got = usage_of(name)
        if got is None:
            missing.append(name)
            continue
        usage, argnames = got
        title = title_of(prose, name)
        body = prose.split("\n", 1)[1] if "\n" in prose else ""
        body = body.strip() or prose.strip()
        rd = [
            f"\\name{{{name}}}",
            f"\\alias{{{name}}}",
            f"\\title{{{escape(title)}}}",
            "\\description{",
            escape(body),
            "}",
            "\\usage{",
            escape(usage),
            "}",
            "\\arguments{",
            *[
                f"\\item{{{a}}}{{{arg_doc(a)}}}"
                for a in (argnames or ["..."])
            ],
            "}",
            "\\keyword{internal}",
            "",
        ]
        (MAN / f"{name}.Rd").write_text("\n".join(rd))
    prune(names)
    print(f"wrote {len(names) - len(missing)} of {len(names)} Rd files")
    if missing:
        print("no roxygen block found for: " + ", ".join(missing), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
