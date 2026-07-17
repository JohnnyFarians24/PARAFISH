"""Running external tools, and deciding whether they actually worked.

The rule this module exists to enforce, learned the expensive way from blastn:

    **Judge an external tool by its OUTPUT, not by its exit code.**

`blastn` exits 0 having printed "(Severe Error) ... CPU usage limit was exceeded" for
every query and written an empty result file. Code that trusted `check=True` therefore
read a server refusal as a successful search that found nothing, and the lie travelled
all the way to `source_status` as `no_hits` — the one thing that column exists to
prevent. Exit codes are a claim the tool makes about itself; the output is evidence.

The second half of the rule, equally important and easy to get wrong in the other
direction:

    **An empty output is not automatically a failure.**

A quality filter that legitimately passes nothing, or a search that genuinely matches
nothing, produce the same empty file as a crash. Emptiness alone cannot tell them apart —
only the tool's own error output can. So the distinction is:

    empty output + error signature   -> FAILURE (it was refused / it died)
    empty output + clean run         -> a real, honest, empty result

Blanket "empty means broken" would turn every honest absence into a false alarm, which is
just the original sin wearing the opposite mask.
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from typing import Optional, Sequence

# Signatures that mean "this tool did not do the work", regardless of exit code. Kept
# generic here; a caller with tool-specific signatures passes its own on top.
GENERIC_FAILURE_SIGNATURES = (
	"severe error",
	"internal_error",
	"segmentation fault",
	"out of memory",
	"cannot allocate",
	"killed",
	"traceback (most recent call last)",
)


@dataclass
class ToolResult:
	"""What an external tool did, judged by evidence rather than by its own claim."""

	command: Sequence[str]
	returncode: int
	stdout: str
	stderr: str
	# The first line of output matching a failure signature, if any. Truthy => the tool
	# reported a real failure even if it exited 0.
	failure_signature: Optional[str] = None
	timed_out: bool = False

	@property
	def output(self) -> str:
		joined = (self.stdout or "") + ("\n" if self.stdout and self.stderr else "") + (self.stderr or "")
		return joined

	@property
	def ok(self) -> bool:
		"""True only if the tool exited cleanly AND said nothing that means it failed."""
		return self.returncode == 0 and not self.timed_out and self.failure_signature is None

	@property
	def killed_by_signal(self) -> bool:
		"""Terminated by the OS — on this box, almost always the OOM killer."""
		return self.returncode is not None and self.returncode < 0


def find_failure_signature(text: str, extra: Sequence[str] = ()) -> Optional[str]:
	"""The first line of ``text`` that means the tool failed, or ``None``.

	Matching is per line and case-insensitive so the caller gets back something worth
	putting in a log, not just a boolean.
	"""
	signatures = tuple(s.lower() for s in (*GENERIC_FAILURE_SIGNATURES, *extra))
	for line in (text or "").splitlines():
		lowered = line.lower()
		if any(signature in lowered for signature in signatures):
			return line.strip()
	return None


def run_tool(
	command: Sequence[str],
	*,
	timeout: Optional[int] = None,
	env: Optional[dict] = None,
	cwd: Optional[str] = None,
	failure_signatures: Sequence[str] = (),
	input_text: Optional[str] = None,
) -> ToolResult:
	"""Run ``command`` and report what it actually did.

	Never raises on a failing tool: the caller gets a ``ToolResult`` and decides. That is
	deliberate — ``check=True`` throws away the tool's own account of what went wrong, and
	that account is the only thing that distinguishes a refusal from an honest empty.
	"""
	try:
		completed = subprocess.run(
			list(command),
			capture_output=True,
			text=True,
			check=False,
			timeout=timeout,
			env=env,
			cwd=cwd,
			input=input_text,
		)
	except subprocess.TimeoutExpired as exc:
		return ToolResult(
			command=command,
			returncode=1,
			stdout=(exc.stdout or "") if isinstance(exc.stdout, str) else "",
			stderr=(exc.stderr or "") if isinstance(exc.stderr, str) else "",
			failure_signature=f"timed out after {timeout}s",
			timed_out=True,
		)
	result = ToolResult(
		command=command,
		returncode=completed.returncode,
		stdout=completed.stdout or "",
		stderr=completed.stderr or "",
	)
	result.failure_signature = find_failure_signature(result.output, failure_signatures)
	return result
