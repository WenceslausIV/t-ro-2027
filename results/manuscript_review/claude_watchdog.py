"""Run a tools-disabled Claude advisory call and record its liveness and result.

This monitor never launches edits or changes system settings. It exits with the
child, or terminates only its own child after the configured timeout.
"""
import argparse
import datetime as dt
import json
from pathlib import Path
import subprocess
import time


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exe', required=True)
    ap.add_argument('--round', required=True)
    ap.add_argument('--timeout', type=float, default=900)
    ap.add_argument('--resume')
    args = ap.parse_args()
    here = Path(__file__).resolve().parent
    prompt = (here / (args.round + '_prompt.md')).read_text(encoding='utf-8')
    status_path = here / (args.round + '_status.json')
    started = time.monotonic()
    with (here / (args.round + '_response.json')).open('w', encoding='utf-8') as out, \
         (here / (args.round + '_stderr.txt')).open('w', encoding='utf-8') as err:
        command = [args.exe, '--print', '--safe-mode', '--tools', '',
                   '--strict-mcp-config', '--output-format', 'json']
        if args.resume:
            command += ['--resume', args.resume]
        proc = subprocess.Popen(
            command,
            stdin=subprocess.PIPE, stdout=out, stderr=err,
            encoding='utf-8', cwd=here.parents[1],
            creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
        proc.stdin.write(prompt)
        proc.stdin.close()
        while True:
            elapsed = time.monotonic() - started
            rc = proc.poll()
            state = 'running' if rc is None else ('completed' if rc == 0 else 'failed')
            if rc is None and elapsed > args.timeout:
                proc.terminate()
                try:
                    proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()
                rc, state = proc.returncode, 'timed_out'
            status = dict(round=args.round, state=state, child_pid=proc.pid,
                          elapsed_s=round(elapsed, 1), exit_code=rc,
                          updated_utc=dt.datetime.now(dt.timezone.utc).isoformat())
            status_path.write_text(json.dumps(status, indent=2), encoding='utf-8')
            print(json.dumps(status), flush=True)
            if rc is not None:
                break
            time.sleep(20)


if __name__ == '__main__':
    main()
