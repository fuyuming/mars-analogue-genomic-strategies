#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""可复现的 MAG FASTA 下载器（纯标准库，CLI 参数，无私有绝对路径硬编码）。

实现 reliable-large-download 技能的下载要求：
  - 全局锁 + 每文件锁（mkdir 原子目录锁）
  - .part 断点续传（HTTP Range，不因网络错误删除进度）
  - 连续 3 次无进度即停止该文件并诊断（不无限重试）
  - 已校验文件不可变、重启即跳过，绝不覆盖
  - 完成前校验：精确字节数 + 官方 MD5 + gzip 完整性 + FASTA 头（'>'）
  - 记录 SHA256；超长/校验失败文件可恢复地隔离到 .quarantine/
  - 输出 receipt（每文件字节/MD5/SHA256/状态）与 status JSON

用法：
  python3 download_mag_fastas.py --manifest download_manifest.tsv \
      --outdir /path/genomes [--parallel 3] [--dry-run]
"""
import argparse
import gzip
import hashlib
import json
import os
import socket
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

UA = "nee56-mag-downloader/1.0 (stdlib)"
CHUNK = 1 << 20


def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def size_of(p):
    try:
        return os.path.getsize(p)
    except OSError:
        return 0


class Logger:
    def __init__(self, path):
        self.path = path

    def __call__(self, *a):
        line = "[%s] %s" % (now(), " ".join(str(x) for x in a))
        print(line, flush=True)
        with open(self.path, "a", encoding="utf-8") as f:
            f.write(line + "\n")


def verify(path, md5_exp, size_exp):
    """返回 (code, tag)。0=通过。"""
    if not os.path.isfile(path):
        return 1, "missing"
    if size_of(path) != int(size_exp):
        return 2, "size_mismatch"
    h = hashlib.md5()
    with open(path, "rb") as f:
        for c in iter(lambda: f.read(CHUNK), b""):
            h.update(c)
    if h.hexdigest() != md5_exp:
        return 3, "md5_mismatch"
    try:
        with gzip.open(path, "rb") as g:
            first = g.read(4096)
            if not first.startswith(b">"):
                return 5, "fasta_header_bad"
            while g.read(CHUNK):
                pass
    except Exception:  # noqa: BLE001
        return 4, "gzip_bad"
    return 0, "ok"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for c in iter(lambda: f.read(CHUNK), b""):
            h.update(c)
    return h.hexdigest()


class Downloader:
    def __init__(self, outdir, parallel, max_rounds, log):
        self.outdir = outdir
        self.parallel = parallel
        self.max_rounds = max_rounds
        self.log = log
        self.locks = os.path.join(outdir, ".locks")
        self.quar = os.path.join(outdir, ".quarantine")
        self.glock = os.path.join(self.locks, "orchestrator.lockdir")
        os.makedirs(self.locks, exist_ok=True)
        os.makedirs(self.quar, exist_ok=True)

    def acquire_global(self):
        if os.path.isdir(self.glock):
            owner = 0
            pidf = os.path.join(self.glock, "pid")
            if os.path.isfile(pidf):
                try:
                    owner = int(open(pidf).read().strip())
                except Exception:  # noqa: BLE001
                    owner = 0
            alive = False
            if owner > 0:
                try:
                    os.kill(owner, 0)
                    alive = True
                except OSError:
                    alive = False
            if alive:
                return False
            os.rename(self.glock, "%s.stale.%d" % (self.glock, int(time.time())))
        try:
            os.mkdir(self.glock)
        except OSError:
            return False
        with open(os.path.join(self.glock, "pid"), "w") as f:
            f.write(str(os.getpid()))
        return True

    def release_global(self):
        try:
            os.remove(os.path.join(self.glock, "pid"))
        except OSError:
            pass
        try:
            os.rmdir(self.glock)
        except OSError:
            pass

    def quarantine(self, path, reason):
        if not os.path.exists(path):
            return
        tgt = os.path.join(self.quar, "%s.%s.%d" % (os.path.basename(path), reason, int(time.time())))
        try:
            os.rename(path, tgt)
            self.log("QUARANTINE", os.path.basename(path), "reason=" + reason, "->", tgt)
        except OSError:
            pass

    def download_one(self, url, md5, size, fname):
        size = int(size)
        out = os.path.join(self.outdir, fname)
        part = out + ".part"
        lock = os.path.join(self.locks, fname + ".lockdir")
        try:
            os.mkdir(lock)
        except OSError:
            self.log("LOCKED", fname)
            return dict(filename=fname, url=url, status="locked")
        t0 = time.time()
        try:
            code, tag = verify(out, md5, size)
            if code == 0:
                self.log("SKIP", fname, "verified")
                return dict(filename=fname, url=url, expected_bytes=size, expected_md5=md5,
                            actual_bytes=size, sha256=sha256(out), status="verified_skip",
                            seconds=round(time.time() - t0, 1))
            if os.path.isfile(out):
                self.quarantine(out, "invalid_final_" + tag)
            if os.path.isfile(part) and size_of(part) > size:
                self.quarantine(part, "oversize")
            if not os.path.exists(part):
                open(part, "a").close()
            noprog = 0
            for rnd in range(1, self.max_rounds + 1):
                before = size_of(part)
                if before == size:
                    code, tag = verify(part, md5, size)
                    if code == 0:
                        os.rename(part, out)
                        self.log("DONE", fname, "bytes=%d" % size, "round=%d" % rnd)
                        return dict(filename=fname, url=url, expected_bytes=size,
                                    expected_md5=md5, actual_bytes=size, sha256=sha256(out),
                                    status="verified", seconds=round(time.time() - t0, 1))
                    self.quarantine(part, "checksum_fail_" + tag)
                    open(part, "a").close()
                    before = 0
                self.log("GET", fname, "round=%d" % rnd, "resume=%d" % before, "expect=%d" % size)
                try:
                    hdrs = {"User-Agent": UA}
                    if before > 0:
                        hdrs["Range"] = "bytes=%d-" % before
                    req = urllib.request.Request(url, headers=hdrs)
                    with urllib.request.urlopen(req, timeout=120) as r:
                        mode = "ab"
                        if before > 0 and r.status != 206:
                            mode = "wb"  # 服务器忽略 Range，重头来
                        with open(part, mode) as f:
                            while True:
                                c = r.read(CHUNK)
                                if not c:
                                    break
                                f.write(c)
                except (urllib.error.URLError, socket.timeout, OSError) as e:
                    self.log("ERR", fname, repr(e))
                after = size_of(part)
                self.log("CURL", fname, "bytes=%d" % after, "gained=%d" % (after - before))
                if after > size:
                    self.quarantine(part, "oversize")
                    open(part, "a").close()
                    noprog = 0
                elif after == size:
                    code, tag = verify(part, md5, size)
                    if code == 0:
                        os.rename(part, out)
                        self.log("DONE", fname, "bytes=%d" % size, "round=%d" % rnd)
                        return dict(filename=fname, url=url, expected_bytes=size,
                                    expected_md5=md5, actual_bytes=size, sha256=sha256(out),
                                    status="verified", seconds=round(time.time() - t0, 1))
                    self.quarantine(part, "checksum_fail_" + tag)
                    open(part, "a").close()
                    noprog = 0
                elif after == before:
                    noprog += 1
                    self.log("STALL", fname, "no_progress=%d" % noprog)
                    if noprog >= 3:
                        self.log("FAIL", fname, "three_no_progress_stop")
                        return dict(filename=fname, url=url, expected_bytes=size,
                                    expected_md5=md5, actual_bytes=after, status="failed_no_progress",
                                    seconds=round(time.time() - t0, 1))
                    time.sleep(30)
                else:
                    noprog = 0
                    time.sleep(5)
            self.log("FAIL", fname, "exhausted_rounds")
            return dict(filename=fname, url=url, expected_bytes=size, expected_md5=md5,
                        actual_bytes=size_of(part), status="failed_rounds",
                        seconds=round(time.time() - t0, 1))
        finally:
            try:
                os.rmdir(lock)
            except OSError:
                pass


def read_manifest(path):
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        if header[:3] != ["url", "md5", "bytes"]:
            raise SystemExit("manifest 头必须是 url<TAB>md5<TAB>bytes[<TAB>filename]")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3 or not p[0]:
                continue
            fname = p[3] if len(p) > 3 and p[3] else p[0].rsplit("/", 1)[-1]
            rows.append((p[0], p[1], p[2], fname))
    return rows


def main():
    ap = argparse.ArgumentParser(description="nee56 MAG FASTA 可复现下载器（纯 stdlib）")
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--parallel", type=int, default=3)
    ap.add_argument("--max-rounds", type=int, default=50)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--status-json", default=None)
    ap.add_argument("--receipt", default=None)
    a = ap.parse_args()

    rows = read_manifest(a.manifest)
    total = sum(int(r[2]) for r in rows)
    print("manifest=%s files=%d bytes=%d outdir=%s parallel=%d"
          % (a.manifest, len(rows), total, a.outdir, a.parallel), flush=True)
    if a.dry_run:
        print("[dry-run] 不下载。")
        return

    os.makedirs(a.outdir, exist_ok=True)
    status_json = a.status_json or os.path.join(a.outdir, "download_status.json")
    receipt = a.receipt or os.path.join(a.outdir, "download_receipt.tsv")
    log = Logger(os.path.join(a.outdir, "download.log"))
    dl = Downloader(a.outdir, a.parallel, a.max_rounds, log)

    if not dl.acquire_global():
        log("FATAL", "另一个下载器持有全局锁，退出")
        sys.exit(6)
    log("START", "files=%d" % len(rows), "bytes=%d" % total, "parallel=%d" % a.parallel)
    results = []
    try:
        with ThreadPoolExecutor(max_workers=a.parallel) as pool:
            futs = {pool.submit(dl.download_one, u, m, s, fn): fn for (u, m, s, fn) in rows}
            for f in as_completed(futs):
                try:
                    results.append(f.result())
                except Exception as e:  # noqa: BLE001
                    results.append(dict(filename=futs[f], status="error", error=repr(e)))
                if len(results) % 10 == 0:
                    with open(status_json, "w", encoding="utf-8") as fh:
                        json.dump(_summ(rows, results), fh, indent=2, ensure_ascii=False)
    finally:
        dl.release_global()

    # 最终复核
    ok = bad = 0
    for (u, m, s, fn) in rows:
        code, _ = verify(os.path.join(a.outdir, fn), m, int(s))
        if code == 0:
            ok += 1
        else:
            bad += 1
    summary = _summ(rows, results)
    summary["final_verified"] = ok
    summary["final_incomplete_or_bad"] = bad
    summary["finished"] = now()
    with open(status_json, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)
    with open(receipt, "w", encoding="utf-8") as fh:
        fh.write("filename\turl\texpected_bytes\texpected_md5\tactual_bytes\tsha256\tstatus\tseconds\n")
        for r in sorted(results, key=lambda x: x.get("filename", "")):
            fh.write("\t".join(str(r.get(k, "")) for k in
                               ("filename", "url", "expected_bytes", "expected_md5",
                                "actual_bytes", "sha256", "status", "seconds")) + "\n")
    log("FINISH", "verified=%d" % ok, "incomplete_or_bad=%d" % bad, "total=%d" % len(rows))
    sys.exit(0 if bad == 0 else 1)


def _summ(rows, results):
    by = {r.get("filename"): r for r in results}
    verified = sum(1 for r in results if r.get("status", "").startswith("verified"))
    return dict(manifest_files=len(rows), planned_bytes=sum(int(r[2]) for r in rows),
                attempted=len(results), verified=verified,
                failed=sum(1 for r in results if "fail" in r.get("status", "") or r.get("status") == "error"),
                updated=now(), rows=sorted(results, key=lambda x: x.get("filename", "")))


if __name__ == "__main__":
    main()
