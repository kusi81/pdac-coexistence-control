"""
중앙디렉토리(EOCD)가 잘린 store(무압축) zip 복구 — 데이터 디스크립터(streaming) 대응.

figshare가 zhou_processed.zip 끝을 잘라 서빙 → EOCD 없음. 게다가 store + 데이터
디스크립터 방식이라 로컬 헤더에 크기가 0. 따라서 각 멤버의 끝은 '다음 로컬 파일 헤더
(PK0304)' 앞의 데이터 디스크립터(PK0708 + crc + comp + uncomp)로 판단한다.
디스크립터의 uncomp 크기가 (경계-본문시작)과 일치하는지로 오탐(HDF5 내부 우연 시그니처)을 배제.
잘린 마지막 멤버는 자동으로 제외된다(앞쪽 온전한 h5ad만 저장).
"""
import sys, os, struct, zlib
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

LFH = b"\x50\x4b\x03\x04"     # local file header
DD = b"\x50\x4b\x07\x08"      # data descriptor
HDF5 = b"\x89HDF\r\n\x1a\n"   # HDF5 파일 시그니처(h5ad 검증용)


def _find_member_end(data, body_start, n):
    """데이터 디스크립터 경계를 찾아 data_end 반환. 멤버는 DEFLATE 압축이므로
    디스크립터의 **압축크기(comp)** 가 (data_end - body_start)와 일치하는지로 검증.
    디스크립터: [PK0708?] crc(4) comp uncomp — 크기 4B 또는 8B(ZIP64)."""
    # (desc_len, has_sig, size_bytes, comp_off_from_end)
    variants = [(16, True, 4, 8), (12, False, 4, 8), (24, True, 8, 16), (20, False, 8, 16)]
    search = body_start
    while True:
        nxt = data.find(LFH, search + 1)
        cand_end = nxt if nxt != -1 else n
        for desc, has_sig, sz, comp_off in variants:
            ds = cand_end - desc
            if ds < body_start:
                continue
            if has_sig and data[ds:ds + 4] != DD:
                continue
            fmt = "<Q" if sz == 8 else "<I"
            comp = struct.unpack(fmt, data[cand_end - comp_off:cand_end - comp_off + sz])[0]
            if comp == (cand_end - desc) - body_start:
                return cand_end - desc, (nxt if nxt != -1 else n)
        if nxt == -1:
            return None, n
        search = nxt


def salvage(zip_path, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    data = open(zip_path, "rb").read()
    n = len(data)
    pos, saved = 0, []
    while pos + 30 <= n and data[pos:pos + 4] == LFH:
        (ver, flag, method, mtime, mdate, crc, comp, uncomp,
         name_len, extra_len) = struct.unpack("<HHHHHIIIHH", data[pos + 4:pos + 30])
        name = data[pos + 30:pos + 30 + name_len].decode("utf-8", "replace")
        body_start = pos + 30 + name_len + extra_len
        if name.endswith("/"):                    # 디렉토리 엔트리 (본문 없음)
            pos = body_start
            continue
        if (flag & 0x08) or comp == 0:            # 데이터 디스크립터(크기 미상)
            data_end, next_pos = _find_member_end(data, body_start, n)
            if data_end is None:
                print(f"  [잘림/마지막] {os.path.basename(name)} - 경계 못찾음, 종료")
                break
        else:
            data_end, next_pos = body_start + comp, body_start + comp
            if data_end > n:
                print(f"  [잘림] {os.path.basename(name)} - 종료")
                break
        if name.endswith(".h5ad"):
            raw = data[body_start:data_end]
            try:
                content = zlib.decompress(raw, -15) if method == 8 else raw
            except Exception as e:  # noqa: BLE001
                print(f"  [실패] {os.path.basename(name)}: 압축해제 오류 {e} (잘림?)")
                pos = next_pos
                continue
            ok = content[:8] == HDF5
            out_path = os.path.join(out_dir, os.path.basename(name))
            with open(out_path, "wb") as f:
                f.write(content)
            saved.append(os.path.basename(name))
            print(f"  [복구{'' if ok else '?'}] {os.path.basename(name)}  "
                  f"({len(content) / 1e6:.0f} MB{'' if ok else ', HDF5헤더 불일치'})")
        pos = next_pos
    return saved


if __name__ == "__main__":
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    zp = os.path.join(here, "data", "zhou_processed.zip")
    od = os.path.join(here, "data", "zhou")
    s = salvage(zp, od)
    print(f"총 {len(s)}개 h5ad 복구 -> {od}")
