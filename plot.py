import re
import numpy as np
import matplotlib.pyplot as plt

FRAME_OPEN_RE = re.compile(r"^<frame\b([^>]*)>$")
FRAME_CLOSE_RE = re.compile(r"^</frame>$")

def parse_frame_attr(attr_text: str) -> dict:
    """
    <frame ...> の ... 部分から、数字っぽいものを拾う用。
    例:
      <frame 0>          -> {"id": 0}
      <frame id=12>      -> {"id": 12}
      <frame step="33">  -> {"id": 33}
    """
    attr_text = attr_text.strip()
    if not attr_text:
        return {"id": None}

    # まず "id=12" や 'id="12"' を探す
    m = re.search(r'\bid\s*=\s*"?(\d+)"?', attr_text)
    if m:
        return {"id": int(m.group(1))}

    # 次に "step=33" とかでもいいことにする
    m = re.search(r'\bstep\s*=\s*"?(\d+)"?', attr_text)
    if m:
        return {"id": int(m.group(1))}

    # 最後に最初に出てくる整数を拾う (<frame 0> とか)
    m = re.search(r'(\d+)', attr_text)
    if m:
        return {"id": int(m.group(1))}

    return {"id": None}


def load_log_frames(path):
    """
    return: frames = [
      {"id": int|None, "waters": (N,3), "soaps": (M,6)},
      ...
    ]
    """
    frames = []
    mode = None            # None / "waters" / "soaps"
    in_frame = False
    cur = None

    def flush_frame():
        nonlocal cur
        if cur is None:
            return
        waters = np.asarray(cur["waters"], dtype=float)
        soaps  = np.asarray(cur["soaps"], dtype=float)

        # 空フレームも許すが、形だけ整える
        if waters.size == 0:
            waters = waters.reshape(0, 3)
        if soaps.size == 0:
            soaps = soaps.reshape(0, 6)

        cur["waters"] = waters
        cur["soaps"] = soaps
        frames.append(cur)
        cur = None

    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue

            # frame open
            m = FRAME_OPEN_RE.match(line)
            if m:
                # もしフレームがネスト/連続してたら前を閉じる
                if in_frame:
                    flush_frame()
                in_frame = True
                mode = None
                attrs = parse_frame_attr(m.group(1))
                cur = {"id": attrs.get("id"), "waters": [], "soaps": []}
                continue

            # frame close
            if FRAME_CLOSE_RE.match(line):
                if in_frame:
                    flush_frame()
                in_frame = False
                mode = None
                continue

            # waters/soaps tags（フレーム外でも一応拾えるようにするなら in_frame を外す）
            if line == "<waters>":
                mode = "waters"
                continue
            if line == "</waters>":
                mode = None
                continue
            if line == "<soaps>":
                mode = "soaps"
                continue
            if line == "</soaps>":
                mode = None
                continue

            # data lines
            if cur is None:
                # frame タグが無いログにも耐性を持たせる（単一フレーム扱い）
                cur = {"id": None, "waters": [], "soaps": []}
                in_frame = True

            vals = line.split()
            if mode == "waters":
                if len(vals) != 3:
                    raise ValueError(f"waters は3列のはず: {line}")
                #for v in vals:
                #    if abs(float(v)) > 300:
                #        print("error")
                cur["waters"].append([float(v) for v in vals])
            elif mode == "soaps":
                if len(vals) != 6:
                    raise ValueError(f"soaps は6列(head xyz tail xyz)のはず: {line}")
                #for v in vals:
                #    if abs(float(v)) > 300:
                #        print("error")
                cur["soaps"].append([float(v) for v in vals])
            else:
                # タグ外の行は無視（必要ならここでエラーにしてもOK）
                pass

    # ファイル末尾で閉じタグ無しでも回収
    if cur is not None:
        flush_frame()

    return frames


def set_equal(ax):
    lims = np.array([ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()], dtype=float)
    c = lims.mean(axis=1)
    r = 0.5 * (lims[:, 1] - lims[:, 0]).max()
    ax.set_xlim(c[0] - r, c[0] + r)
    ax.set_ylim(c[1] - r, c[1] + r)
    ax.set_zlim(c[2] - r, c[2] + r)


# --- 実行部 ---
frames = load_log_frames("./debug_500x500_for_stb.log")
print(f"loaded frames: {len(frames)}")
if len(frames) == 0:
    raise SystemExit("frame が0件。タグ形式を確認して。")

# 全フレーム共通の軸範囲を作る（アニメで視点が揺れない）
all_pts = []
for fr in frames:
    w = fr["waters"]
    s = fr["soaps"]
    if w.size:
        all_pts.append(w)
    if s.size:
        all_pts.append(s[:, 0:3])  # head
        all_pts.append(s[:, 3:6])  # tail
all_pts = np.vstack(all_pts) if all_pts else np.zeros((0,3))

mins = all_pts.min(axis=0) if all_pts.size else np.array([-1,-1,-1], float)
maxs = all_pts.max(axis=0) if all_pts.size else np.array([ 1, 1, 1], float)
center = 0.5 * (mins + maxs)
radius = 0.5 * np.max(maxs - mins) if all_pts.size else 1.0
radius = max(radius, 1e-6)
print(center, radius)
#radius = 60.0
#center = np.array([0.0, 0.0, 0.0])

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection="3d")

def draw_frame(idx: int):
    fr = frames[idx]
    ax.cla()

    w = fr["waters"]
    s = fr["soaps"]

    # waters
    if w.size:
        ax.scatter(w[:,0], w[:,1], w[:,2], s=4, alpha=0.35)

    # soaps
    if s.size:
        H = s[:, 0:3]
        T = s[:, 3:6]
        ax.scatter(H[:,0], H[:,1], H[:,2], s=16)
        ax.scatter(T[:,0], T[:,1], T[:,2], s=16)

        # head-tail bond
        for h, t in zip(H, T):
            ax.plot([h[0], t[0]], [h[1], t[1]], [h[2], t[2]], linewidth=1.0)

    fid = fr["id"]
    ax.set_title(f"frame {idx}" + (f" (id={fid})" if fid is not None else ""))

    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")

    # 固定スケール
    ax.set_xlim(center[0]-radius, center[0]+radius)
    ax.set_ylim(center[1]-radius, center[1]+radius)
    ax.set_zlim(center[2]-radius, center[2]+radius)

# 簡易アニメ（スペースで止めたいなら event を足す）
for i in range(len(frames)):
    draw_frame(i)
    plt.pause(0.03)

plt.show()
