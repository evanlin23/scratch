"""Training grid of PREREG.md (variant names for gc.py run)."""
GRID = (["raw:mcl:2", "raw:mcl:3", "raw:mcl:6"] + ["es4:mcl:%s" % i for i in (2, 3, 4, 6)]
        + ["%s:leidmod:%s" % (g, r) for g in ("raw", "es4") for r in (10, 30, 100)]
        + ["%s:leidcpm:%s" % (g, r) for g in ("raw", "es4") for r in (0.003, 0.01, 0.02)]
        + ["raw:louvain:30", "es4:louvain:30"]
        + ["raw:cc:%d" % t for t in (6, 8, 10)]
        + ["raw:agglo:%d" % t for t in (1, 2, 4)] + ["es4:agglo:1"]
        + ["raw:lpa", "es4:lpa"])
if __name__ == "__main__":
    print(" ".join(["raw:mcl:4"] + GRID))
